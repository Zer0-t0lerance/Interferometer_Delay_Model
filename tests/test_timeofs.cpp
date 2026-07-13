// test_timeofs.cpp
//
// Проверка формулы TIMEOFS: даунлинк-задержка сигнала RASTRON -> наземный пункт приёма
// (Пущино) на 13:00:00 против эталона из cfx TIMEOFS00 = -6.965338706732580e-001 с.
// R_RASTRON(J2000) из орбиты; R_пущино из геодезии (ITRF) -> J2000 через r2000(t)
// (пункт приёма ДВИЖЕТСЯ с вращением Земли). TIMEOFS = -|R_RASTRON - R_пущино|/c.
// Проверяем также вариант с ретардацией (позиция RASTRON в момент излучения t - tau).

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif
#include "../src/functions.h"
#include <cstdio>
#include <cmath>
#include <fstream>

using namespace ariadna;

static std::string find_eph() {
    const char* c[] = {"external/dephem-master/linux_p1550p2650.440t", "../external/dephem-master/linux_p1550p2650.440t"};
    for (auto p : c) { std::ifstream f(p, std::ios::binary); if (f.good()) return p; } return c[0];
}
static double dms(double d, double m, double s) { return (d + m / 60.0 + s / 3600.0) * cnst::CDEGRAD; }

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    printf("=====================================================================\n");
    printf("  TIMEOFS: даунлинк RASTRON->Пущино на 13:00 vs эталон -0.6965339 с\n");
    printf("---------------------------------------------------------------------\n");
    try { init_ephemeris(find_eph()); } catch (const std::exception& e) { printf("SKIP: %s\n", e.what()); return 2; }

    std::vector<SpaceStation> orbit;
    if (!read_scf_orbit("example/RA140423_1200_v02.scf", orbit)) { printf("FAIL: scf\n"); return 1; }

    const int mjd = 56770; const double utc = 46800.0 / cnst::SECDAY; const double jd0 = cnst::MJD_OFFSET + mjd;
    const double rows[7][4] = {
        {56767, 0.060944, 0.436953, -0.2295574}, {56768, 0.062152, 0.437984, -0.2308630},
        {56769, 0.063753, 0.438706, -0.2322999}, {56770, 0.065273, 0.439403, -0.2338988},
        {56771, 0.066868, 0.439919, -0.2356298}, {56772, 0.068557, 0.440457, -0.2374814},
        {56773, 0.070301, 0.441033, -0.2394680}};
    std::vector<EOPData> eop(7);
    for (int i = 0; i < 7; ++i) { eop[i].mjd = rows[i][0]; eop[i].x = rows[i][1]; eop[i].y = rows[i][2]; eop[i].ut1_utc = rows[i][3]; eop[i].ut1_tai = rows[i][3] - 35.0; }

    double TAI, TT; tai_time(mjd, utc, TAI, TT);
    Observation o{}; o.mjd = mjd; o.utc = utc;
    double u1; Eigen::VectorXd ei(5), dei(5), ar(8); Eigen::MatrixXd dd(3, 2), dl(3, 2);
    interp_eop(0, o, TT, u1, ei, dei, ar, dd, dl, eop);
    double xp = ei(1), yp = ei(2), ut1_frac = utc + ei(0) / cnst::SECDAY;
    double jd_tt = jd0 + TT, jd_ut1 = jd0 + ut1_frac;
    Eigen::Matrix3d R, dR, d2R; get_r2000_matrices(jd_tt, jd_ut1, xp, yp, R, dR, d2R);

    // Пущино: geodetic -> ITRF XYZ.
    double lon = dms(37, 37, 41.8), lat = dms(54, 49, 14.2), h = 239.7;
    double a = cnst::AE, f = 1.0 / cnst::F, e2 = f * (2.0 - f);
    double N = a / std::sqrt(1.0 - e2 * std::sin(lat) * std::sin(lat));
    Eigen::Vector3d pu_itrf((N + h) * std::cos(lat) * std::cos(lon),
                            (N + h) * std::cos(lat) * std::sin(lon),
                            (N * (1.0 - e2) + h) * std::sin(lat));
    Eigen::Vector3d pu_j2000 = R * pu_itrf; // пункт приёма в J2000 (движется с Землёй)

    // RASTRON в J2000.
    Eigen::Vector3d x, v, aa; orbit_interp(orbit, mjd + utc, x, v, aa);

    double dist = (x - pu_j2000).norm();
    double timeofs = -dist / cnst::C;
    // С ретардацией: RASTRON в момент излучения t + timeofs (сигнал шёл timeofs).
    Eigen::Vector3d xr, vr, ar2; orbit_interp(orbit, mjd + utc + timeofs / cnst::SECDAY, xr, vr, ar2);
    double timeofs_ret = -((xr - pu_j2000).norm()) / cnst::C;

    const double ref = -6.965338706732580e-001;
    printf("  RASTRON |x|=%.1f км;  Пущино |R_itrf|=%.1f км\n", x.norm() / 1e3, pu_itrf.norm() / 1e3);
    printf("  простой   TIMEOFS = % .12e с  |Δэт|=%.2e с (%.1f м)\n", timeofs, std::fabs(timeofs - ref), std::fabs(timeofs - ref) * 3e8);
    printf("  ретардация TIMEOFS = % .12e с  |Δэт|=%.2e с (%.1f м)\n", timeofs_ret, std::fabs(timeofs_ret - ref), std::fabs(timeofs_ret - ref) * 3e8);
    printf("  эталон             = % .12e с\n", ref);
    printf("---------------------------------------------------------------------\n");
    return 0;
}
