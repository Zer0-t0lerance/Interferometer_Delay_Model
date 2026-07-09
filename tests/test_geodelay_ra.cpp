// test_geodelay_ra.cpp
//
// Геоцентрическая задержка КОСМИЧЕСКОГО телескопа RASTRON на 13:00:00 (позиция из орбиты
// .scf через orbit_interp) против эталона example/RA_L.TXT P0 (-4.74842393132049e-1).
// Геоцентр = station1 (is_space, r=0), RASTRON = station2 (is_space, x_orbit из орбиты).
// Покажет: чисто геоцентрическая задержка или в неё входит даунлинк на движущийся пункт
// приёма (TIMEOFS ~ -0.7 с).

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

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    printf("=====================================================================\n");
    printf("  Геоцентрическая задержка RASTRON на 13:00:00 vs эталон RA_L.TXT\n");
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
    Observation o0{}; o0.mjd = mjd; o0.utc = utc;
    double ut1o; Eigen::VectorXd ei(5), dei(5), ar(8); Eigen::MatrixXd dd(3, 2), dl(3, 2);
    interp_eop(0, o0, TT, ut1o, ei, dei, ar, dd, dl, eop);
    double xp = ei(1), yp = ei(2), ut1_utc = ei(0);
    double ut1_frac = utc + ut1_utc / cnst::SECDAY;
    double ct = t_eph40(jd0, TT, ut1_frac, 0, 0, 0);
    double jd_tt = jd0 + TT, jd_ut1 = jd0 + ut1_frac;
    double cent = (jd_tt - cnst::JD2000) / cnst::JUL_CENT, ut1_sec = utc * cnst::SECDAY + ut1_utc;
    Eigen::Matrix3d Earth, Sun, Moon; get_celestial_bodies(jd0, ct, Earth, Sun, Moon);
    Eigen::Matrix3d Eg; Eigen::MatrixXd sm, mm; jpl_eph(jd0, ct, Eg, sm, mm);
    Eigen::Matrix<double, 3, 2> sun_geo = sm, moon_geo = mm;
    Eigen::Matrix3d R, dR, d2R; get_r2000_matrices(jd_tt, jd_ut1, xp, yp, R, dR, d2R);
    Eigen::VectorXd f(5), fd(5); double c2; fund_arg(jd0, ct, c2, f, fd);
    Eigen::Vector2d gast = gast_iau2006(jd_tt, jd_ut1);

    std::vector<Source> src(1); src[0].ra = 137.2920483250 * cnst::CDEGRAD; src[0].dec = 1.3598938083 * cnst::CDEGRAD;
    std::vector<Eigen::Vector3d> kv; source_vec(src, mjd + utc, kv);

    // RASTRON: положение из орбиты (J2000, м).
    Eigen::Vector3d x, v, a; orbit_interp(orbit, mjd + utc, x, v, a);
    printf("  RASTRON @13:00: |x|=%.1f км  |v|=%.4f км/с\n", x.norm() / 1e3, v.norm() / 1e3);
    SitePrep geo; geo.is_space = true;
    SitePrep ra; ra.is_space = true; ra.x_orbit = x; ra.v_orbit = v; ra.a_orbit = a;

    Observation obs{}; obs.sta1 = 0; obs.sta2 = 1;
    double tau, dtau;
    compute_delay_obs(geo, ra, kv[0], obs, mjd, utc, jd0, ct, cent, ut1_sec, f, fd, gast,
                      Earth, Sun, Moon, sun_geo, moon_geo, xp, yp, 0, 0, R, dR, d2R, tau, dtau, nullptr, false);

    const double P0 = -4.74842393132049e-1;
    printf("  наша tau (геоцентр->RASTRON) = % .12e с\n", tau);
    printf("  эталон P0 (RA_L.TXT)         = % .12e с\n", P0);
    printf("  |tau - P0| = %.3e с (%.1f м)\n", std::fabs(tau - P0), std::fabs(tau - P0) * 3e8);
    printf("  (для справки: TIMEOFS00 = -0.6965 с — даунлинк космос->пункт приёма)\n");
    printf("---------------------------------------------------------------------\n");
    bool close = std::fabs(tau - P0) < 1.0; // в пределах метров-км -> геоцентрическая
    printf("  РЕЗУЛЬТАТ: %s (|Δ|=%.3e с)\n", close ? "близко к геоцентру" : "далеко (нужен даунлинк/пункт приёма)", std::fabs(tau - P0));
    return 0; // диагностический тест
}
