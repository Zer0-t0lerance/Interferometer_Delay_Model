// test_geodelay.cpp
//
// КРИТИЧНАЯ ПРОВЕРКА подхода «полиномы для коррелятора»: геоцентрическая задержка
// наземной станции (BADARY) на начало скана (2014-04-23 13:00:00) против эталонного
// P0 из example/BADARY.TXT (-1.33142861057677e-2). Геоцентр моделируем как station1
// (is_space в начале координат: положение/скорость/ускорение = 0, поправки минуются),
// BADARY — station2. Задержка = compute_delay_obs(геоцентр, BADARY).
//
// Проверяем, входит ли в эталонный P0 дефолтная CLOCK DELAY (212.437e-6): печатаем
// |tau - P0| и |tau + clock - P0|.

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
    const char* c[] = {"external/dephem-master/linux_p1550p2650.440t",
                       "../external/dephem-master/linux_p1550p2650.440t"};
    for (auto p : c) { std::ifstream f(p, std::ios::binary); if (f.good()) return p; }
    return c[0];
}

// Подготовка SitePrep наземной станции: xsta = xyz + vel*dyear, геодезия, vw.
static void prep_ground(SitePrep& s, const Eigen::Vector3d& cat, const Eigen::Vector3d& vel,
                        double dyear, double offs, const std::string& mount) {
    s.xsta_itrf = cat + vel * dyear;
    double r = s.xsta_itrf.norm();
    s.lat_gcen = std::asin(s.xsta_itrf.z() / r);
    s.lon_gcen = std::atan2(s.xsta_itrf.y(), s.xsta_itrf.x()); if (s.lon_gcen < 0) s.lon_gcen += cnst::TWOPI;
    double req = std::sqrt(s.xsta_itrf.x() * s.xsta_itrf.x() + s.xsta_itrf.y() * s.xsta_itrf.y());
    GEOD(req, s.xsta_itrf.z(), s.lat_geod, s.h_geod);
    double cla = std::cos(s.lat_geod), sla = std::sin(s.lat_geod), clo = std::cos(s.lon_gcen), slo = std::sin(s.lon_gcen);
    s.vw_i.col(0) << cla * clo, cla * slo, sla;
    s.vw_i.col(1) << -slo, clo, 0.0;
    s.vw_i.col(2) << -sla * clo, -sla * slo, cla;
    s.tide_data.amplitudes.setZero(); s.tide_data.phases.setZero(); s.atm_load.coef.setZero();
    s.axsty = mount; s.offs = offs;
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    printf("=====================================================================\n");
    printf("  Геоцентрическая задержка BADARY на 13:00:00 vs эталон BADARY.TXT\n");
    printf("---------------------------------------------------------------------\n");

    try { init_ephemeris(find_eph()); }
    catch (const std::exception& e) { printf("SKIP: эфемериды: %s\n", e.what()); return 2; }

    // Эпоха: MJD 56770, 13:00:00 UTC.
    const int mjd = 56770;
    const double utc = 46800.0 / cnst::SECDAY; // 13h
    const double jd0 = cnst::MJD_OFFSET + mjd;

    // EOP-узлы EOPC04 (56767..56773): x", y", UT1-UTC[s]; leap 2014 = 35.
    const double rows[7][4] = {
        {56767, 0.060944, 0.436953, -0.2295574}, {56768, 0.062152, 0.437984, -0.2308630},
        {56769, 0.063753, 0.438706, -0.2322999}, {56770, 0.065273, 0.439403, -0.2338988},
        {56771, 0.066868, 0.439919, -0.2356298}, {56772, 0.068557, 0.440457, -0.2374814},
        {56773, 0.070301, 0.441033, -0.2394680}};
    std::vector<EOPData> eop(7);
    for (int i = 0; i < 7; ++i) {
        eop[i].mjd = rows[i][0]; eop[i].x = rows[i][1]; eop[i].y = rows[i][2];
        eop[i].ut1_utc = rows[i][3]; eop[i].ut1_tai = rows[i][3] - 35.0;
    }

    // Время.
    double TAI, TT; tai_time(mjd, utc, TAI, TT);
    Observation obs0{}; obs0.mjd = mjd; obs0.utc = utc;
    double ut1_out; Eigen::VectorXd eop_int(5), deop_int(5), arg_oc(8); Eigen::MatrixXd dd(3, 2), dl(3, 2);
    interp_eop(0, obs0, TT, ut1_out, eop_int, deop_int, arg_oc, dd, dl, eop);
    const double xp = eop_int(1), yp = eop_int(2), ut1_utc = eop_int(0);
    const double ut1_frac = utc + ut1_utc / cnst::SECDAY;
    const double ct = t_eph40(jd0, TT, ut1_frac, 0, 0, 0);
    const double jd_tt = jd0 + TT, jd_ut1 = jd0 + ut1_frac;
    const double cent = (jd_tt - cnst::JD2000) / cnst::JUL_CENT;
    const double ut1_sec = utc * cnst::SECDAY + ut1_utc;

    // Эфемериды, ориентация.
    Eigen::Matrix3d Earth, Sun, Moon; get_celestial_bodies(jd0, ct, Earth, Sun, Moon);
    Eigen::Matrix3d Eg; Eigen::MatrixXd sunM, moonM; jpl_eph(jd0, ct, Eg, sunM, moonM);
    Eigen::Matrix<double, 3, 2> sun_geo = sunM, moon_geo = moonM;
    Eigen::Matrix3d R, dR, d2R; get_r2000_matrices(jd_tt, jd_ut1, xp, yp, R, dR, d2R);
    Eigen::VectorXd f(5), fd(5); double c2; fund_arg(jd0, ct, c2, f, fd);
    Eigen::Vector2d gast = gast_iau2006(jd_tt, jd_ut1);

    // Источник 0906+015 (RA/Dec из .cfx, в рад).
    std::vector<Source> src(1); src[0].ra = 137.2920483250 * cnst::CDEGRAD; src[0].dec = 1.3598938083 * cnst::CDEGRAD;
    std::vector<Eigen::Vector3d> k_star; source_vec(src, mjd + utc, k_star);

    // Станции: 0 = геоцентр (is_space в нуле), 1 = BADARY.
    SitePrep geo; geo.is_space = true; // x/v/a_orbit = 0 по умолчанию
    SitePrep bad;
    double dyear = (mjd + utc - 54466.0) / cnst::DAYS_PER_YEAR; // эпоха координат BADARY = 54466
    prep_ground(bad, Eigen::Vector3d(-838200.93240, 3865751.56640, 4987670.90800),
                Eigen::Vector3d(-0.026960, -0.001420, -0.003500), dyear, 0.004, "AZEL");
    bad.pres = 1013.25; bad.tC = 0.0; // стандартная атмосфера (метео в .cfx нет)

    Observation obs{}; obs.sta1 = 0; obs.sta2 = 1;
    obs.p1 = 1013.25; obs.t1 = 0.0; obs.e1 = 50.0;
    obs.p2 = 1013.25; obs.t2 = 0.0; obs.e2 = 50.0;

    double tau, dtau;
    // Вакуумная геометрия (with_tropo=false) — как полином коррелятора.
    compute_delay_obs(geo, bad, k_star[0], obs, mjd, utc, jd0, ct, cent, ut1_sec, f, fd, gast,
                      Earth, Sun, Moon, sun_geo, moon_geo, xp, yp, 0, 0, R, dR, d2R, tau, dtau,
                      nullptr, false);

    const double P0 = -1.33142861057677e-2;
    const double clock = 212.437e-6;
    printf("  наша tau (геоцентр->BADARY) = % .12e с\n", tau);
    printf("  эталон P0                   = % .12e с\n", P0);
    printf("  |tau        - P0| = %.3e с  (%.3f мкс, %.1f м)\n", std::fabs(tau - P0), std::fabs(tau - P0) * 1e6, std::fabs(tau - P0) * 3e8);
    printf("  |tau+clock  - P0| = %.3e с  (%.3f мкс, %.1f м)\n", std::fabs(tau + clock - P0), std::fabs(tau + clock - P0) * 1e6, std::fabs(tau + clock - P0) * 3e8);
    printf("  |tau-clock  - P0| = %.3e с  (%.3f мкс)\n", std::fabs(tau - clock - P0), std::fabs(tau - clock - P0) * 1e6);
    printf("---------------------------------------------------------------------\n");
    // Успех, если какая-то комбинация попадает в разумную близость (< 1 мкс = 300 м).
    double best = std::min({std::fabs(tau - P0), std::fabs(tau + clock - P0), std::fabs(tau - clock - P0)});
    bool ok = best < 1e-6;
    printf("  РЕЗУЛЬТАТ: %s (лучшее совпадение %.3e с)\n", ok ? "PASS" : "разошлось", best);
    return ok ? 0 : 1;
}
