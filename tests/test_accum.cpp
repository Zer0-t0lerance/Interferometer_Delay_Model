// test_accum.cpp
//
// Проверка НАКОПЛЕНИЯ ошибки: считаем геоцентрическую вакуумную задержку BADARY
// НЕЗАВИСИМО (свежий расчёт с нуля) в начале (13:00), середине (13:30) и конце (13:59)
// сеанса и сравниваем с эталоном BADARY.TXT. Если ошибка растёт ДАЖЕ в независимых
// расчётах (не накапливается от прогона) — значит дело в каталогах/эпохе, а не в
// накоплении/утечке. Если бы росла только при последовательном прогоне — был бы баг.

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

// Свежий (независимый) расчёт вакуумной геоцентрической задержки BADARY на момент.
static double fresh_badary_delay(int mjd, double utc, const std::vector<EOPData>& eop) {
    double jd0 = cnst::MJD_OFFSET + mjd;
    double TAI, TT; tai_time(mjd, utc, TAI, TT);
    Observation o{}; o.mjd = mjd; o.utc = utc;
    double u1; Eigen::VectorXd ei(5), dei(5), ar(8); Eigen::MatrixXd dd(3, 2), dl(3, 2);
    interp_eop(0, o, TT, u1, ei, dei, ar, dd, dl, eop);
    double xp = ei(1), yp = ei(2), ut1_utc = ei(0), ut1_frac = utc + ut1_utc / cnst::SECDAY;
    double ct = t_eph40(jd0, TT, ut1_frac, 0, 0, 0), jd_tt = jd0 + TT, jd_ut1 = jd0 + ut1_frac;
    double cent = (jd_tt - cnst::JD2000) / cnst::JUL_CENT, ut1_sec = utc * cnst::SECDAY + ut1_utc;
    Eigen::Matrix3d Earth, Sun, Moon; get_celestial_bodies(jd0, ct, Earth, Sun, Moon);
    Eigen::Matrix3d Eg; Eigen::MatrixXd sm, mm; jpl_eph(jd0, ct, Eg, sm, mm);
    Eigen::Matrix<double, 3, 2> sun_geo = sm, moon_geo = mm;
    Eigen::Matrix3d R, dR, d2R; get_r2000_matrices(jd_tt, jd_ut1, xp, yp, R, dR, d2R);
    Eigen::VectorXd f(5), fd(5); double c2; fund_arg(jd0, ct, c2, f, fd);
    Eigen::Vector2d gast = gast_iau2006(jd_tt, jd_ut1);
    std::vector<Source> src(1); src[0].ra = 137.2920483250 * cnst::CDEGRAD; src[0].dec = 1.3598938083 * cnst::CDEGRAD;
    std::vector<Eigen::Vector3d> kv; source_vec(src, mjd + utc, kv);
    SitePrep geo; geo.is_space = true;
    SitePrep bad;
    double dyear = (mjd + utc - 54466.0) / cnst::DAYS_PER_YEAR;
    bad.xsta_itrf = Eigen::Vector3d(-838200.93240, 3865751.56640, 4987670.90800) + Eigen::Vector3d(-0.026960, -0.001420, -0.003500) * dyear;
    double r = bad.xsta_itrf.norm();
    bad.lat_gcen = std::asin(bad.xsta_itrf.z() / r);
    bad.lon_gcen = std::atan2(bad.xsta_itrf.y(), bad.xsta_itrf.x()); if (bad.lon_gcen < 0) bad.lon_gcen += cnst::TWOPI;
    double req = std::sqrt(bad.xsta_itrf.x() * bad.xsta_itrf.x() + bad.xsta_itrf.y() * bad.xsta_itrf.y());
    GEOD(req, bad.xsta_itrf.z(), bad.lat_geod, bad.h_geod);
    double cla = std::cos(bad.lat_geod), sla = std::sin(bad.lat_geod), clo = std::cos(bad.lon_gcen), slo = std::sin(bad.lon_gcen);
    bad.vw_i.col(0) << cla * clo, cla * slo, sla; bad.vw_i.col(1) << -slo, clo, 0.0; bad.vw_i.col(2) << -sla * clo, -sla * slo, cla;
    bad.tide_data.amplitudes.setZero(); bad.tide_data.phases.setZero(); bad.atm_load.coef.setZero();
    bad.axsty = "AZEL"; bad.offs = 0.004; bad.pres = 1013.25;
    Observation obs{}; obs.sta1 = 0; obs.sta2 = 1; obs.p2 = 1013.25; obs.e2 = 50.0;
    double tau, dtau;
    compute_delay_obs(geo, bad, kv[0], obs, mjd, utc, jd0, ct, cent, ut1_sec, f, fd, gast,
                      Earth, Sun, Moon, sun_geo, moon_geo, xp, yp, 0, 0, R, dR, d2R, tau, dtau, nullptr, false);
    return tau;
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    printf("=====================================================================\n");
    printf("  Накопление? Свежий (независимый) расчёт BADARY: начало/середина/конец\n");
    printf("---------------------------------------------------------------------\n");
    try { init_ephemeris(find_eph()); } catch (const std::exception& e) { printf("SKIP: %s\n", e.what()); return 2; }

    const double rows[7][4] = {
        {56767, 0.060944, 0.436953, -0.2295574}, {56768, 0.062152, 0.437984, -0.2308630},
        {56769, 0.063753, 0.438706, -0.2322999}, {56770, 0.065273, 0.439403, -0.2338988},
        {56771, 0.066868, 0.439919, -0.2356298}, {56772, 0.068557, 0.440457, -0.2374814},
        {56773, 0.070301, 0.441033, -0.2394680}};
    std::vector<EOPData> eop(7);
    for (int i = 0; i < 7; ++i) { eop[i].mjd = rows[i][0]; eop[i].x = rows[i][1]; eop[i].y = rows[i][2]; eop[i].ut1_utc = rows[i][3]; eop[i].ut1_tai = rows[i][3] - 35.0; }

    struct E { const char* lbl; double utc_sec; double ref; };
    E ep[] = {
        {"начало  13:00", 46800.0, -1.33142861057677e-2},
        {"середина 13:30", 48600.0, -1.28661815051242e-2},
        {"конец   13:59", 50340.0, -1.22286584492648e-2},
    };
    for (auto& e : ep) {
        double tau = fresh_badary_delay(56770, e.utc_sec / cnst::SECDAY, eop);
        double d = std::fabs(tau - e.ref);
        printf("  %s: наш=% .12e эталон=% .12e |Δ|=%.2e с (%.2f м)\n", e.lbl, tau, e.ref, d, d * 3e8);
    }
    printf("---------------------------------------------------------------------\n");
    printf("  Если |Δ| растёт от начала к концу ДАЖЕ в независимых расчётах — это\n");
    printf("  каталоги/эпоха (эфемериды/EOP), а НЕ накопление в реализации.\n");
    return 0;
}
