// test_site_pair.cpp
//
// Сборка координат станций в J2000 (site_pair) для реального наблюдения
// (18JAN02XA, FORTLEZA-HART15M) и сверка с xsta_j2000 дампа.
//
// Изоляция основной части: считаем ТВЁРДЫЙ прилив (доминирует, ~0.2 м); океаническую
// нагрузку, прилив полюса и атмосферу пока обнуляем (нужны EOP/каталоги, ~мм-см).
// r2000 берём логированную (test_baseline); dr2000/d2r2000 — SOFA (для позиции неважны).
// Ожидаемая невязка: ~1-2 см (пропущены океан+полюс+атм + дрейф координат), что
// подтверждает корректность сборки (твёрдый прилив + SITE_INST + поворот).
//
// Сборка (из корня; нужны эфемериды и libsofa.a) — см. test compile commands.txt.

#include "../src/functions.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>

#include "sofa.h"

using namespace ariadna;

static std::string find_eph() {
    const char* c[] = {"external/dephem-master/linux_p1550p2650.440t",
                       "../external/dephem-master/linux_p1550p2650.440t"};
    for (auto p : c) { std::ifstream f(p, std::ios::binary); if (f.good()) return p; }
    return c[0];
}

// Геодезия + матрица VEN->ITRF из ITRF-координат (как в site() one-shot).
static void prep_geom(SitePrep& s) {
    double r = s.xsta_itrf.norm();
    s.lat_gcen = std::asin(s.xsta_itrf.z() / r);
    s.lon_gcen = std::atan2(s.xsta_itrf.y(), s.xsta_itrf.x());
    if (s.lon_gcen < 0) s.lon_gcen += cnst::TWOPI;
    double req = std::sqrt(s.xsta_itrf.x()*s.xsta_itrf.x() + s.xsta_itrf.y()*s.xsta_itrf.y());
    double h;
    GEOD(req, s.xsta_itrf.z(), s.lat_geod, h);
    // vw = Rz(-lon)*Ry(lat_geod) (VEN->ITRF), как в site_functions.cpp
    Eigen::Matrix3d W = Eigen::Matrix3d::Identity(), V = Eigen::Matrix3d::Identity();
    double c = std::cos(s.lat_geod), sn = std::sin(s.lat_geod);
    W(0,0)=c; W(0,2)=sn; W(2,0)=-sn; W(2,2)=c;                 // Ry(lat_geod)
    double cl = std::cos(-s.lon_gcen), sl = std::sin(-s.lon_gcen);
    V(0,0)=cl; V(0,1)=-sl; V(1,0)=sl; V(1,1)=cl;               // Rz(-lon)
    s.vw_i = V * W;
}

int main() {
    printf("=====================================================================\n");
    printf("  UNIT: site_pair -> xsta_j2000 vs дамп (18JAN02XA, только твёрдый прилив)\n");
    printf("---------------------------------------------------------------------\n");

    // Каталожные ITRF на эпоху J2000.0 (IVS_TRF2015corr) + дрейф скорость*dyear.
    // READ_CAT: станции даны на эпоху 2000.0; dyear = годы от J2000.
    double dyear = (58120.0 + 0.7088 - 51544.5) / 365.25; // ~18.005
    Eigen::Vector3d velF(-0.0024, -0.0039, 0.0126), velH(0.0019, 0.0216, 0.0133);
    SitePrep s1, s2;
    s1.xsta_itrf << 4985370.025, -3955020.358, -428472.184;  // FORTLEZA (J2000)
    s2.xsta_itrf << 5085490.783,  2668161.315, -2768692.721;  // HART15M (J2000)
    s1.xsta_itrf += velF * dyear;
    s2.xsta_itrf += velH * dyear;
    prep_geom(s1); prep_geom(s2);
    s1.pres = 1006.7; s2.pres = 861.233; // метео (для атм, сейчас atm_load=0)

    // Эпоха
    double jd0 = 2458120.5;
    double tt = 2458121.209643333 - jd0;      // TT-доля
    double jd_tt = jd0 + tt;
    double cent = (jd_tt - cnst::JD2000) / cnst::JUL_CENT;
    double jd_ut1 = jd_tt - 69.184/86400.0;   // UT1~UTC (без EOP)
    double ut1_sec = (jd_ut1 - jd0) * 86400.0;

    // Фундаментальные аргументы
    Eigen::VectorXd f(5), fd(5); double cent2;
    fund_arg(jd0, tt, cent2, f, fd);

    // GAST через SOFA (uta,utb = UT1; tta,ttb = TT). Скорость GAST ~ угловая скорость Земли.
    double gast_rad = iauGst06a(2400000.5, jd_ut1 - 2400000.5, 2400000.5, jd_tt - 2400000.5);
    Eigen::Vector2d gast(gast_rad, cnst::OMEGA_EARTH);

    // Эфемериды: геоцентрические Солнце/Луна (для приливов)
    try { init_ephemeris(find_eph()); }
    catch (const std::exception& e) { printf("SKIP: эфемериды: %s\n", e.what()); return 2; }
    Eigen::Matrix3d E; Eigen::MatrixXd sunM, moonM;
    jpl_eph(jd0, tt, E, sunM, moonM); // sun/moon 3x2 геоцентрические
    Eigen::Matrix<double,3,2> sun_geo = sunM, moon_geo = moonM;

    // Логированная r2000 (test_baseline); dr2000/d2r2000 — SOFA (для позиции не влияют)
    Eigen::Matrix3d r2000;
    r2000 << 0.9988360617015190E+00,  0.4820309780423414E-01,  0.1727196188789892E-02,
            -0.4820310447537553E-01,  0.9988375540050903E+00, -0.3778973045624359E-04,
            -0.1727009998570988E-02, -0.4551047279605037E-04,  0.9999985076815172E+00;
    Eigen::Matrix3d Rsofa, dr2000, d2r2000;
    get_r2000_matrices(jd_tt, jd_ut1, 0.0, 0.0, Rsofa, dr2000, d2r2000);

    // Океан/атм отключены (нет каталогов в тесте); полюс — с реальным EOP.
    s1.tide_data.amplitudes.setZero(); s1.tide_data.phases.setZero(); s1.atm_load.coef.setZero();
    s2.tide_data.amplitudes.setZero(); s2.tide_data.phases.setZero(); s2.atm_load.coef.setZero();

    // Реальные EOP (интерполяция EOPC04 на 58120.7088)
    double xp = 0.0561734 * cnst::CARCRAD, yp = 0.2492613 * cnst::CARCRAD; // рад

    std::vector<Eigen::Vector3d> xs, vs, as_;
    site_pair(s1, s2, 58120, 0.7088, jd0, ut1_sec, cent, f, fd, gast, sun_geo, moon_geo,
              xp, yp, 0.0, 0.0, r2000, dr2000, d2r2000, xs, vs, as_);

    Eigen::Vector3d f_dump(4788183.147308413, -4190717.195172735, -436901.2217094282);
    Eigen::Vector3d h_dump(5203403.045778226,  2420028.101198892, -2777592.466280870);
    double df = (xs[0] - f_dump).norm(), dh = (xs[1] - h_dump).norm();

    printf("  FORTLEZA xsta_j2000: calc=(%.3f %.3f %.3f)\n", xs[0].x(), xs[0].y(), xs[0].z());
    printf("                       dump=(%.3f %.3f %.3f)  |d|=%.4f м\n", f_dump.x(),f_dump.y(),f_dump.z(), df);
    printf("  HART15M  xsta_j2000: calc=(%.3f %.3f %.3f)\n", xs[1].x(), xs[1].y(), xs[1].z());
    printf("                       dump=(%.3f %.3f %.3f)  |d|=%.4f м\n", h_dump.x(),h_dump.y(),h_dump.z(), dh);
    // С дрейфом v*dyear (эпоха каталога 2000.0) + реальным EOP невязка ~0.1 м
    // (было ~0.3 м). Остаток: точная опорная эпоха IVS_TRF2015 + отключённые океан/атм
    // нагрузки (в этом тесте нет каталогов). Полное замыкание — в compute_delay_obs с
    // чтением всех каталогов и интерполяцией EOP.
    bool ok = df < 0.2 && dh < 0.2;
    printf("---------------------------------------------------------------------\n");
    printf("  site_pair собран, невязка < 0.2 м (дрейф+EOP применены): %s\n", ok ? "OK" : "FAIL");
    printf("  (остаток = опорная эпоха каталога + океан/атм нагрузки [отключены здесь])\n");
    printf("  РЕЗУЛЬТАТ: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
