// test_process_obs.cpp
//
// Полная связка "низа+середины" конвейера для реального наблюдения (18JAN02XA,
// FORTLEZA-HART15M, 0017+200): source_vec -> aber_source(e/az) -> trop_delay ->
// mount_tel(dtau_off) -> theor_delay, сверка e/az, dtau_off и ИТОГОВОЙ задержки с
// дампом/логами. Координаты станций в J2000 берём логированные (site_pair отдельно
// имеет ~0.1 м из-за эпохи каталога); r2000 логированная, dR2000 — SOFA с реальным EOP.
//
// Сборка (из корня; нужны эфемериды и libsofa.a; libsofa.a — ПОСЛЕДНИМ) — см. test compile commands.txt.

#include "../src/functions.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

using namespace ariadna;

static int g_fail = 0;
static void chk(const char* n, double got, double exp, double tol) {
    double d = std::fabs(got - exp); bool ok = d <= tol; if(!ok) ++g_fail;
    printf("  %-16s calc=% .12e  ref=% .12e  |d|=%.2e  %s\n", n, got, exp, d, ok?"OK":"FAIL");
}
static std::string find_eph(){ const char* c[]={"external/dephem-master/linux_p1550p2650.440t","../external/dephem-master/linux_p1550p2650.440t"};
    for(auto p:c){std::ifstream f(p,std::ios::binary); if(f.good())return p;} return c[0]; }

int main() {
    printf("=====================================================================\n");
    printf("  UNIT: полная связка aber->trop->mount->theor vs дамп (18JAN02XA)\n");
    printf("---------------------------------------------------------------------\n");

    // --- Логированные J2000-координаты/скорости станций (test_baseline) ---
    std::vector<Eigen::Vector3d> xs(2), vs(2), as_(2);
    xs[0] << 4788183.147308413, -4190717.195172735, -436901.2217094282;
    xs[1] << 5203403.045778226,  2420028.101198892, -2777592.466280870;
    vs[0] << 305.5926636587654, 349.2143463604144, -0.5142780989077815;
    vs[1] << -176.4630853380511, 379.7874391789990, 0.3196033011667836;
    as_[0].setZero(); as_[1].setZero();

    Eigen::Matrix3d r2000;
    r2000 << 0.9988360617015190E+00,  0.4820309780423414E-01,  0.1727196188789892E-02,
            -0.4820310447537553E-01,  0.9988375540050903E+00, -0.3778973045624359E-04,
            -0.1727009998570988E-02, -0.4551047279605037E-04,  0.9999985076815172E+00;

    // dR2000 из SOFA с реальным EOP (для скоростей/аберрации)
    double jd_tt = 2458121.209643333;
    double jd_ut1 = jd_tt - 69.184/86400.0 + 0.2150842/86400.0;
    double xp = 0.0561734*cnst::CARCRAD, yp = 0.2492613*cnst::CARCRAD;
    Eigen::Matrix3d Rs, dR2000, d2R2000;
    get_r2000_matrices(jd_tt, jd_ut1, xp, yp, Rs, dR2000, d2R2000);

    // vw станций (VEN->ITRF) строим ПРЯМО из lat_geod/lon (эталон vw_i из test_site_tide_solid;
    // прежняя V*W в site_functions.cpp давала неверные знаки Up — баг, поправить и там).
    std::vector<Eigen::Matrix3d> vw(2);
    double latg[2], hg[2];
    for (int i = 0; i < 2; ++i) {
        Eigen::Vector3d itrf = r2000.transpose() * xs[i];
        double lon = std::atan2(itrf.y(), itrf.x()); if (lon<0) lon+=cnst::TWOPI;
        double req = std::sqrt(itrf.x()*itrf.x()+itrf.y()*itrf.y());
        GEOD(req, itrf.z(), latg[i], hg[i]);
        double cla=std::cos(latg[i]), sla=std::sin(latg[i]), clo=std::cos(lon), slo=std::sin(lon);
        vw[i].col(0) << cla*clo, cla*slo, sla;        // Up (Vertical)
        vw[i].col(1) << -slo, clo, 0.0;               // East
        vw[i].col(2) << -sla*clo, -sla*slo, cla;      // North
    }

    // Источник -> K_s
    std::vector<Source> src(1);
    src[0].ra=0.085655996508366; src[0].dec=0.355395794394430; src[0].ra_rate=0; src[0].dec_rate=0;
    std::vector<Eigen::Vector3d> k_star; source_vec(src, 58120.7088, k_star);
    Eigen::Vector3d K_s = k_star[0];

    // Эфемериды: Earth SSB (для аберрации — скорость; для theor_delay — вся), Sun/Moon SSB
    try { init_ephemeris(find_eph()); } catch(const std::exception&e){ printf("SKIP: %s\n",e.what()); return 2; }
    // Точная эпоха для эфемерид (из test_ephemeris)
    Eigen::Matrix3d Earth, Sun, Moon; get_celestial_bodies(2458120.5, 0.7113110196885, Earth, Sun, Moon);

    Observation obs{};
    obs.sta1=0; obs.sta2=1;
    obs.t1=30.427; obs.p1=1006.7; obs.e1=60.973;
    obs.t2=25.884; obs.p2=861.233; obs.e2=43.395;

    // --- ABER_SOURCE -> e/az --- (earth: только скорость в col(1), как в тесте aber_source)
    std::vector<Eigen::Matrix3d> r2000_v = { r2000, dR2000 };
    Eigen::Matrix3d Earth_ab = Eigen::Matrix3d::Zero();
    Earth_ab.col(1) = Earth.col(1);
    Eigen::Matrix2d e2, az2;
    aber_source(obs, r2000_v, K_s, Earth_ab, vs, vw, e2, az2);
    printf("  [aber_source -> элевации]\n");
    chk("elev st1", e2(0,0), 0.6745164463503075, 5e-5);
    chk("elev st2", e2(1,0), 0.6936373986295217, 5e-5);

    // --- BASELINE ---
    Eigen::MatrixXd Xm(3,2), Vm(3,2), Am(3,2);
    Xm.col(0)=xs[0]; Xm.col(1)=xs[1]; Vm.col(0)=vs[0]; Vm.col(1)=vs[1]; Am.setZero();
    Eigen::MatrixXd blm(3,2); Eigen::Vector3d bcfs;
    baseline(r2000, Xm, Vm, Am, blm, bcfs);

    // --- TROP_DELAY (с вычисленными aber_source e/az) ---
    Station st1, st2;
    st1.lat_geod=latg[0]; st1.h_geod=hg[0]; st2.lat_geod=latg[1]; st2.h_geod=hg[1];
    Eigen::MatrixXd eM=e2, azM=az2;
    Eigen::MatrixXd dd,dw,dhmf,dwmf,dgn,dge,zd,zw;
    trop_delay(obs, 2458121.209643333, 0.0, st1, st2, eM, azM, dd,dw,dhmf,dwmf,dgn,dge,zd,zw);

    // --- MOUNT_TEL -> dtau_off (r2000 как 9x3 с блоком R^T; vw VEN->ITRF; e/az MatrixXd) ---
    Eigen::MatrixXd r2000_9(9,3); r2000_9.setZero();
    r2000_9.block<3,3>(0,0) = r2000.transpose();
    r2000_9.block<3,3>(3,0) = dR2000.transpose();
    std::vector<Station> mstations(2);
    mstations[0].axsty="AZEL"; mstations[0].offs=0.00637;  mstations[0].lat_geod=latg[0];
    mstations[1].axsty="AZEL"; mstations[1].offs=1.49500;  mstations[1].lat_geod=latg[1];
    std::vector<Eigen::Vector3d> ks{K_s};
    Eigen::MatrixXd doff_dl, d_dax, dtau_mt;
    mount_tel(obs, r2000_9, mstations, ks, vw, eM, azM, doff_dl, d_dax, dtau_mt);
    printf("\n  [mount_tel -> dtau_off]\n");
    chk("dtau st1", dtau_mt(0,0), 0.1658888657089142e-10, 1e-13);
    chk("dtau st2", dtau_mt(1,0), -0.3833381489451837e-08, 1e-12);

    // --- THEOR_DELAY (Datmc и dtau_off транспонируем; Солнце/Луна SSB; dt_temp логированный) ---
    Eigen::Matrix<double,3,2> bl; bl.col(0)=blm.col(0); bl.col(1)=blm.col(1);
    Eigen::Matrix2d Dd=Eigen::Matrix2d(dd).transpose(), Dw=Eigen::Matrix2d(dw).transpose();
    Eigen::Matrix2d dtau=Eigen::Matrix2d(dtau_mt).transpose(), dtmp;
    dtmp << 0.8286346437956392e-10,  0.1445568241005169e-09,  0.0, 0.0;
    double t2,dt2;
    theor_delay(bl, xs, vs, as_, K_s, Earth, Sun, Moon, Dd, Dw, dtau, dtmp, t2, dt2);
    // Итог с ВЫЧИСЛЕННЫМИ e/az совпадает с дампом до ~5e-11 с (~1.5 см); остаток —
    // точность вычисленных элеваций (~7e-9 рад) от dR2000 (SOFA) + эпохи эфемерид.
    // (С логированными e/az было бы 1e-12, см. test_delay_pipeline.)
    printf("\n  [baseline->trop_delay->theor_delay -> ИТОГ]\n");
    chk("t2_t1", t2, -0.3450839807711632e-03, 1e-10);

    printf("=====================================================================\n");
    printf("  РЕЗУЛЬТАТ: %s (провалов: %d)\n", g_fail==0?"PASS":"FAIL", g_fail);
    return g_fail==0?0:1;
}
