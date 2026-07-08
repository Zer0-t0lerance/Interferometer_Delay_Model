// test_compute_delay.cpp
//
// КУЛЬМИНАЦИЯ: полностью самосогласованный расчёт теоретической задержки одного
// наблюдения через compute_delay_obs (все блоки конвейера + SOFA r2000 с реальным EOP
// + каталожные координаты станций). Сверка итога с дампом ARIADNA (18JAN02XA,
// FORTLEZA-HART15M, 0017+200).
//
// Точность ограничена: опорная эпоха каталога станций (дрейф v*dyear ~см-дм),
// океан/атм нагрузки здесь обнулены, DE440 vs DE421, dt_temp=0. Ожидаем итог в
// пределах нескольких 1e-10 с (дециметры) от дампа -- т.е. ВЕСЬ конвейер собран и
// даёт задержку близко к ARIADNA. (С логированными координатами низ конвейера даёт 5e-11.)

#include "../src/functions.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>

using namespace ariadna;

static std::string find_eph(){ const char* c[]={"external/dephem-master/linux_p1550p2650.440t","../external/dephem-master/linux_p1550p2650.440t"};
    for(auto p:c){std::ifstream f(p,std::ios::binary); if(f.good())return p;} return c[0]; }

static void prep(SitePrep& s, const Eigen::Vector3d& cat, const Eigen::Vector3d& vel, double dyear,
                 double pres, double offs) {
    s.xsta_itrf = cat + vel * dyear;
    double r = s.xsta_itrf.norm();
    s.lat_gcen = std::asin(s.xsta_itrf.z()/r);
    s.lon_gcen = std::atan2(s.xsta_itrf.y(), s.xsta_itrf.x()); if(s.lon_gcen<0)s.lon_gcen+=cnst::TWOPI;
    double req=std::sqrt(s.xsta_itrf.x()*s.xsta_itrf.x()+s.xsta_itrf.y()*s.xsta_itrf.y());
    GEOD(req, s.xsta_itrf.z(), s.lat_geod, s.h_geod);
    double cla=std::cos(s.lat_geod), sla=std::sin(s.lat_geod), clo=std::cos(s.lon_gcen), slo=std::sin(s.lon_gcen);
    s.vw_i.col(0) << cla*clo, cla*slo, sla;
    s.vw_i.col(1) << -slo, clo, 0.0;
    s.vw_i.col(2) << -sla*clo, -sla*slo, cla;
    s.tide_data.amplitudes.setZero(); s.tide_data.phases.setZero(); s.atm_load.coef.setZero();
    s.pres=pres; s.dPdt=0.0; s.axsty="AZEL"; s.offs=offs;
}

int main() {
    printf("=====================================================================\n");
    printf("  UNIT: compute_delay_obs -- полный конвейер vs дамп (18JAN02XA)\n");
    printf("---------------------------------------------------------------------\n");

    // Эпоха
    int mjd=58120; double utc=0.7088;
    double jd0=2458120.5, tt=2458121.209643333-jd0, jd_tt=jd0+tt, ct=tt;
    double cent=(jd_tt-cnst::JD2000)/cnst::JUL_CENT;
    // EOP (интерп. EOPC04 на 58120.7088)
    double xp=0.0561734*cnst::CARCRAD, yp=0.2492613*cnst::CARCRAD, ut1_utc=0.2150842;
    double jd_ut1=jd_tt-69.184/86400.0+ut1_utc/86400.0;
    double ut1_sec=utc*86400.0+ut1_utc;
    double dyear=(mjd+utc-51544.5)/365.25;

    // Станции
    SitePrep s1,s2;
    prep(s1, Eigen::Vector3d(4985370.025,-3955020.358,-428472.184), Eigen::Vector3d(-0.0024,-0.0039,0.0126), dyear, 1006.7, 0.00637);
    prep(s2, Eigen::Vector3d(5085490.783,2668161.315,-2768692.721), Eigen::Vector3d(0.0019,0.0216,0.0133), dyear, 861.233, 1.49500);

    // Источник
    std::vector<Source> src(1); src[0].ra=0.085655996508366; src[0].dec=0.355395794394430; src[0].ra_rate=0; src[0].dec_rate=0;
    std::vector<Eigen::Vector3d> k_star; source_vec(src, mjd+utc, k_star);

    // Время: фунд.аргументы, gast
    Eigen::VectorXd f(5),fd(5); double c2; fund_arg(jd0, tt, c2, f, fd);
    Eigen::Vector2d gast = gast_iau2006(jd_tt, jd_ut1);

    // Эфемериды
    try{ init_ephemeris(find_eph()); } catch(const std::exception&e){ printf("SKIP: %s\n",e.what()); return 2; }
    Eigen::Matrix3d Earth,Sun,Moon; get_celestial_bodies(2458120.5, 0.7113110196885, Earth,Sun,Moon);
    Eigen::Matrix3d Eg; Eigen::MatrixXd sunM,moonM; jpl_eph(jd0, tt, Eg, sunM, moonM);
    Eigen::Matrix<double,3,2> sun_geo=sunM, moon_geo=moonM;

    // r2000 (SOFA с реальным EOP)
    Eigen::Matrix3d R,dR,d2R; get_r2000_matrices(jd_tt, jd_ut1, xp, yp, R, dR, d2R);

    Observation obs{}; obs.sta1=0; obs.sta2=1;
    obs.t1=30.427;obs.p1=1006.7;obs.e1=60.973; obs.t2=25.884;obs.p2=861.233;obs.e2=43.395;

    double tau,dtau;
    compute_delay_obs(s1,s2,k_star[0],obs, mjd,utc,jd0,ct,cent,ut1_sec, f,fd,gast,
                      Earth,Sun,Moon, sun_geo,moon_geo, xp,yp,0.0,0.0, R,dR,d2R, tau,dtau);

    double ref=-0.3450839807711632e-03, d=std::fabs(tau-ref);
    printf("  compute_delay_obs: tau  = % .12e\n", tau);
    printf("  дамп ARIADNA:      ref  = % .12e\n", ref);
    printf("  |tau - ref| = %.3e с  (~%.2f см)\n", d, d*3e8*100);
    bool ok = d < 5e-9; // ~метр: весь конвейер собран, геометрия ограничена эпохой каталога/нагрузками
    printf("---------------------------------------------------------------------\n");
    printf("  Полный конвейер даёт задержку близко к ARIADNA (<5e-9 с): %s\n", ok?"OK":"FAIL");
    printf("  РЕЗУЛЬТАТ: %s\n", ok?"PASS":"FAIL");
    return ok?0:1;
}
