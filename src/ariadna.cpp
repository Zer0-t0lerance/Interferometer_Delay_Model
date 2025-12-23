#include "functions.h"
#include <fstream>
#include <iostream>

namespace ariadna {

void process_ariadna(const std::vector<Station>& stations, const std::vector<Source>& sources,
                     const std::vector<Observation>& observations,
                     const std::vector<SpaceStation>& space_stations,
                     const std::vector<OrbitData>& orbit_data,
                     int n_segm, int k_ch_c, int k_ch_z, double delta_sec,
                     const std::string& output_path,
                     const std::vector<EOPData>& eop_data,
                     double mjd_mean, double utc_mean, double t_mean)
{
    // Инициализация параметров
    int n_stations = stations.size();
    // int n_sources = sources.size(); // Используется размер входного вектора sources
    int nobs = observations.size();
    int i_choice = 2; // IAU 2000 Theory of Nutation
    double dyear = (mjd_mean + utc_mean - (cnst::JD2000 - 2400000.5)) / 365.25;

    // Выделение массивов
    std::vector<Eigen::Vector3d> site_xyz(n_stations), site_vel(n_stations);
    std::vector<double> lat_geod(n_stations), h_geod(n_stations), lat_gcen(n_stations), lon_gcen(n_stations), sph_rad(n_stations), u_site(n_stations), v_site(n_stations);
    std::vector<Eigen::Matrix3d> vw(n_stations);
    std::vector<Eigen::Vector3d> site_meteo(n_stations);
    
    // Переменные для главного цикла
    int j1, j2, j3; // Индексы станций и источника
    double mjd, utc, tai, tt, ut1, ct, dtaidct;
    double gast = 0.0;
    
    Eigen::VectorXd f(5), fd(5); // Фундаментальные аргументы
    Eigen::Vector3d earth, sun, moon; // Эфемериды
    Eigen::Matrix3d r2000[3]; // Матрица поворота J2000.0 (3x3x3)
    Eigen::VectorXd eop_int(7), deop_int(7); // Интерполированные EOP и их производные
    Eigen::MatrixXd deop_diu(3, 2), deop_lib(3, 2); // Приливные поправки EOP
    Eigen::VectorXd arg_oc_tide(8); // Аргументы океанских приливов
    
    // Векторы коррекций для каждой станции (3x2: pos/vel для 2 станций)
    Eigen::MatrixXd dxtide(3, 2), dvtide(3, 2); 
    Eigen::MatrixXd dx_octide(3, 2), dv_octide(3, 2);
    Eigen::MatrixXd dx_poltide(3, 2), dv_poltide(3, 2);
    Eigen::MatrixXd dx_atm(3, 2), dv_atm(3, 2);
    Eigen::MatrixXd dx_temp(3, 2), dv_temp(3, 2); // Температурные поправки (если есть)

    Eigen::MatrixXd xsta_j2000t(3, 2), vsta_j2000t(3, 2), asta_j2000(3, 2); // Итоговые координаты/скорости станций
    Eigen::Vector3d base_line; // Вектор базовой линии
    Eigen::VectorXd b_cfs(6); // Коэффициенты базовой линии
    std::vector<Eigen::Vector3d> k_star; // Векторы источников

    // *** 1. ПОДГОТОВИТЕЛЬНЫЕ РАСЧЕТЫ (до цикла) ***

    // 1.1 Расчет геодезических и геоцентрических координат (SITE)
    // Вызов функции для расчета стационарных параметров станций
    site(stations, space_stations, n_stations, dyear, site_xyz, site_vel, lat_geod, h_geod, lat_gcen, lon_gcen, sph_rad, u_site, v_site, vw);

    // 1.2 Расчет единичных векторов источников (SOURCE_VEC40)
    source_vec(sources, t_mean, k_star);

    // *** 2. ОСНОВНОЙ ЦИКЛ ОБРАБОТКИ (по наблюдениям) ***

    for (int k = 0; k < nobs; ++k) {
        const auto& obs = observations[k];
        j1 = obs.sta1; // Индекс 1-й станции (0-based)
        j2 = obs.sta2; // Индекс 2-й станции (0-based)
        j3 = obs.sou;  // Индекс источника (0-based)
        mjd = (double)obs.mjd;
        utc = obs.utc;

        // 2.1 ВРЕМЯ (TAITIME40, FUND_ARG40)
        tai_time(mjd, utc, tai, tt); //
        fund_arg(mjd, utc, ct, f, fd); //

        // 2.2 ЭФЕМЕРИДЫ (JPLEPH_421)
        double jd = mjd + 2400000.5;
        double ct = tt;
        double cent = (jd - cnst::JD2000) / cnst::JUL_CENT;
        jpleph(jd, ct, earth, sun, moon);

        // 2.3 EOP ИНТЕРПОЛЯЦИЯ И ТИДАЛЬНЫЕ ПОПРАВКИ (INTERP_EOP40, UT1R_2010)
        // В Фортране INTERP_EOP40 часто вызывает TERMS_71/TERMS_Lib.
        int k_int = 0;
        interp_eop40(k_int, mjd, utc, tt, ut1, eop_int, deop_int, arg_oc_tide, deop_diu, deop_lib, eop_data);
        
        // UT1R_2010 для расчета поправки UT1 от приливных вариаций
        double dut, dlod, domega;
        ut1r_2010(f, dut, dlod, domega); //

        // 2.4 МАТРИЦА ПРЕОБРАЗОВАНИЯ R2000
        r2000_matrix(mjd, ut1, eop_int, deop_int, i_choice, r2000, gast); // Предполагается C++ функция

        // 2.5 ПОПРАВКИ КООРДИНАТ СТАНЦИЙ (Приливы и нагрузки)
        // SITE_TIDE_SOLID, SITE_TIDE_OC, POLE_TIDE
        Eigen::Vector3d dx_tide1, dv_tide1;
        dx_tide1 = dxtide.col(0); dv_tide1 = dvtide.col(0);
        SITE_TIDE_SOLID(site_xyz[j1], lat_gcen[j1], lon_gcen[j1], sun, moon, f, fd, vw[j1], gast, r2000[0], dx_tide1, dv_tide1);
        dxtide.col(0) = dx_tide1; dvtide.col(0) = dv_tide1;

        Eigen::Vector3d dx_tide2, dv_tide2;
        dx_tide2 = dxtide.col(1); dv_tide2 = dvtide.col(1);
        SITE_TIDE_SOLID(site_xyz[j2], lat_gcen[j2], lon_gcen[j2], sun, moon, f, fd, vw[j2], gast, r2000[0], dx_tide2, dv_tide2);
        dxtide.col(1) = dx_tide2; dvtide.col(1) = dv_tide2;

        // Заглушка для tide_data
        OceanTideData tide_data;
        tide_data.amplitudes.setZero();
        tide_data.phases.setZero();
        Eigen::Vector3d dx_oc1, dv_oc1;
        dx_oc1 = dx_octide.col(0); dv_oc1 = dv_octide.col(0);
        SITE_TIDE_OC(mjd, ut1, tide_data, vw[j1], r2000[0], dx_oc1, dv_oc1);
        dx_octide.col(0) = dx_oc1; dv_octide.col(0) = dv_oc1;

        Eigen::Vector3d dx_oc2, dv_oc2;
        dx_oc2 = dx_octide.col(1); dv_oc2 = dv_octide.col(1);
        SITE_TIDE_OC(mjd, ut1, tide_data, vw[j2], r2000[0], dx_oc2, dv_oc2);
        dx_octide.col(1) = dx_oc2; dv_octide.col(1) = dv_oc2;

        Eigen::Vector3d dx_pol1, dv_pol1;
        dx_pol1 = dx_poltide.col(0); dv_pol1 = dv_poltide.col(0);
        POLE_TIDE(cent, lat_geod[j1], lon_gcen[j1], eop_int(1), eop_int(2), deop_int(1), deop_int(2), vw[j1], r2000[0], dx_pol1, dv_pol1);
        dx_poltide.col(0) = dx_pol1; dv_poltide.col(0) = dv_pol1;

        Eigen::Vector3d dx_pol2, dv_pol2;
        dx_pol2 = dx_poltide.col(1); dv_pol2 = dv_poltide.col(1);
        POLE_TIDE(cent, lat_geod[j2], lon_gcen[j2], eop_int(1), eop_int(2), deop_int(1), deop_int(2), vw[j2], r2000[0], dx_pol2, dv_pol2);
        dx_poltide.col(1) = dx_pol2; dv_poltide.col(1) = dv_pol2;
        
        // SITE_ATM40 (расчет атмосферной нагрузки, обычно для пары станций)
        Eigen::Vector2d dPdt = Eigen::Vector2d::Zero(); // Заглушка
        Eigen::Vector3d dx_atm1, dv_atm1;
        SITE_ATM40(dPdt(0), vw[j1], r2000[0], dx_atm1, dv_atm1);
        dx_atm.col(0) = dx_atm1; dv_atm.col(0) = dv_atm1;
        // Для второй станции - заглушка
        dx_atm.col(1).setZero(); dv_atm.col(1).setZero();

        // 2.6 ОБЪЕДИНЕНИЕ ПОПРАВОК (SITE_INST)
        Eigen::Vector3d xsta_out1, vsta_out1, asta_out1;
        SITE_INST(stations[j1].xyz, Eigen::Vector3d::Zero(), r2000[0], dx_tide1, dv_tide1, dx_oc1, dv_oc1, dx_pol1, dv_pol1, dx_atm1, dv_atm1, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), xsta_out1, vsta_out1, asta_out1);
        xsta_j2000t.col(0) = xsta_out1; vsta_j2000t.col(0) = vsta_out1; asta_j2000.col(0) = asta_out1;
        Eigen::Vector3d dx_atm2 = dx_atm.col(1), dv_atm2 = dv_atm.col(1);
        Eigen::Vector3d xsta_out2, vsta_out2, asta_out2;
        SITE_INST(stations[j2].xyz, Eigen::Vector3d::Zero(), r2000[0], dx_tide2, dv_tide2, dx_oc2, dv_oc2, dx_pol2, dv_pol2, dx_atm2, dv_atm2, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), xsta_out2, vsta_out2, asta_out2);
        xsta_j2000t.col(1) = xsta_out2; vsta_j2000t.col(1) = vsta_out2; asta_j2000.col(1) = asta_out2;
        dx_atm.col(1) = dx_atm2; dv_atm.col(1) = dv_atm2;

        // 2.7 БАЗОВАЯ ЛИНИЯ (BASELINE)
        Eigen::MatrixXd base_line_vec(3, 2);
        Eigen::Vector3d b_cfs_vec;
        baseline(r2000[0], xsta_j2000t, vsta_j2000t, asta_j2000, base_line_vec, b_cfs_vec);
        base_line = base_line_vec.col(0) - base_line_vec.col(1); // Вектор базовой линии
        b_cfs = b_cfs_vec;

        // 2.8 АБЕРРАЦИЯ И ТРОПОСФЕРА
        AberrationResult aber_res;
        Eigen::MatrixXd earth_mat(3, 3);
        earth_mat.col(0) = earth; // position
        earth_mat.col(1) = Eigen::Vector3d::Zero(); // velocity, заглушка
        earth_mat.col(2) = Eigen::Vector3d::Zero(); // acceleration, заглушка
        std::vector<Eigen::Vector3d> vsta_vec = {xsta_j2000t.col(0), xsta_j2000t.col(1)};
        aber_source(observations[k], r2000[0], lat_geod[j1], lat_geod[j2], h_geod[j1], h_geod[j2], k_star[j3], earth_mat, vsta_vec, vw, jd, ct, aber_res.elevation, aber_res.azimuth);

        Eigen::MatrixXd datmc_d_mat(2, 2), datmc_w_mat(2, 2), datmp_hmf(2, 2), datmp_wmf(2, 2), dgrad_n(2, 2), dgrad_e(2, 2), zen_dry_mat(2, 2), zen_wet_mat(2, 2);
        trop_delay(observations[k], jd, ct, stations[j1], stations[j2], aber_res.elevation, aber_res.azimuth,
                   datmc_d_mat, datmc_w_mat, datmp_hmf, datmp_wmf, dgrad_n, dgrad_e, zen_dry_mat, zen_wet_mat);
        
        // 2.9 ТЕОРЕТИЧЕСКАЯ ЗАДЕРЖКА (THEOR_DELAY4_10)
        double t2_t1, dt2_t1;
        double dtau_off = 0.0; 
        double dt_temp = 0.0;
        Eigen::VectorXd datmc_d_vec = datmc_d_mat.col(0);
        Eigen::VectorXd datmc_w_vec = datmc_w_mat.col(0);
        
        Eigen::MatrixXd sun_mat(3, 3), moon_mat(3, 3), dtau_off_mat(2, 2);
        sun_mat.col(0) = sun; sun_mat.col(1) = Eigen::Vector3d::Zero(); sun_mat.col(2) = Eigen::Vector3d::Zero();
        moon_mat.col(0) = moon; moon_mat.col(1) = Eigen::Vector3d::Zero(); moon_mat.col(2) = Eigen::Vector3d::Zero();
        theor_delay({base_line}, {xsta_j2000t.col(0), xsta_j2000t.col(1)}, {vsta_j2000t.col(0), vsta_j2000t.col(1)}, {asta_j2000.col(0), asta_j2000.col(1)}, k_star[j3], earth, sun_mat, moon_mat, datmc_d_mat, datmc_w_mat, dtau_off_mat, t2_t1, dt2_t1);

        // 2.10 ПОПРАВКА ОСИ ТЕЛЕСКОПА (MOUNT_TEL)
        Eigen::MatrixXd doff_dl(2, 2), d_dax(2, 2);
        mount_tel(observations[k], r2000[0], stations, k_star, vw, aber_res.elevation, aber_res.azimuth, doff_dl, d_dax, dtau_off_mat);
    }
}
}