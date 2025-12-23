#pragma once
#include "structures.h"

namespace ariadna {
// Time conversions
void ut1r_2010(const Eigen::VectorXd& f, double& dut, double& dlod, double& domega);
void t_eph(const Observation& obs, double tai, double ut1, double tt, double lon_gcen, double u_site, double v_site, double& ct, double& dtaidct);
void interp_eop(int k_int, const Observation& obs, double tt, double& ut1,
                Eigen::VectorXd& eop_int, Eigen::VectorXd& deop_int,
                Eigen::VectorXd& arg_oc_tide, Eigen::MatrixXd& deop_diu,
                Eigen::MatrixXd& deop_lib, const std::vector<EOPData>& eop_data);
void interp_eop40(int k_int, double mjd, double utc, double tt, double& ut1,
                  Eigen::VectorXd& eop_int, Eigen::VectorXd& deop_int,
                  Eigen::VectorXd& arg_oc_tide, Eigen::MatrixXd& deop_diu,
                  Eigen::MatrixXd& deop_lib, const std::vector<EOPData>& eop_data);
void tai_time(double mjd, double UTC, double &TAI, double &TT);
void fund_arg(double jd, double ct, double& cent, Eigen::VectorXd& f, Eigen::VectorXd& fd);
void terms_71(double cent, const Eigen::VectorXd& f, const Eigen::VectorXd& fd, Eigen::MatrixXd& dEOP_diu, Eigen::VectorXd& arg_oc_tide);
void terms_lib(double cent, const Eigen::VectorXd& f, const Eigen::VectorXd& fd, Eigen::MatrixXd& dEOP_lib);
void nsec(double mjd, double& idelt);
void interp_iers(const Observation& obs, Eigen::VectorXd& eop_int);

// Coordinate and transformation calculations
void site(const std::vector<Station>& stations, int j1, int j2, const Observation& observation,
          double cent, const Eigen::VectorXd& f, const Eigen::VectorXd& fd,
          double gast, const Eigen::MatrixXd& r2000_full,
          const Eigen::VectorXd& eop_int, const Eigen::MatrixXd& deop_diu,
          const Eigen::MatrixXd& deop_lib,
          std::vector<Eigen::Vector3d>& xsta_j2000t, std::vector<Eigen::Vector3d>& vsta_j2000t,
          std::vector<Eigen::Vector3d>& asta_j2000);
void site_atm40(int j1, int j2,
    const std::vector<Station>& stations,
    const Observation& observation,
    const Eigen::Vector2d& dPdt,
    const std::vector<Eigen::Matrix3d>& vw,
    const Eigen::MatrixXd& r2000,
    Eigen::Matrix<double, 3, 2>& dx_atm,
    Eigen::Matrix<double, 3, 2>& dv_atm);
void site_inst(const Eigen::Vector3d& xsta_itrf, const Eigen::Vector3d& vsta_itrf,
               const Eigen::MatrixXd& r2000,
               const Eigen::Vector3d& dx_tide, const Eigen::Vector3d& dv_tide,
               const Eigen::Vector3d& dx_octide, const Eigen::Vector3d& dv_octide,
               const Eigen::Vector3d& dx_poltide, const Eigen::Vector3d& dv_poltide,
               const Eigen::Vector3d& dx_atm, const Eigen::Vector3d& dv_atm,
               const Eigen::Vector3d& dx_temp, const Eigen::Vector3d& dv_temp,
               Eigen::Vector3d& xsta_j2000t, Eigen::Vector3d& vsta_j2000t,
               Eigen::Vector3d& asta_j2000);

void site(const std::vector<Station>& stations, const std::vector<SpaceStation>& space_stations, int n_stations, double dyear,
          std::vector<Eigen::Vector3d>& site_xyz, std::vector<Eigen::Vector3d>& site_vel, std::vector<double>& lat_geod,
          std::vector<double>& h_geod, std::vector<double>& lat_gcen,
          std::vector<double>& lon_gcen, std::vector<double>& sph_rad,
          std::vector<double>& u_site, std::vector<double>& v_site,
          std::vector<Eigen::Matrix3d>& vw);
void source_vec(const std::vector<Source>& sources, double t_mean, std::vector<Eigen::Vector3d>& k_star);
void jpl_eph(double jd, double ct, Eigen::Matrix3d& earth, Eigen::MatrixXd& sun, Eigen::MatrixXd& moon);
void jpleph(double jd, double ct, Eigen::Vector3d& earth, Eigen::Vector3d& sun, Eigen::Vector3d& moon);
void fund_arg(double jd, double ct, double& cent, Eigen::VectorXd& f, Eigen::VectorXd& fd);
void eps_a06(double jd, double ct, double& eps2000, double& eps_p03_2000, Eigen::VectorXd& e_mn);
void prec_matrix(double jd, double ct, double eps_p03_2000, Eigen::MatrixXd& pr, Eigen::MatrixXd& dpdp_ls, Eigen::MatrixXd& dpdp_pl, Eigen::MatrixXd& dpdp_om);
void bias(double eps2000, Eigen::Matrix3d& bias_matr);
void iau_2000_2006(double jd, double ct, const Eigen::VectorXd& f, const Eigen::VectorXd& fd, double eop_int_x, double deop_int_x, double eop_int_y, double deop_int_y, double eps_p03_2000, const Eigen::VectorXd& e_mn, Eigen::VectorXd& dpsir, Eigen::VectorXd& depsr, Eigen::VectorXd& eps, Eigen::MatrixXd& rn, Eigen::MatrixXd& dndpsi, Eigen::MatrixXd& dndeps);
void gas_time(double jd, double ut1, double ct, const Eigen::VectorXd& f, const Eigen::VectorXd& fd, const Eigen::VectorXd& dpsir, const Eigen::VectorXd& depsr, const Eigen::VectorXd& e_mn, const Eigen::VectorXd& eps, double dtaidct, double deop_int_ut1, double& diurnv, Eigen::VectorXd& gast, Eigen::MatrixXd& rs, Eigen::MatrixXd& drsdp_ls, Eigen::MatrixXd& drsdp_pl);
void wobble(double cent, double eop_x, double eop_y, double deop_x, double deop_y, Eigen::MatrixXd& ryx, Eigen::Matrix3d& ydxdx, Eigen::Matrix3d& dydyx, Eigen::Matrix3d& ddxdyx, Eigen::Matrix3d& ddydyx);
void r2000_matrix(double mjd, double ut1, const Eigen::VectorXd& eop_int, const Eigen::VectorXd& deop_int, int i_choice, Eigen::Matrix3d r2000[3], double& gast);

// Tidal and atmospheric effects
void therm_def(const Station& station, const Observation& obs, double dtdt, const Eigen::Matrix3d& vw, const Eigen::MatrixXd& r2000, Eigen::Vector3d& dx_temp, Eigen::Vector3d& dv_temp);
     /**
     * @brief Преобразование декартовых (экваториальная и полярная компоненты) 
     * координат в геодезические (широта и высота) по методу Борковского (1989).
     * Соответствует SUBROUTINE GEOD (SITE_CORR.FOR.TXT).
     * @param equatorial_radius_r Экваториальная компонента (sqrt(X^2 + Y^2)) [m]
     * @param z_polar Полярная компонента Z [m]
     * @param geodetic_latitude_fi Геодезическая широта [rad] (Выходной параметр)
     * @param geodetic_height_h Геодезическая высота [m] (Выходной параметр)
     */
    void GEOD(double equatorial_radius_r, double z_polar, 
              double& geodetic_latitude_fi, double& geodetic_height_h);
    /**
     * @brief Вычисляет смещение и скорость станции, вызванные приливом полюса
     * (Pole Tide) в системе J2000.0 (IERS Conventions 2000).
     * * @param cent Эпоха в юлианских столетиях с J2000.0.
     * @param lat_geod Геодезическая широта станции [rad].
     * @param lon_geod Восточная долгота станции [rad].
     * @param xp Фактическая координата полюса X [rad].
     * @param yp Фактическая координата полюса Y [rad].
     * @param xp_rate Скорость движения полюса X [rad/day].
     * @param yp_rate Скорость движения полюса Y [rad/day].
     * @param vw_i Матрица преобразования VEN (локальная) -> ITRF (краст-фиксированная).
     * @param r2000 Матрица преобразования ITRF -> J2000.0 и ее производные (r2000(3,3,0)=R, r2000(3,3,1)=R_dot, r2000(3,3,2)=R_ddot).
     * @param dx_poltide Выход: Вектор смещения (X,Y,Z) в J2000.0 [m].
     * @param dv_poltide Выход: Вектор скорости (Vx,Vy,Vz) в J2000.0 [m/s].
     */
    void POLE_TIDE(double cent, double lat_geod, double lon_geod,
                   double xp, double yp, double xp_rate, double yp_rate,
                   const Eigen::Matrix3d& vw_i, const Eigen::MatrixXd& r2000,
                   Eigen::Vector3d& dx_poltide, Eigen::Vector3d& dv_poltide);

    /**
     * @brief Вычисляет смещение и скорость станции, вызванные твердыми земными приливами (Solid Earth Tide) в J2000.0 системе.
     * Соответствует SUBROUTINE SITE_TIDE_SOLID (SITE_TIDE_SOLID.for.txt).
     * @param xsta Геоцентрическое положение станции в ITRF [m].
     * @param lat_gcen Геоцентрическая широта станции [rad].
     * @param lon_gcen Восточная долгота станции [rad].
     * @param sun Вектор положения Солнца в ITRF [m].
     * @param moon Вектор положения Луны в ITRF [m].
     * @param f Вектор угловых аргументов приливов (f_i) [rad].
     * @param fd Вектор производных угловых аргументов приливов (df_i/dt) [rad/s].
     * @param vw_i Матрица преобразования VEN (локальная) -> ITRF (краст-фиксированная).
     * @param gast Звёздное время Гринвича GAST [rad].
     * @param r2000 Матрица преобразования ITRF -> J2000.0 и ее производные.
     * @param dxtide Выход: Вектор смещения (X,Y,Z) в J2000.0 [m].
     * @param dvtide Выход: Вектор скорости (Vx,Vy,Vz) в J2000.0 [m/s].
     */
    void SITE_TIDE_SOLID(const Eigen::Vector3d& xsta, double lat_gcen, double lon_gcen,
                         const Eigen::Vector3d& sun, const Eigen::Vector3d& moon,
                         const Eigen::VectorXd& f, const Eigen::VectorXd& fd,
                         const Eigen::Matrix3d& vw_i, double gast,
                         const Eigen::MatrixXd& r2000,
                         Eigen::Vector3d& dxtide, Eigen::Vector3d& dvtide);

    /**
     * @brief Вычисляет смещение и скорость станции, вызванные океаническими приливами (Ocean Tide Loading) в J2000.0 системе.
     * Соответствует SUBROUTINE SITE_TIDE_OC (SITE_TIDE_OC.FOR.TXT).
     * @param mjd_start MJD в 0 часов UTC.
     * @param ut1_sec UT1 время дня [s].
     * @param tide_data Амплитуды и фазы 11 приливных волн для станции.
     * @param vw_i Матрица преобразования VEN (локальная) -> ITRF (краст-фиксированная).
     * @param r2000 Матрица преобразования ITRF -> J2000.0 и ее производные.
     * @param dx_octide Выход: Вектор смещения (X,Y,Z) в J2000.0 [m].
     * @param dv_octide Выход: Вектор скорости (Vx,Vy,Vz) в J2000.0 [m/s].
     */
    void SITE_TIDE_OC(int mjd_start, double ut1_sec, const ariadna::OceanTideData& tide_data,
                      const Eigen::Matrix3d& vw_i, const Eigen::MatrixXd& r2000,
                      Eigen::Vector3d& dx_octide, Eigen::Vector3d& dv_octide);

    // Вспомогательная функция для расчета аргументов приливов
    void calc_tide_angles(int mjd_start, double ut1_sec, Eigen::VectorXd& angle, Eigen::VectorXd& speed_angle);

    /**
     * @brief Вычисляет смещение и скорость станции, вызванные атмосферной нагрузкой (Atmospheric Loading) в J2000.0 системе.
     * Соответствует SUBROUTINE SITE_ATM40 (SITE_ATM40.FOR.TXT).
     * @param dPdt Скорость изменения давления на станции [mbar/s].
     * @param vw_i Матрица преобразования VEN (локальная) -> ITRF (краст-фиксированная).
     * @param r2000 Матрица преобразования ITRF -> J2000.0 и ее производные.
     * @param dx_atm Выход: Вектор смещения (X,Y,Z) в J2000.0 [m].
     * @param dv_atm Выход: Вектор скорости (Vx,Vy,Vz) в J2000.0 [m/s].
     */
    void SITE_ATM40(double dPdt, 
                    const Eigen::Matrix3d& vw_i, const Eigen::MatrixXd& r2000,
                    Eigen::Vector3d& dx_atm, Eigen::Vector3d& dv_atm);

    /**
     * @brief Вычисляет итоговую позицию и скорость станции в системе J2000.0 
     * путем суммирования всех поправок (приливы, полюс, атмосфера, температура).
     * Соответствует SUBROUTINE SITE_INST (SITE_INST.FOR.TXT).
     * @param xsta_itrf Геоцентрическое положение станции в ITRF (краст-фиксированная) [m].
     * @param vsta_itrf Геоцентрическая скорость станции в ITRF [m/s].
     * @param r2000 Матрица преобразования ITRF -> J2000.0 и ее производные (R(3,3,0), R_dot(3,3,1)).
     * @param dx_tide Вектор смещения твердых приливов в J2000.0 [m].
     * @param dv_tide Вектор скорости твердых приливов в J2000.0 [m/s].
     * @param dx_octide Вектор смещения океанических приливов в J2000.0 [m].
     * @param dv_octide Вектор скорости океанических приливов в J2000.0 [m/s].
     * @param dx_poltide Вектор смещения прилива полюса в J2000.0 [m].
     * @param dv_poltide Вектор скорости прилива полюса в J2000.0 [m/s].
     * @param dx_atm Вектор смещения атмосферной нагрузки в J2000.0 [m].
     * @param dv_atm Вектор скорости атмосферной нагрузки в J2000.0 [m/s].
     * @param dx_temp Вектор смещения температурной деформации в J2000.0 [m] (заглушка).
     * @param dv_temp Вектор скорости температурной деформации в J2000.0 [m/s] (заглушка).
     * @param xsta_j2000t Выход: Итоговая позиция станции в J2000.0 [m].
     * @param vsta_j2000t Выход: Итоговая скорость станции в J2000.0 [m/s].
     * @param asta_j2000 Выход: Ускорение станции в J2000.0 [m/s^2] (заглушка).
     */
    void SITE_INST(const Eigen::Vector3d& xsta_itrf, const Eigen::Vector3d& vsta_itrf,
                   const Eigen::MatrixXd& r2000,
                   const Eigen::Vector3d& dx_tide, const Eigen::Vector3d& dv_tide,
                   const Eigen::Vector3d& dx_octide, const Eigen::Vector3d& dv_octide,
                   const Eigen::Vector3d& dx_poltide, const Eigen::Vector3d& dv_poltide,
                   const Eigen::Vector3d& dx_atm, const Eigen::Vector3d& dv_atm,
                   const Eigen::Vector3d& dx_temp, const Eigen::Vector3d& dv_temp,
                   Eigen::Vector3d& xsta_j2000t, Eigen::Vector3d& vsta_j2000t,
                   Eigen::Vector3d& asta_j2000);

// Delay and derivative calculations
void dmeteo1_dt(const std::vector<Observation>& observations, const std::vector<Station>& stations, double t_mean, std::vector<Eigen::Vector3d>& site_meteo, int& ndeg, std::vector<Eigen::VectorXd>& t_coef, std::vector<Eigen::VectorXd>& p_coef, std::vector<Eigen::VectorXd>& hum_coef);
void dmeteo2_dt(const Station& station, int n_stations, int ndeg, const Observation& obs, double t_mean, const std::vector<Eigen::VectorXd>& t_coef, const std::vector<Eigen::VectorXd>& p_coef, const std::vector<Eigen::VectorXd>& hum_coef, double& dtdt, double& dpdt, double& dhumdt);
//void aber_source(const Observation& obs, const Eigen::Matrix3d& r2000, double lat_geod1, double lat_geod2, double h_geod1, double h_geod2, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& vsta_j2000t, const std::vector<Eigen::Matrix3d>& vw, double jd, double ct, Eigen::MatrixXd& e, Eigen::MatrixXd& az);

/**
 * @brief Вычисляет положение источника, скорректированное за годовую и суточную аберрацию.
 * Соответствует SUBROUTINE ABER_SOURCE (ABER_SOURCE.FOR.TXT).
 * Рассчитывает видимые углы возвышения (Elevation) и азимута (Azimuth), а также их производные 
 * для двух станций наблюдения в топоцентрической системе координат (VEN).
 * * @param obs Структура с данными наблюдения (индексы станций и временные метки).
 * @param r2000 Вектор матриц вращения (3x3): 
 * [0] - Матрица перехода ITRF -> GCRS (J2000).
 * [1] - Первая производная матрицы (1/sec).
 * @param k_s Единичный вектор на источник в системе J2000.0.
 * @param earth Матрица характеристик Земли, где столбец index=1 (второй) — барицентрическая скорость Земли [m/s].
 * @param vsta_j2000t Вектор скоростей станций в системе J2000.0 (геоцентрические) [m/s].
 * @param vw Вектор матриц перехода (3x3) для каждой станции: VEN (локальная) -> ITRF (краст-фиксированная).
 * @param e Выход: Матрица 2x2 углов возвышения: e(i,0) - угол [rad], e(i,1) - производная [rad/s].
 * @param az Выход: Матрица 2x2 азимутов: az(i,0) - угол [rad], az(i,1) - производная [rad/s].
 */
void aber_source(
    const Observation& obs,
    const std::vector<Eigen::Matrix3d>& r2000, 
    const Eigen::Vector3d& k_s,
    const Eigen::Matrix<double, 3, 3>& earth,
    const std::vector<Eigen::Vector3d>& vsta_j2000t,
    const std::vector<Eigen::Matrix3d>& vw,
    Eigen::Matrix2d& e,
    Eigen::Matrix2d& az
);

void mount_tel(const Observation& obs, const Eigen::Matrix3d& r2000, const std::vector<Station>& stations, const std::vector<Eigen::Vector3d>& k_star, const std::vector<Eigen::Matrix3d>& vw, const Eigen::MatrixXd& e, const Eigen::MatrixXd& az, Eigen::MatrixXd& doff_dl, Eigen::MatrixXd& d_dax, Eigen::MatrixXd& dtau_off);
double sbend(double el_rad, double temp_k, double humid_f, double press_hg);

void baseline(const Eigen::Matrix3d& r2000, const Eigen::MatrixXd& xsta_j2000t, const Eigen::MatrixXd& vsta_j2000t, const Eigen::MatrixXd& asta_j2000, Eigen::MatrixXd& base_line, Eigen::Vector3d& b_cfs);
void uv_plane(const Source& source, const std::vector<Eigen::Vector3d>& base_line, const std::vector<Eigen::Vector3d>& xsta_j2000t, double scale, Eigen::Vector3d& uv_coor);

void trop_delay(const Observation& obs, double jd, double ct, const Station& sta1, const Station& sta2, const Eigen::MatrixXd& e, const Eigen::MatrixXd& az, Eigen::MatrixXd& datmc_d, Eigen::MatrixXd& datmc_w, Eigen::MatrixXd& datmp_hmf, Eigen::MatrixXd& datmp_wmf, Eigen::MatrixXd& dgrad_n, Eigen::MatrixXd& dgrad_e, Eigen::MatrixXd& zen_dry, Eigen::MatrixXd& zen_wet);
void nhmf2(double epoch, double latitude, double height, double elev, Eigen::Vector2d& hmf);
void nwmf2(double latitude, double elev, Eigen::Vector2d& wmf);
void sast_dry(double pres, double dot_pres, double lat_geod, double height, double dpdh, double& z_d, double& dot_z_d, double& dz_ddh);
void sast_wet(double rel_hum, double tc, double dot_rel_hum, double dot_tc, double& z_w, double& dot_z_w);

void theor_delay(const std::vector<Eigen::Vector3d>& base_line, const std::vector<Eigen::Vector3d>& xsta_j2000t, const std::vector<Eigen::Vector3d>& vsta_j2000t, const std::vector<Eigen::Vector3d>& asta_j2000, const Eigen::Vector3d& k_s, const Eigen::Vector3d& earth, const Eigen::MatrixXd& sun, const Eigen::MatrixXd& moon, const Eigen::MatrixXd& datmc_d, const Eigen::MatrixXd& datmc_w, const Eigen::MatrixXd& dtau_off, double& t2_t1, double& dt2_t1);
void theor_delay_orb(const std::vector<Eigen::Vector3d>& base_line, const std::vector<Eigen::Vector3d>& xsta_j2000t, const std::vector<Eigen::Vector3d>& vsta_j2000t, const std::vector<Eigen::Vector3d>& asta_j2000, const Eigen::Vector3d& k_s, const Eigen::Vector3d& earth, const Eigen::MatrixXd& sun, const Eigen::MatrixXd& moon, const Eigen::MatrixXd& datmc_d, const Eigen::MatrixXd& datmc_w, const Eigen::MatrixXd& dtau_off, double& t2_t1, double& dt2_t1);
void der_star(double jd, double ct, double dyear, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, Eigen::MatrixXd& dstar, Eigen::MatrixXd& dstar_rate);
void der_site(double dyear, const Eigen::MatrixXd& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, std::vector<Eigen::MatrixXd>& dsite, std::vector<Eigen::MatrixXd>& dsite_v);
void der_polar(const Eigen::MatrixXd& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, const Eigen::Vector3d& b_cfs, std::vector<Eigen::MatrixXd>& pr, std::vector<Eigen::MatrixXd>& rn, std::vector<Eigen::MatrixXd>& rs, std::vector<Eigen::MatrixXd>& ryx, std::vector<Eigen::MatrixXd>& ydxdx, std::vector<Eigen::MatrixXd>& dydyx, std::vector<Eigen::MatrixXd>& ddxdyx, std::vector<Eigen::MatrixXd>& ddydyx, std::vector<Eigen::VectorXd>& dx_pol_dx, std::vector<Eigen::VectorXd>& dx_pol_dy, const Eigen::VectorXd& arg_oc_tide, Eigen::MatrixXd& dwob, std::vector<Eigen::MatrixXd>& dx_aj, std::vector<Eigen::MatrixXd>& dx_bj, std::vector<Eigen::MatrixXd>& dy_aj, std::vector<Eigen::MatrixXd>& dy_bj);
void der_ut1(const Eigen::MatrixXd& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, const Eigen::Vector3d& b_cfs, const Eigen::VectorXd& gast, std::vector<Eigen::MatrixXd>& pr, std::vector<Eigen::MatrixXd>& rn, std::vector<Eigen::MatrixXd>& rs, std::vector<Eigen::MatrixXd>& ryx, const Eigen::VectorXd& arg_oc_tide, double diurnv, Eigen::VectorXd& dut1_tai, std::vector<Eigen::MatrixXd>& dut1_aj, std::vector<Eigen::MatrixXd>& dut1_bj);
void der_nut(const Eigen::MatrixXd& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, const Eigen::Vector3d& b_cfs, const Eigen::VectorXd& gast, std::vector<Eigen::MatrixXd>& pr, std::vector<Eigen::MatrixXd>& rn, std::vector<Eigen::MatrixXd>& rs, std::vector<Eigen::MatrixXd>& ryx, const Eigen::VectorXd& e_mn, const Eigen::MatrixXd& dndpsi, const Eigen::MatrixXd& dndeps, Eigen::MatrixXd& dnut);
void der_prec(const Eigen::MatrixXd& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, const Eigen::Vector3d& b_cfs, const Eigen::MatrixXd& pr, const Eigen::MatrixXd& rn, const Eigen::MatrixXd& rs, const Eigen::MatrixXd& ryx, const Eigen::MatrixXd& dpdp_ls, const Eigen::MatrixXd& dpdp_pl, const Eigen::MatrixXd& drsdp_ls, const Eigen::MatrixXd& drsdp_pl, Eigen::MatrixXd& dpr_lspl);
void der_love_number(const Station& sta1, const Station& sta2, const Eigen::MatrixXd& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, const std::vector<Eigen::Matrix3d>& drdh0_3, const std::vector<Eigen::Matrix3d>& drdh02_3, const std::vector<Eigen::Matrix3d>& drdrl0_3, const std::vector<Eigen::Matrix3d>& drdl02_3, const std::vector<Eigen::Matrix3d>& drdh3_3, const std::vector<Eigen::Matrix3d>& drdl3_3, const std::vector<Eigen::Matrix3d>& drdl1_1_2000_3, const std::vector<Eigen::Matrix3d>& drdl1_2_2000_3, const std::vector<Eigen::Matrix3d>& drdhi_1_2000_3, const std::vector<Eigen::Matrix3d>& drdhi_2_2000_3, const std::vector<Eigen::Matrix3d>& drdli_1_2000_3, const std::vector<Eigen::Matrix3d>& drdli_2_2000_3, Eigen::VectorXd& d_dh0, Eigen::VectorXd& d_dh02, Eigen::VectorXd& d_dl0, Eigen::VectorXd& d_dl02, Eigen::VectorXd& d_dh3, Eigen::VectorXd& d_dl3, Eigen::VectorXd& d_dl1_1, Eigen::VectorXd& d_dl1_2, Eigen::VectorXd& d_dhi_1, Eigen::VectorXd& d_dhi_2, Eigen::VectorXd& d_dli_1, Eigen::VectorXd& d_dli_2);
void der_atm_load(const Station& sta1, const Station& sta2, const Eigen::MatrixXd& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, const std::vector<Eigen::Matrix3d>& dr1_da2000, const std::vector<Eigen::Matrix3d>& dr1_db2000, const std::vector<Eigen::Matrix3d>& dv1_da2000, const std::vector<Eigen::Matrix3d>& dv1_db2000, const std::vector<Eigen::Matrix3d>& dr2_da2000, const std::vector<Eigen::Matrix3d>& dr2_db2000, const std::vector<Eigen::Matrix3d>& dv2_da2000, const std::vector<Eigen::Matrix3d>& dv2_db2000, const std::vector<Eigen::Vector3d>& dr_dreg2000, const std::vector<Eigen::Vector3d>& dv_dreg2000, std::vector<Eigen::MatrixXd>& dt_da, std::vector<Eigen::MatrixXd>& dt_db, std::vector<Eigen::MatrixXd>& df_da, std::vector<Eigen::MatrixXd>& df_db, Eigen::VectorXd& dt_dreg, Eigen::VectorXd& df_dreg);
void create_matr(const Observation& obs, const Station& sta1, const Station& sta2, int i_good, int nobs, int n_sources, int n_stations, double t_mean, double t2_t1, double dt2_t1, const Eigen::MatrixXd& dwob, const Eigen::VectorXd& dut1_tai, const Eigen::MatrixXd& dnut, const std::vector<Eigen::MatrixXd>& dsite, const std::vector<Eigen::MatrixXd>& dsite_v, const Eigen::MatrixXd& datmp_hmf, const Eigen::MatrixXd& datmp_wmf, const Eigen::MatrixXd& dgrad_n, const Eigen::MatrixXd& dgrad_e, const Eigen::MatrixXd& dstar, const Eigen::MatrixXd& dstar_rate, const Eigen::MatrixXd& dpr_lspl, const Eigen::VectorXd& d_dh0, const Eigen::VectorXd& d_dh02, const Eigen::VectorXd& d_dl0, const Eigen::VectorXd& d_dl02, const Eigen::VectorXd& d_dh3, const Eigen::VectorXd& d_dl3, const Eigen::VectorXd& d_dl1_1, const Eigen::VectorXd& d_dl1_2, const Eigen::VectorXd& d_dhi_1, const Eigen::VectorXd& d_dhi_2, const Eigen::VectorXd& d_dli_1, const Eigen::VectorXd& d_dli_2, const Eigen::MatrixXd& d_dax, const std::vector<Eigen::MatrixXd>& dt_da, const std::vector<Eigen::MatrixXd>& dt_db, const std::vector<Eigen::MatrixXd>& df_da, const std::vector<Eigen::MatrixXd>& df_db, const Eigen::VectorXd& dt_dreg, const Eigen::VectorXd& df_dreg, const std::vector<Eigen::MatrixXd>& dx_aj, const std::vector<Eigen::MatrixXd>& dx_bj, const std::vector<Eigen::MatrixXd>& dy_aj, const std::vector<Eigen::MatrixXd>& dy_bj, const std::vector<Eigen::MatrixXd>& dut1_aj, const std::vector<Eigen::MatrixXd>& dut1_bj, const Eigen::MatrixXd& zen_dry, const Eigen::MatrixXd& zen_wet, Eigen::MatrixXd& m_pd, Eigen::VectorXd& y, Eigen::VectorXd& w);
void integr8_asc(const Station& station, const Observation& obs, int nobs, int l_segm, int n_wr_tot, double delta_sec, const std::string& track_site, int mjd_beg, double utc_beg, double ct_beg, double dyear, const Eigen::MatrixXd& r2000, std::vector<Eigen::Vector3d>& xsta_j2000t, std::vector<Eigen::Vector3d>& vsta_j2000t, std::vector<Eigen::Vector3d>& asta_j2000, Eigen::MatrixXd& r_sat_pr, const Eigen::Vector3d& k_s, double& phase1);

// Main processing function
void process_ariadna(const std::vector<Station>& stations, const std::vector<Source>& sources, const std::vector<Observation>& observations, const std::vector<SpaceStation>& space_stations, const std::vector<OrbitData>& orbit_data, int n_segm, int k_ch_c, int k_ch_z, double delta_sec, const std::string& output_path, const std::vector<EOPData>& eop_data, double mjd_mean, double utc_mean, double t_mean);
} // namespace ariadna