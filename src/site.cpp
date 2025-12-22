#include "functions.h"

namespace ariadna {
void site(const std::vector<Station>& stations, int j1, int j2, const Observation& observation,
          double cent, const Eigen::VectorXd& f, const Eigen::VectorXd& fd,
          double gast, const Eigen::MatrixXd& r2000_full,
          const Eigen::VectorXd& eop_int, const Eigen::MatrixXd& deop_diu,
          const Eigen::MatrixXd& deop_lib,
          std::vector<Eigen::Vector3d>& xsta_j2000t, std::vector<Eigen::Vector3d>& vsta_j2000t,
          std::vector<Eigen::Vector3d>& asta_j2000) {

    using Eigen::Vector3d;
    using Eigen::Matrix3d;
    
    // Предполагается, что r2000_full (r2000(3,3,3) в Fortran) представлен как 3x9 матрица: R | R_dot | R_ddot
    Matrix3d r2000_m    = r2000_full.block<3, 3>(0, 0);
    Matrix3d r2000_dot_m = r2000_full.block<3, 3>(0, 3);
    Matrix3d r2000_ddot_m = r2000_full.block<3, 3>(0, 6);

    // Время, прошедшее с эпохи ITRF (предположительно J2000.0) в годах,
    // для обновления ITRF-координат (sta.vel в м/год)
    // cent (юлианских столетий) * 100.0 = лет
    double delta_t_years = cent * 100.0;
    
    // Итерация по двум станциям, участвующим в наблюдении (индексы 0 и 1 в выходных векторах)
    std::vector<int> site_indices = {j1, j2};

    for (size_t i = 0; i < site_indices.size(); ++i) {
        int site_idx = site_indices[i];
        const Station& sta = stations[site_idx];
        
        // 1. ITRF позиция и скорость на момент наблюдения
        // x_itrf(t) = x_itrf(t0) + v_itrf * (t - t0)
        Vector3d xsta_itrf = sta.xyz + sta.vel * delta_t_years;
        // Скорость в ITRF (м/с). sta.vel в м/год.
        Vector3d vsta_itrf = sta.vel / (cnst::SECDAY * 365.25);
        
        // 2. Инициализация векторов поправок (в системе ITRF)
        Vector3d dxtide = Vector3d::Zero();     // Solid Earth Tide Position
        Vector3d dvtide = Vector3d::Zero();     // Solid Earth Tide Velocity
        Vector3d dx_octide = Vector3d::Zero();  // Ocean Tide Position
        Vector3d dv_octide = Vector3d::Zero();  // Ocean Tide Velocity
        Vector3d dx_poltide = Vector3d::Zero(); // Pole Tide Position
        Vector3d dv_poltide = Vector3d::Zero(); // Pole Tide Velocity
        Vector3d dx_atm = Vector3d::Zero();     // Atmospheric Loading Position
        Vector3d dv_atm = Vector3d::Zero();     // Atmospheric Loading Velocity
        Vector3d dx_temp = Vector3d::Zero();    // Thermal Deformation Position (Placeholder, т.к. файл не предоставлен, но упоминается в SITE_INST)
        Vector3d dv_temp = Vector3d::Zero();    // Thermal Deformation Velocity

        // ==========================================================
        // 3. Вызов подпрограмм коррекции (ПОКА ЗАГЛУШКИ!)
        // ==========================================================

        // a. Solid Earth Tide (SITE_TIDE_SOLID)
        // site_tide_solid(xsta_itrf, ..., dxtide, dvtide, ...);

        // b. Pole Tide (POLE_TIDE)
        // pole_tide(cent, sta.lat_geod, sta.lon_geod, xp, yp, dxdt, dydt, ..., dx_poltide, dv_poltide, ...);

        // c. Ocean Tide Loading (SITE_TIDE_OC)
        // site_tide_oc(sta_index, ..., dx_octide, dv_octide);

        // d. Atmospheric Loading (SITE_ATM40)
        
        site_atm40(j1, j2, stations, observation,
                    dPdt, vw, r2000, dx_atm, dv_atm);

    Eigen::Matrix<double, 3, 2>& dv_atm);

        // e. Thermal Deformation (необходима реализация)
        // thermal_deformation(..., dx_temp, dv_temp);

        // ==========================================================
        
        // 4. Суммирование поправок и вращение в J2000.0 (SITE_INST)
        site_inst(xsta_itrf, vsta_itrf,
                  r2000_m, r2000_dot_m, r2000_ddot_m,
                  dxtide, dvtide, dx_octide, dv_octide, dx_poltide, dv_poltide,
                  dx_atm, dv_atm, dx_temp, dv_temp,
                  xsta_j2000t[i], vsta_j2000t[i], asta_j2000[i]);
    }
}
} // namespace ariadna