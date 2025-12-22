// SITE_TIDE_OC.cpp
#include "functions.h"

namespace ariadna {

void SITE_TIDE_OC(int mjd_start, double ut1_sec, const OceanTideData& tide_data,
                  const Eigen::Matrix3d& vw_i, const Eigen::MatrixXd& r2000,
                  Eigen::Vector3d& dx_octide, Eigen::Vector3d& dv_octide) {
    
    using namespace cnst;
    using Eigen::Vector3d;
    using Eigen::Matrix3d;
    using Eigen::VectorXd;
    using std::sin;
    using std::cos;
    
    // --- 1. Расчет угловых аргументов приливов ---
    VectorXd angle;      // Угловые аргументы [rad]
    VectorXd speed_angle; // Скорости приливных аргументов [rad/s]
    calc_tide_angles(mjd_start, ut1_sec, angle, speed_angle);

    // --- 2. Инициализация смещения и скорости в ITRF (VEN) ---
    // VEN-система: (Up, North, East)
    Vector3d dr_ven = Vector3d::Zero();
    Vector3d dv_ven = Vector3d::Zero();

    // --- 3. Суммирование вклада 11 волн ---
    for (int k = 0; k < NUM_TIDES; ++k) {
        // Текущий угол и его скорость
        double current_angle = angle(k);
        double current_speed = speed_angle(k);

        // Вклад k-ой волны в смещение (в VEN-системе)
        // dx_ven = Amp * cos(Angle - Phase)
        
        // i=0: Up, i=1: North, i=2: East
        for (int i = 0; i < 3; ++i) {
            double amplitude = tide_data.amplitudes(i, k); // [m]
            double phase = tide_data.phases(i, k);       // [rad]
            
            // Смещение dr_ven(i) = Amplitude * cos(Angle - Phase)
            dr_ven(i) += amplitude * cos(current_angle - phase);
            
            // Скорость dv_ven(i) = d/dt (Amplitude * cos(Angle - Phase))
            // dv/dt = Amplitude * (-sin(Angle - Phase)) * d(Angle)/dt
            dv_ven(i) -= amplitude * sin(current_angle - phase) * current_speed;
        }
    }

    // --- 4. Преобразование из VEN (Up, North, East) в ITRF (краст-фиксированная X,Y,Z) ---
    // dr_itrf = VW_i * dr_ven
    Vector3d dx_itrf = vw_i * dr_ven;
    Vector3d dv_itrf = vw_i * dv_ven;

    // --- 5. Преобразование из ITRF в J2000.0 ---
    
    // R - матрица ITRF -> J2000.0 (R2000(3,3,0))
    const Matrix3d R = r2000.block<3, 3>(0, 0); 
    dx_octide = R * dx_itrf;

    // R_dot - производная R (R2000(3,3,1))
    const Matrix3d R_dot = r2000.block<3, 3>(0, 3);
    
    // dv_j2000 = R_dot * dx_itrf + R * dv_itrf
    dv_octide = R_dot * dx_itrf + R * dv_itrf;
}

} // namespace ariadna