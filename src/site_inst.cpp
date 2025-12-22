// SITE_INST.cpp

#include "functions.h"
#include "constants.h"
#include <Eigen/Dense>

namespace ariadna {

void SITE_INST(const Eigen::Vector3d& xsta_itrf, const Eigen::Vector3d& vsta_itrf,
               const Eigen::MatrixXd& r2000,
               const Eigen::Vector3d& dx_tide, const Eigen::Vector3d& dv_tide,
               const Eigen::Vector3d& dx_octide, const Eigen::Vector3d& dv_octide,
               const Eigen::Vector3d& dx_poltide, const Eigen::Vector3d& dv_poltide,
               const Eigen::Vector3d& dx_atm, const Eigen::Vector3d& dv_atm,
               const Eigen::Vector3d& dx_temp, const Eigen::Vector3d& dv_temp,
               Eigen::Vector3d& xsta_j2000t, Eigen::Vector3d& vsta_j2000t,
               Eigen::Vector3d& asta_j2000) {

    using Eigen::Vector3d;
    using Eigen::Matrix3d;

    // --- 1. Позиция (Xsta) и Скорость (Vsta) ITRF -> J2000.0 ---
    // Формулы перехода:
    // X_J2000 = R * X_ITRF
    // V_J2000 = R_dot * X_ITRF + R * V_ITRF
    
    // Извлечение матриц R и R_dot из r2000
    const Matrix3d R = r2000.block<3, 3>(0, 0);       // R2000(3,3,0)
    const Matrix3d R_dot = r2000.block<3, 3>(0, 3);   // R2000(3,3,1)
    //const Matrix3d R_ddot = r2000.block<3, 3>(0, 6); // R2000(3,3,2) (для ускорения)

    // A. Преобразование базовой ITRF-позиции и скорости в J2000.0
    Vector3d x_base_j2000 = R * xsta_itrf;
    Vector3d v_base_j2000 = R_dot * xsta_itrf + R * vsta_itrf;

    // --- 2. Суммирование всех поправок к позиции (Xsta) ---
    // Fortran-код: xsta_j2000(i) = x_base_j2000(i) + dxtide(i,j1) + dx_octide(i,j1) + ...
    
    xsta_j2000t = x_base_j2000 
                  + dx_tide 
                  + dx_octide 
                  + dx_poltide 
                  + dx_atm 
                  + dx_temp;

    // --- 3. Суммирование всех поправок к скорости (Vsta) ---
    // Fortran-код: vsta_j2000(i) = v_base_j2000(i) + dvtide(i,j1) + dv_octide(i,j1) + ...

    vsta_j2000t = v_base_j2000 
                  + dv_tide 
                  + dv_octide 
                  + dv_poltide 
                  + dv_atm 
                  + dv_temp;
    
    // --- 4. Ускорение (Asta) ---
    // Fortran-код SITE_INST вычисляет ускорение (по умолчанию 0),
    // но его точный расчет требует R_ddot и производных всех поправок.
    // Поскольку ускорение не передается как входная поправка, устанавливаем 0.
    // Если требуются поправки ускорения, их необходимо передать и просуммировать.
    asta_j2000 = Vector3d::Zero();
}

} // namespace ariadna