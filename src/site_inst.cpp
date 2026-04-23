#include "functions.h"

namespace ariadna {

void site_inst(const Eigen::Vector3d& xsta_itrf, const Eigen::Vector3d& vsta_itrf,
               const Eigen::MatrixXd& r2000,
               const Eigen::Vector3d& dx_tide, const Eigen::Vector3d& dv_tide,
               const Eigen::Vector3d& dx_octide, const Eigen::Vector3d& dv_octide,
               const Eigen::Vector3d& dx_poltide, const Eigen::Vector3d& dv_poltide,
               const Eigen::Vector3d& dx_atm, const Eigen::Vector3d& dv_atm,
               const Eigen::Vector3d& dx_temp, const Eigen::Vector3d& dv_temp,
               Eigen::Vector3d& xsta_j2000t, Eigen::Vector3d& vsta_j2000t,
               Eigen::Vector3d& asta_j2000) {

    // A. Преобразование базовой ITRF-позиции и скорости в J2000.0
    Eigen::Vector3d x_base_j2000 = r2000.block<3, 3>(0, 0) * xsta_itrf;
    Eigen::Vector3d v_base_j2000 = r2000.block<3, 3>(3, 0) * xsta_itrf + r2000.block<3, 3>(0, 0) * vsta_itrf;

    // B. Суммирование всех поправок (Приливы Земли, Океан, Полюс, Атмосфера, Температура)
    xsta_j2000t = x_base_j2000 + dx_tide + dx_octide + dx_poltide + dx_atm + dx_temp;
    vsta_j2000t = v_base_j2000 + dv_tide + dv_octide + dv_poltide + dv_atm + dv_temp;
    
    // C. Ускорение (заглушка, обычно используется только для спутников)
    asta_j2000.setZero();
}
}