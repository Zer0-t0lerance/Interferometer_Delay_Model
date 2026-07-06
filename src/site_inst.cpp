#include "functions.h"

namespace ariadna {

void SITE_INST(const Eigen::Vector3d& xsta_itrf,
               const Eigen::Matrix3d& r2000, const Eigen::Matrix3d& dr2000, const Eigen::Matrix3d& d2r2000,
               const Eigen::Vector3d& dx_tide, const Eigen::Vector3d& dv_tide,
               const Eigen::Vector3d& dx_octide, const Eigen::Vector3d& dv_octide,
               const Eigen::Vector3d& dx_poltide, const Eigen::Vector3d& dv_poltide,
               const Eigen::Vector3d& dx_atm, const Eigen::Vector3d& dv_atm,
               const Eigen::Vector3d& dx_temp, const Eigen::Vector3d& dv_temp,
               Eigen::Vector3d& xsta_j2000t, Eigen::Vector3d& vsta_j2000t,
               Eigen::Vector3d& asta_j2000) 
{
    // 1. Все поправки УЖЕ вычислены в системе J2000.0 в соответствующих подпрограммах.
    // Поэтому мы просто суммируем их.
    Eigen::Vector3d dx_total = dx_tide + dx_octide + dx_poltide + dx_atm + dx_temp;
    Eigen::Vector3d dv_total = dv_tide + dv_octide + dv_poltide + dv_atm + dv_temp;

    // 2. Вращаем только базовые координаты станции из ITRF в систему J2000.0
    Eigen::Vector3d xsta_j2000 = r2000 * xsta_itrf;
    Eigen::Vector3d vsta_j2000 = dr2000 * xsta_itrf;
    Eigen::Vector3d asta_base  = d2r2000 * xsta_itrf;

    // 3. Формируем итоговые векторы
    xsta_j2000t = xsta_j2000 + dx_total;
    
    vsta_j2000t = vsta_j2000 + dv_total;
    
    // В радиоинтерферометрии ускорения от приливных деформаций ничтожно малы,
    // поэтому остается только базовое центростремительное ускорение Земли.
    asta_j2000  = asta_base;
}

} // namespace ariadna