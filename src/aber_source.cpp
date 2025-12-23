#include "functions.h"

namespace ariadna {

void aber_source(
    const Observation& obs,
    const std::vector<Eigen::Matrix3d>& r2000, 
    const Eigen::Vector3d& k_s,
    const Eigen::Matrix<double, 3, 3>& earth,
    const std::vector<Eigen::Vector3d>& vsta_j2000t,
    const std::vector<Eigen::Matrix3d>& vw,
    Eigen::Matrix2d& e,
    Eigen::Matrix2d& az) 
{
    int sites[2] = {obs.sta1, obs.sta2};

    // Проходим по обеим станциям (базис)
    for (int j = 0; j < 2; ++j) {
        int idx = sites[j];

        // --- 1. РАСЧЕТ ПОЛНОЙ СКОРОСТИ И АБЕРРАЦИИ ---
        // Суммируем барицентрическую скорость Земли и геоцентрическую скорость станции
        Eigen::Vector3d v_total = earth.col(1) + vsta_j2000t[j];
        
        // Математика согласно IERS / формуле (6-66)
        double v_inv_c = v_total.norm() / cnst::C;
        double v_dot_s_c = k_s.dot(v_total) / cnst::C;
        
        // Релятивистский множитель Лоренца
        double rel_p = 1.0 / std::sqrt(1.0 - v_inv_c * v_inv_c);
        // Коэффициент коррекции направления
        double w2 = 1.0 / (1.0 + v_dot_s_c);
        // Вспомогательный член для сохранения точности при малых v/c
        double f = rel_p / (1.0 + rel_p);

        // Получаем единичный аберрированный вектор источника (K_unit_aber в Fortran)
        Eigen::Vector3d k_unit_aber = w2 * (k_s / rel_p + v_total / cnst::C + f * v_dot_s_c * (v_total / cnst::C));
        k_unit_aber.normalize();

        // --- 2. ПЕРЕХОД В ТОПОЦЕНТРИЧЕСКУЮ СИСТЕМУ (VEN) ---
        // Шаг А: Из GCRS(J2000) в ITRF (Crust-Fixed). Используем r2000[0]
        Eigen::Matrix3d rt2000 = r2000[0].transpose();
        Eigen::Vector3d crust_star = rt2000 * k_unit_aber;
        
        // Шаг Б: Из ITRF в локальную систему VEN (Vertical, East, North)
        Eigen::Matrix3d vw_tr = vw[idx].transpose();
        Eigen::Vector3d ven_star = vw_tr * crust_star;

        // --- 3. ВЫЧИСЛЕНИЕ УГЛОВ ---
        // Возвышение (Elevation) = asin(Vertical_component)
        e(j, 0) = std::asin(ven_star(0));
        
        // Азимут (Azimuth). В Fortran: atan2(East, North)
        az(j, 0) = std::atan2(ven_star(1), ven_star(2));
        if (az(j, 0) < 0.0) az(j, 0) += cnst::TWOPI;

        // --- 4. ПРОИЗВОДНЫЕ ПО ВРЕМЕНИ ---
        // Используем r2000[1] (первая производная матрицы вращения)
        if (r2000.size() > 1) {
            Eigen::Matrix3d rt2000_der = r2000[1].transpose();
            Eigen::Vector3d ven_star_der = vw_tr * (rt2000_der * k_unit_aber);

            // dE/dt = V_der / cos(E)
            e(j, 1) = ven_star_der(0) / std::cos(e(j, 0));
            
            // dAz/dt для atan2(y, x)
            double x = ven_star(2); // North
            double y = ven_star(1); // East
            az(j, 1) = (ven_star_der(1) * x - y * ven_star_der(2)) / (x * x + y * y);
        }
    }
}

} // namespace ariadna