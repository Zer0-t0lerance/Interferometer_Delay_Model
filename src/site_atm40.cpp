// SITE_ATM40.cpp (Обновленная версия)

#include "functions.h"

namespace ariadna {

void SITE_ATM40(double dPdt,
                const Eigen::Matrix3d& vw_i, const Eigen::Matrix3d& r2000,
                Eigen::Vector3d& dx_atm, Eigen::Vector3d& dv_atm) {

    using namespace cnst;
    using Eigen::Vector3d;
    using Eigen::Matrix3d;

    // Скорость изменения давления (dP/dt) [Pa/s]
    double dPdt_Pa = dPdt * ATM_CONV_MBAR_TO_PA;

    // Коэффициенты атмосферной нагрузки K_atm [m/Pa] (Up, North, East)
    const Vector3d K_atm_VEN(ATM_LOADING_COEFF_ARRAY[0], ATM_LOADING_COEFF_ARRAY[1], ATM_LOADING_COEFF_ARRAY[2]);

    // A. Смещение dr_ven (Up, North, East) - заглушка
    Vector3d dr_ven = Vector3d::Zero();

    // B. Скорость dv_ven (Up, North, East)
    Vector3d dv_ven = K_atm_VEN * dPdt_Pa;

    // Преобразование из VEN в ITRF
    Vector3d dx_itrf = vw_i * dr_ven;
    Vector3d dv_itrf = vw_i * dv_ven;

    // Преобразование в J2000.0
    dx_atm = r2000 * dx_itrf;
    dv_atm = r2000 * dv_itrf; // Упрощение, без производной
}

} // namespace ariadna