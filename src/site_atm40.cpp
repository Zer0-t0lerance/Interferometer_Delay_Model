// site_atm40.cpp
//
// Порт подпрограммы SITE_ATM40 (FORTRAN_SOURCE/SITE_ATM40.FOR): смещение и
// скорость станции из-за атмосферной нагрузки в системе J2000.0.
//
// Модель (каталог VLBI_atmload4_12.cat): по каждой компоненте VEN
//   dr = Sum_{i=1..3}[A_i cos(omega_i t) + B_i sin(omega_i t)],
// для вертикали дополнительно + b0 + b1*(P - p_0);
//   t = MJD - ATM_EPOCH_MJD [сут], omega_i = TWOPI/ATM_PERIODS[i] [рад/сут].
// Амплитуды каталога — в мм; перевод в метры и скорость мм/сут -> м/с здесь же.
// Затем VEN -> ITRF (vw_i) -> J2000 (r2000), для скорости учитывается вращение.
//
// Частные производные (dr1_dA2000 и т.п.) из Фортрана НЕ вычисляются — они нужны
// только для слоя оценивания (МНК), который в задаче геометрических задержек не участвует.

#include "functions.h"

namespace ariadna {

void SITE_ATM40(int mjd, double utc, double pres, double dPdt,
                const AtmLoadData& atm,
                const Eigen::Matrix3d& vw_i,
                const Eigen::Matrix3d& r2000, const Eigen::Matrix3d& dr2000_dt,
                Eigen::Vector3d& dx_atm, Eigen::Vector3d& dv_atm)
{
    // Время от опорной эпохи ряда (1984-01-01), сутки.
    const double Dt = (static_cast<double>(mjd) - cnst::ATM_EPOCH_MJD) + utc;

    // Угловые скорости и аргументы трёх гармонических членов.
    double omega[3], cosA[3], sinA[3];
    for (int n = 0; n < 3; ++n) {
        omega[n] = cnst::TWOPI / cnst::ATM_PERIODS[n];    // рад/сут
        double arg = std::fmod(omega[n] * Dt, cnst::TWOPI);
        cosA[n] = std::cos(arg);
        sinA[n] = std::sin(arg);
    }

    // Топоцентрическое смещение (dr) [м] и скорость (dv) [м/с] по компонентам VEN.
    Eigen::Vector3d dr_ven = Eigen::Vector3d::Zero();
    Eigen::Vector3d dv_ven = Eigen::Vector3d::Zero();

    for (int c = 0; c < 3; ++c) {              // c = 0(Vertical), 1(East), 2(North)
        double dr = 0.0; // мм
        double dv = 0.0; // мм/сут
        for (int n = 0; n < 3; ++n) {
            const double A = atm.coef(c, 2 * n);      // A_{n+1}
            const double B = atm.coef(c, 2 * n + 1);  // B_{n+1}
            dr += A * cosA[n] + B * sinA[n];
            dv += -A * omega[n] * sinA[n] + B * omega[n] * cosA[n];
        }
        // Вертикаль: постоянный член b0 и регрессия по локальному давлению b1*(P - p_0).
        if (c == 0) {
            const double b0 = atm.coef(0, 6);
            const double b1 = atm.coef(0, 7);
            dr += b0 + b1 * (pres - atm.p_0);
            dv += b1 * dPdt;
        }
        dr_ven(c) = dr * 1.0e-3;                       // мм -> м
        dv_ven(c) = dv / cnst::SECDAY * 1.0e-3;        // мм/сут -> м/с
    }

    // VEN -> ITRF -> J2000 (та же конвенция, что в POLE_TIDE / SITE_TIDE_OC).
    Eigen::Vector3d dx_itrf = vw_i * dr_ven;
    Eigen::Vector3d dv_itrf = vw_i * dv_ven;

    dx_atm = r2000 * dx_itrf;
    dv_atm = r2000 * dv_itrf + dr2000_dt * dx_itrf;
}

} // namespace ariadna
