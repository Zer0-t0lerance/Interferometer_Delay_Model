#include "functions.h"

namespace ariadna {

void SITE_TIDE_OC(double jd, double ut1_sec, const OceanTideData& tide_data,
                  const Eigen::Matrix3d& vw_i, 
                  const Eigen::Matrix3d& r2000, 
                  const Eigen::Matrix3d& dr2000_dt,
                  Eigen::Vector3d& dx_octide, Eigen::Vector3d& dv_octide) 
{
    // --- 1. Вычисление фундаментальных аргументов (OC_ARG) ---
    double icapd = jd - 2442412.5; // Дни от 0 Jan 1975
    double capt = (27392.500528 + 1.000000035 * icapd) / cnst::JUL_CENT;

    double h0 = (279.69668 + (36000.768930485 + 3.03e-4 * capt) * capt) * cnst::CDEGRAD;
    double s0 = (((1.9e-6 * capt - 0.001133) * capt + 481267.88314137) * capt + 270.434358) * cnst::CDEGRAD;
    double p0 = (((-1.2e-5 * capt - 0.010325) * capt + 4069.0340329577) * capt + 334.329653) * cnst::CDEGRAD;

    double angle[11];
    for (int j = 0; j < 11; ++j) {
        angle[j] = cnst::SPEED_TIDES[j] * ut1_sec + 
                   cnst::ANGFAC[j][0] * h0 + 
                   cnst::ANGFAC[j][1] * s0 + 
                   cnst::ANGFAC[j][2] * p0 + 
                   cnst::ANGFAC[j][3] * cnst::TWOPI;
                   
        angle[j] = std::fmod(angle[j], cnst::TWOPI);
        if (angle[j] < 0.0) angle[j] += cnst::TWOPI;
    }

    // --- 2. Вычисление топоцентрических смещений (Vertical, West, South) ---
    Eigen::Vector3d drsw   = Eigen::Vector3d::Zero();
    Eigen::Vector3d drsw_v = Eigen::Vector3d::Zero();

    // k = 0 (Up), 1 (West), 2 (South)
    for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < 11; ++j) {
            // ИЗМЕНЕНО: теперь читаем (Ось, Волна) -> (k, j)
            double amp = tide_data.amplitudes(k, j);
            double phs = tide_data.phases(k, j) * cnst::CDEGRAD;

            drsw(k)   += amp * std::cos(angle[j] - phs);
            drsw_v(k) -= amp * std::sin(angle[j] - phs) * cnst::SPEED_TIDES[j];
        }
    }

    // --- 3. Переход от (Up, West, South) к (Up, East, North) - VEN система ---
    drsw(1)   = -drsw(1);
    drsw(2)   = -drsw(2);
    drsw_v(1) = -drsw_v(1);
    drsw_v(2) = -drsw_v(2);

    // --- 4. Преобразование VEN -> ITRF -> J2000 ---
    // Координаты: X_J2000 = R2000 * (VW * dX_VEN)
    Eigen::Vector3d work = vw_i * drsw;
    dx_octide = r2000 * work;

    // Скорости: Производная сложной функции
    Eigen::Vector3d dwork = vw_i * drsw_v;
    Eigen::Vector3d v1 = r2000 * dwork;
    Eigen::Vector3d v2 = dr2000_dt * work;
    dv_octide = v1 + v2;
}

} // namespace ariadna