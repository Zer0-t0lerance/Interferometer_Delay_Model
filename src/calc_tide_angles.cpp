// calc_tide_angles.cpp (расчет угловых аргументов приливов)
#include "functions.h"

namespace ariadna {

// Вспомогательная функция для расчета угловых аргументов приливов
void calc_tide_angles(int mjd_start, double ut1_sec, Eigen::VectorXd& angle, Eigen::VectorXd& speed_angle) {
    using namespace cnst;
    using std::sin;
    using std::cos;
    using std::fmod;

    // Векторы для угла и скорости
    angle.resize(NUM_TIDES);
    speed_angle.resize(NUM_TIDES);

    // Эпоха в юлианских столетиях с J2000.0
    double MJD_J2000 = JD2000 - 2400000.5; // MJD2000 = 51544.5
    double CAPT = (mjd_start - MJD_J2000) / JUL_CENT; // T_0 [Julian centuries]

    // Mean Longitude of the Moon (H0) at beginning of day [rad]
    double H0 = ((( -1.135e-6 * CAPT + 0.000305e0 ) * CAPT + 481267.88177531e0)
               * CAPT + 218.316479e0 ) * CDEGRAD;
    
    // Mean Longitude of the Sun (S0) at beginning of day [rad]
    double S0 = ((( 1.9e-6 * CAPT - 0.001133e0 ) * CAPT + 481267.88314137e0)
               * CAPT + 270.434358e0 ) * CDEGRAD;

    // Mean Longitude of Lunar Perigee (P0) at beginning of day [rad]
    double P0 = ((( -1.2e-5 * CAPT - .010325e0 ) * CAPT + 4069.0340329577e0)
               * CAPT + 334.329653e0 ) * CDEGRAD;

    // Расчет угловых аргументов и их скоростей для 11 волн
    for (int k = 0; k < NUM_TIDES; ++k) {
        // Угловой аргумент для времени 0 UT1
        double ang_k = ANGFAC[0][k] * H0 + ANGFAC[1][k] * S0 
                     + ANGFAC[2][k] * P0 + ANGFAC[3][k] * TWOPI;
        
        // Угловой аргумент в UT1 (sec)
        angle(k) = SPEED_TIDES[k] * ut1_sec + ang_k;

        // Нормализация к [0, 2*PI)
        angle(k) = fmod(angle(k), TWOPI);
        if (angle(k) < 0.0) {
            angle(k) += TWOPI;
        }

        // Скорость приливного аргумента [rad/s]
        speed_angle(k) = SPEED_TIDES[k];
    }
}

} // namespace ariadna