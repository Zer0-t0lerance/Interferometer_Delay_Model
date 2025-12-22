// GEOD.cpp

#include "functions.h"
#include <algorithm>
#include <iostream>

namespace ariadna {

void GEOD(double equatorial_radius_r, double z_polar, double& geodetic_latitude_fi, double& geodetic_height_h) {
    using namespace cnst;
    using std::sqrt;
    using std::pow;
    using std::atan;
    using std::copysign;
    using std::acos;
    using std::cos;
    using std::fabs;
    using std::min;
    using std::max;

    // Константы IERS
    const double A = A_SEMI_MAJOR_AXIS;
    const double FR = FR_INV_FLATTENING;

    // Расчет полярной полуоси b (b = a * (1 - 1/fr)) со знаком Z
    double b_polar_semi_minor = copysign(A - A / FR, z_polar);
    
    // Переменные Борковского E и F
    // E = ((z + b)*b/a - a)/r
    double E = ((z_polar + b_polar_semi_minor) * b_polar_semi_minor / A - A) / equatorial_radius_r;
    // F = ((z - b)*b/a + a)/r
    double F = ((z_polar - b_polar_semi_minor) * b_polar_semi_minor / A + A) / equatorial_radius_r;

    // Параметры для решения t**4 + 2*E*t**3 + 2*F*t - 1 = 0
    double P = (E * F + 1.0) * 4.0 / 3.0;
    double Q = (E * E - F * F) * 2.0;
    double D = P * P * P + Q * Q;

    double t; // Корень t

    if (D >= 0.0) {
        // Случай с действительными корнями
        double s = sqrt(D);
        double y = copysign(pow(s, 2.0) / 3.0, pow(Q, 3.0));
        
        // Аргумент для acos: (P-y) / sqrt(y^2 + (y-P)^2)
        double cos_arg = (P - y) / sqrt(y * y + pow(y - P, 2.0));
        
        // Ограничение аргумента acos в диапазоне [-1, 1] для безопасности
        cos_arg = min(1.0, max(-1.0, cos_arg));

        t = sqrt(copysign(fabs((s + Q) / 2.0), 1.0) + (y - P) * cos(acos(cos_arg) / 3.0));
    } else {
        // Случай с комплексными корнями (три действительных)
        double R3_root = sqrt(-D);
        
        double P_cubed = P * P * P;
        
        double R2 = sqrt(fabs((P_cubed + R3_root) / 2.0));
        double R1 = sqrt(fabs((P_cubed - R3_root) / 2.0));
        
        R2 = copysign(R2, Q);
        R1 = copysign(R1, Q);
        
        double R3 = R1 + R2; 
        double s = copysign(pow(R3, 2.0) / 3.0, Q); // Q^3 sign simplified to Q sign for double
        double y = sqrt(fabs(pow(s, 2.0) / 3.0));

        t = copysign(sqrt(fabs((s + Q) / 2.0)), 1.0); 
        
        // Аргумент для acos: (P-s) / sqrt(y^2 + (s-P)^2)
        double cos_arg = (P - s) / sqrt(y * y + pow(s - P, 2.0));
        
        // Ограничение аргумента acos в диапазоне [-1, 1] для безопасности
        cos_arg = min(1.0, max(-1.0, cos_arg));
        
        t = t + y * cos(acos(cos_arg) / 3.0);
    }
    
    // Лямбда-функция для выполнения одного шага итерации (Ньютон-Рафсон)
    // t = t - (E*t*t + F*t - E*F) / (2.d0*t*t + 2.d0*E*t + F - F)
    auto newton_step = [&](double t_old) {
        double numerator = E * t_old * t_old + F * t_old - E * F;
        // Знаменатель: 2.0 * t_old * (t_old + E)
        double denominator = 2.0 * t_old * t_old + 2.0 * E * t_old;
        
        if (fabs(denominator) < 1e-15) { // Проверка на ноль (для избежания краха)
            return t_old; 
        }
        return t_old - numerator / denominator;
    };

    // Три итерации уточнения (согласно Fortran-коду)
    t = newton_step(t);
    t = newton_step(t);
    t = newton_step(t);

    // Расчет геодезической широты (fi) и высоты (h)
    // num_fi_h = r*t + z*F*F/a*E/E
    double num_fi_h = equatorial_radius_r * t + z_polar * F * F / A; // E/E = 1, опущено
    // den_fi_h = z*t - r*F*F/a*E/E
    double den_fi_h = z_polar * t - equatorial_radius_r * F * F / A; // E/E = 1, опущено

    // Геодезическая широта [rad]
    geodetic_latitude_fi = atan(num_fi_h / den_fi_h);
    
    // Геодезическая высота [m]
    geodetic_height_h = sqrt(pow(num_fi_h, 2.0) + pow(den_fi_h, 2.0)) - A;

    // Коррекция знака широты (if (z .lt. 0.d0) fi = -fi)
    if (z_polar < 0.0) {
        geodetic_latitude_fi = -geodetic_latitude_fi; 
    }
}

} // namespace ariadna