#include "functions.h"
#include "constants.h"
#include <cmath>

namespace ariadna {

void GEOD(double r, double z, double& fi, double& h) {
    // Константы эллипсоида WGS84 / IERS
    double a = cnst::AE; 
    double fr = cnst::F; 

    // Защита от точки строго на оси Z (полюс), чтобы избежать деления на ноль при расчете E и F
    if (r < 1e-10) {
        fi = (z > 0.0) ? cnst::HALFPI : -cnst::HALFPI;
        h = std::abs(z) - (a - a / fr);
        return;
    }

    // Аналог фортрановского dsign(a - a/fr, z)
    double b = std::copysign(a - a / fr, z);
    
    double E = ((z + b) * b / a - a) / r;
    double F_val = ((z - b) * b / a + a) / r;

    double P = (E * F_val + 1.0) * 4.0 / 3.0;
    double Q = (E * E - F_val * F_val) * 2.0;

    double D = P * P * P + Q * Q;

    double v, s;
    if (D >= 0.0) {
        s = std::sqrt(D) + Q;
        // Аналог фортрановского dsign(dexp(dlog(dabs(s))/3.d0), s)
        s = std::copysign(std::exp(std::log(std::abs(s)) / 3.0), s);
        v = P / s - s;
        // Строго по Фортрану и твоему коду (деление на 3.0 * P)
        v = -(Q + Q + v * v * v) / (3.0 * P);
    } else {
        // Исправлена опечатка cos*acos
        v = 2.0 * std::sqrt(-P) * std::cos(std::acos(Q / P / std::sqrt(-P)) / 3.0);
    }

    double G = 0.5 * (E + std::sqrt(E * E + v));
    double t = std::sqrt(G * G + (F_val - v * G) / (G + G - E)) - G;
    
    // Строго по Фортрану и твоему коду
    fi = std::atan((1.0 - t * t) * a / (2.0 * b * t));
    h = (r - a * t) * std::cos(fi) + (z - b) * std::sin(fi);
}

} // namespace ariadna