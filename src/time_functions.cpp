#include "functions.h"
#include "constants.h"
#include <cmath>

namespace ariadna {

// Функция для расчета leap seconds (NSEC)
void nsec(double mjd, int& idelt) {
    // Используем таблицу leap seconds из constants.h
    idelt = cnst::LEAP_SECONDS_BASE;
    for (size_t i = 0; i < sizeof(cnst::LEAP_SECOND_MJD_INTERVALS) / sizeof(double); ++i) {
        if (mjd >= cnst::LEAP_SECOND_MJD_INTERVALS[i]) {
            idelt++;
        } else {
            break;
        }
    }
}

// Функция для расчета TAI и TT
void tai_time(double mjd, double UTC, double& TAI, double& TT) {
    int idelt;
    nsec(mjd, idelt);
    TAI = UTC + static_cast<double>(idelt) / cnst::SECDAY;
    TT = TAI + 32.184 / cnst::SECDAY;
}

// Функция для расчета фундаментальных аргументов
void fund_arg(double jd, double ct, double& cent, Eigen::VectorXd& f, Eigen::VectorXd& fd) {
    cent = (jd - cnst::JD2000) / cnst::JUL_CENT;
    double t = cent;

    // Фундаментальные аргументы (в угловых секундах)
    f(0) = cnst::F1[0] + cnst::F1[1] * t + cnst::F1[2] * t * t + cnst::F1[3] * t * t * t + cnst::F1[4] * t * t * t * t; // l
    f(1) = cnst::F2[0] + cnst::F2[1] * t + cnst::F2[2] * t * t + cnst::F2[3] * t * t * t + cnst::F2[4] * t * t * t * t; // l'
    f(2) = cnst::F3[0] + cnst::F3[1] * t + cnst::F3[2] * t * t + cnst::F3[3] * t * t * t + cnst::F3[4] * t * t * t * t; // F
    f(3) = cnst::F4[0] + cnst::F4[1] * t + cnst::F4[2] * t * t + cnst::F4[3] * t * t * t + cnst::F4[4] * t * t * t * t; // D
    f(4) = cnst::F5[0] + cnst::F5[1] * t + cnst::F5[2] * t * t + cnst::F5[3] * t * t * t + cnst::F5[4] * t * t * t * t; // Omega

    // Производные (в угловых секундах за век)
    fd(0) = cnst::F1[1] + 2 * cnst::F1[2] * t + 3 * cnst::F1[3] * t * t + 4 * cnst::F1[4] * t * t * t;
    fd(1) = cnst::F2[1] + 2 * cnst::F2[2] * t + 3 * cnst::F2[3] * t * t + 4 * cnst::F2[4] * t * t * t;
    fd(2) = cnst::F3[1] + 2 * cnst::F3[2] * t + 3 * cnst::F3[3] * t * t + 4 * cnst::F3[4] * t * t * t;
    fd(3) = cnst::F4[1] + 2 * cnst::F4[2] * t + 3 * cnst::F4[3] * t * t + 4 * cnst::F4[4] * t * t * t;
    fd(4) = cnst::F5[1] + 2 * cnst::F5[2] * t + 3 * cnst::F5[3] * t * t + 4 * cnst::F5[4] * t * t * t;
}

// Реализация t_eph
void t_eph(const Observation& obs, double tai, double ut1, double tt, double lon_gcen, double u_site, double v_site, double& ct, double& dtaidct) {
    // Relativistic correction for coordinate time
    // CT = TT + (u * lambda + v) / c^2, where u, v in km, lambda in rad
    double correction = (u_site * 1000.0 * lon_gcen + v_site * 1000.0) / (cnst::C * cnst::C);
    ct = tt + correction;
    dtaidct = 1.0; // Approximation
}

// Функция для расчета приливных поправок EOP
void terms_71(double cent, const Eigen::VectorXd& f, const Eigen::VectorXd& fd, Eigen::MatrixXd& dEOP_diu, Eigen::VectorXd& arg_oc_tide) {
    double t = cent;
    double t2 = t * t;
    double t3 = t * t2;
    double dcent = 1.0 / (cnst::JUL_CENT * cnst::SECDAY);

    // Расчет Theta (GMST)
    double theta = (cnst::THETA_COEFFS[4] * t3 + cnst::THETA_COEFFS[3] * t2 +
                    cnst::THETA_COEFFS[2] * t + cnst::THETA_COEFFS[1] * t + cnst::THETA_COEFFS[0]) * 15.0 + cnst::SEC360 / 2.0;
    double dtheta = ((3.0 * cnst::THETA_COEFFS[4] * t2 + 2.0 * cnst::THETA_COEFFS[3] * t) * dcent +
                     cnst::THETA_COEFFS[2] * dcent + cnst::THETA_COEFFS[1] * dcent) * 15.0;

    double dut = 0.0, dx = 0.0, dy = 0.0, dutdt = 0.0, dxdt = 0.0, dydt = 0.0;
    int k = 0;

    for (int j = 0; j < 71; ++j) {  // 71 terms
        double arg = 0.0, argdot = 0.0;
        for (int i = 0; i < 5; ++i) {
            arg += cnst::amp_XY[j][i] * f(i);
            argdot += cnst::amp_XY[j][i] * fd(i) * dcent;
        }
        arg += cnst::amp_XY[j][5] * theta;
        argdot += cnst::amp_XY[j][5] * dtheta;

        arg = fmod(arg, cnst::SEC360) * cnst::CARCRAD;
        argdot *= cnst::CARCRAD;

        // Сохраняем аргументы для определенных членов
        if (j == 6 || j == 11 || j == 22 || j == 26 || j == 48 || j == 55 || j == 62 || j == 65) {
            arg_oc_tide(k++) = arg;
        }

        dut += cnst::amp_UT[j][6] * sin(arg) + cnst::amp_UT[j][7] * cos(arg);
        dx += cnst::amp_XY[j][6] * sin(arg) + cnst::amp_XY[j][7] * cos(arg);
        dy += cnst::amp_XY[j][8] * sin(arg) + cnst::amp_XY[j][9] * cos(arg);

        dutdt += (cnst::amp_UT[j][6] * cos(arg) - cnst::amp_UT[j][7] * sin(arg)) * argdot;
        dxdt += (cnst::amp_XY[j][6] * cos(arg) - cnst::amp_XY[j][7] * sin(arg)) * argdot;
        dydt += (cnst::amp_XY[j][8] * cos(arg) - cnst::amp_XY[j][9] * sin(arg)) * argdot;
    }

    dEOP_diu(0, 0) = dut * 1e-6;
    dEOP_diu(1, 0) = dx * 1e-6;
    dEOP_diu(2, 0) = dy * 1e-6;
    dEOP_diu(0, 1) = dutdt * 1e-6;
    dEOP_diu(1, 1) = dxdt * 1e-6;
    dEOP_diu(2, 1) = dydt * 1e-6;
}

} // namespace ariadna