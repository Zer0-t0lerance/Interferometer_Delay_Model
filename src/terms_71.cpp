#include "functions.h"

namespace ariadna {

void terms_71(double cent, const Eigen::VectorXd& f, const Eigen::VectorXd& fd, Eigen::MatrixXd& dEOP_diu, Eigen::VectorXd& arg_oc_tide) {
    dEOP_diu.setZero(3, 2);
    arg_oc_tide.setZero(8);

    double dut = 0.0, dx = 0.0, dy = 0.0;
    double dutdt = 0.0, dxdt = 0.0, dydt = 0.0;

    double t = cent;
    double t2 = t * t;
    double t3 = t2 * t;

    // Константа перевода из столетий в секунды
    const double dcent = 1.0 / (cnst::JUL_CENT * cnst::SECDAY);

    // 1. Расчет Theta (арксекунды)
    // Формула из Фортрана: (th1 + (th2+th3)*t + th4*t^2 + th5*t^3) * 15 + 648000
    double theta = (cnst::THETA_COEFFS[0] + 
                    (cnst::THETA_COEFFS[1] + cnst::THETA_COEFFS[2]) * t + 
                    cnst::THETA_COEFFS[3] * t2 + 
                    cnst::THETA_COEFFS[4] * t3) * 15.0 + 648000.0;

    // 2. Расчет dTheta (арксекунды/сек)
    // Производная: ( (th2+th3) + 2*th4*t + 3*th5*t^2 ) * 15 * dcent
    double dtheta = ((cnst::THETA_COEFFS[1] + cnst::THETA_COEFFS[2]) + 
                     2.0 * cnst::THETA_COEFFS[3] * t + 
                     3.0 * cnst::THETA_COEFFS[4] * t2) * 15.0 * dcent;

    int k = 0;
    for (int j = 0; j < 71; ++j) {
        double arg_as = 0.0;
        double argdot_as_sec = 0.0;

        for (int i = 0; i < 5; ++i) {
            arg_as += cnst::amp_UT[j][i] * f(i);
            // fd в arcsec/century -> переводим в arcsec/sec
            argdot_as_sec += cnst::amp_UT[j][i] * fd(i) * dcent;
        }

        arg_as += cnst::amp_UT[j][5] * theta;
        argdot_as_sec += cnst::amp_UT[j][5] * dtheta;

        // Нормализация аргумента и перевод в радианы
        double arg = std::fmod(arg_as, 1296000.0) * cnst::CARCRAD;
        if (arg < 0.0) arg += cnst::TWOPI;

        // Производная аргумента в рад/сек
        double argdot = argdot_as_sec * cnst::CARCRAD;

        // Сохранение выбранных аргументов (как в Фортране)
        if (j == 6 || j == 11 || j == 22 || j == 26 || j == 48 || j == 55 || j == 62 || j == 65) {
            if (k < 8) arg_oc_tide[k++] = arg;
        }

        double s = std::sin(arg);
        double c = std::cos(arg);

        // Суммирование эффектов (амплитуды в микро-единицах)
        dut += cnst::amp_UT[j][6] * s + cnst::amp_UT[j][7] * c;
        dx  += cnst::amp_XY[j][6] * s + cnst::amp_XY[j][7] * c;
        dy  += cnst::amp_XY[j][8] * s + cnst::amp_XY[j][9] * c;

        // Производные
        dutdt += (cnst::amp_UT[j][6] * c - cnst::amp_UT[j][7] * s) * argdot;
        dxdt  += (cnst::amp_XY[j][6] * c - cnst::amp_XY[j][7] * s) * argdot;
        dydt  += (cnst::amp_XY[j][8] * c - cnst::amp_XY[j][9] * s) * argdot;
    }

    // Перевод из микро-единиц в СИ (секунды и арксекунды)
    dEOP_diu(0, 0) = dut * 1e-6;
    dEOP_diu(1, 0) = dx * 1e-6;
    dEOP_diu(2, 0) = dy * 1e-6;

    dEOP_diu(0, 1) = dutdt * 1e-6;
    dEOP_diu(1, 1) = dxdt * 1e-6;
    dEOP_diu(2, 1) = dydt * 1e-6;
}

} // namespace ariadna