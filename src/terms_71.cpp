#include "functions.h"

namespace ariadna {
void terms_71(double cent, const Eigen::VectorXd& f, const Eigen::VectorXd& fd, Eigen::MatrixXd& dEOP_diu, Eigen::VectorXd& arg_oc_tide) {

// Инициализация сумм
    double dut = 0.0, dx = 0.0, dy = 0.0;
    double dutdt = 0.0, dxdt = 0.0, dydt = 0.0;
    
    double arg = 0.0;
    double argdot = 0.0;

// Вычисление Theta и dTheta (GMST-like)
    double t = cent;
    double t2 = t * t;
    double t3 = t * t2;
//    double t4 = t * t3;

    double dcent = 1.0 / (cnst::JUL_CENT * cnst::SECDAY); // Век в дни

//    printf("dcent: %.16le\n", dcent);

    double theta = (cnst::THETA_COEFFS[4] * t3 + cnst::THETA_COEFFS[3] * t2 + cnst::THETA_COEFFS[2] * t + cnst::THETA_COEFFS[1] * t + cnst::THETA_COEFFS[0]) * 15.0 + cnst::SEC360 / 2.0;
    double dtheta = (( 3.0 * cnst::THETA_COEFFS[4] * t2 + 2.0 * cnst::THETA_COEFFS[3] * t ) * dcent + cnst::THETA_COEFFS[2] * dcent + cnst::THETA_COEFFS[1] * dcent ) * 15.0;

//    printf("theta (arcsec): %.16le\n", theta);
//    printf("dtheta (arcsec): %.16le\n", dtheta);

    int k = 0; // Счетчик для arg_oc_tide

    constexpr int n = 71;
    for (int j = 0; j < n; j++)
    {
        arg = 0.0;
        argdot = 0.0;

        if (j == 6 || j == 11 || j == 22 || j == 26 || j == 48 || j == 55 || j == 62 || j == 65)
        {
            for (int i = 0; i < 5; i++)
            {
//              Argument of tidal term(Arcsecond)
                arg += cnst::amp_XY[j][i] * f[i];
//              Derivatives of the arguments(Arcsec / Sec)
                argdot += cnst::amp_XY[j][i] * fd[i] * dcent;
            }
        }
        
        arg += cnst::amp_XY[j][5] * theta;
        argdot += cnst::amp_XY[j][5] * dtheta;

//        printf("arg = %.16le\n", arg);
//        printf("fmod(arg) = %.16le\n", fmod(arg, cnst::SEC360));
//        Conversion of argument to radians
        arg = fmod(arg, cnst::SEC360) * cnst::CARCRAD;
        argdot = argdot * cnst::CARCRAD;

//        printf("arg_rad = %.16le\n", arg);

        if (j == 6 || j == 11 || j == 22 || j == 26 || j == 48 || j == 55 || j == 62 || j == 65)
        {
//            Copy and store the arguments for calculation of partials (in DER_UT1, DER_POLAR)
            arg_oc_tide[k] = arg;
            k++;
        }

//      Total tides effect
        dut += cnst::amp_UT[j][6] * sin(arg) + cnst::amp_UT[j][7] * cos(arg);
        dx += cnst::amp_XY[j][6] * sin(arg) + cnst::amp_XY[j][7] * cos(arg);
        dy += cnst::amp_XY[j][8] * sin(arg) + cnst::amp_XY[j][9] * cos(arg);

//      The time derivatives of the EOP tidal corrections
        dutdt += (cnst::amp_UT[j][6] * cos(arg) - cnst::amp_UT[j][7] * sin(arg)) * argdot;
        dxdt += (cnst::amp_XY[j][6] * cos(arg) - cnst::amp_XY[j][7] * sin(arg)) * argdot;
        dydt += (cnst::amp_XY[j][8] * cos(arg) - cnst::amp_XY[j][9] * sin(arg)) * argdot;

        if(j == 0)
        {
            printf("j, arg, argdot \ndut, dx, dy \ndutdt, dxdt, dydt\n%d %.16le %.16le\n%.16le %.16le %.16le\n%.16le %.16le %.16le\n\n", j, arg, argdot, dut, dx, dy, dutdt, dxdt, dydt);
        }
    }

    // Конвертация единиц
    dEOP_diu(0, 0) = dut * 1e-6;    // UT1 sec
    dEOP_diu(1, 0) = dx * 1e-6;     // X arcsec
    dEOP_diu(2, 0) = dy * 1e-6;     // Y arcsec
    dEOP_diu(0, 1) = dutdt * 1e-6;  // dUT1/dt sec/sec
    dEOP_diu(1, 1) = dxdt * 1e-6;   // dX/dt arcsec/sec
    dEOP_diu(2, 1) = dydt * 1e-6;   // dY/dt arcsec/sec
}
} // namespace ariadna