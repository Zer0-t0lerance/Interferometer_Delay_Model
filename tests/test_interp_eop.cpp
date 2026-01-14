#include "../src/functions.h"
#include <iostream>
#include <iomanip>

using namespace ariadna;

void print_row(const std::string& name, double calc, double exp, double tol) {
    double diff = std::abs(calc - exp);
    std::cout << std::left << std::setw(30) << name << " | "
              << std::scientific << std::setprecision(14) << calc << " | "
              << exp << " | "
              << (diff < tol ? "OK" : "FAIL") << " (diff: " << diff << ")" << std::endl;
}

int main() {
    Observation obs;
    obs.mjd = 57402;
    obs.utc = 0.2705150462934398;
    double tt = 0.2713042129601065;

    // Те самые данные из твоего блока (Fortran Reference)
    double f_arg0    = 5.5548599903828810e-01;
    double f_argdot0 = 7.2921158553754072e-05;
    double f_dut     = 1.4297097582344920e-01;
    double f_dx      = 7.6467982137447488e-01;
    double f_dy      = -5.5958508447481958e-01;
    double f_ut1_tai = -3.5934987846776899e+01;

    // Входные данные (узлы)
    std::vector<EOPData> eop_nodes = {
        {57400, 0.0589932, -35.9410068, 0.031105, 0.274426, -0.00009, -0.00014},
        {57401, 0.0568475, -35.9431525, 0.029487, 0.275674, -0.00006, -0.00011},
        {57402, 0.0546433, -35.9453567, 0.028220, 0.276671, -0.00015, -0.00013},
        {57403, 0.0524926, -35.9475074, 0.026690, 0.278122, -0.00017, -0.00012},
        {57404, 0.0504060, -35.9495940, 0.024996, 0.279528, -0.00015, -0.00009},
        {57405, 0.0483894, -35.9516106, 0.023229, 0.281340, -0.00014, -0.00006},
        {57406, 0.0464968, -35.9535032, 0.021297, 0.283184, -0.00013, -0.00003}
    };

    double ut1_out;
    Eigen::VectorXd eop_int(5), deop_int(5), arg_oc(8);
    Eigen::MatrixXd diu(3,2), lib(3,2);

    // Запуск
    interp_eop(0, obs, tt, ut1_out, eop_int, deop_int, arg_oc, diu, lib, eop_nodes);

    std::cout << std::string(110, '=') << "\n";
    std::cout << std::left << std::setw(30) << "Parameter" << " | " << std::setw(22) << "Calculated" << " | " << std::setw(22) << "Fortran Ref" << " | Status\n";
    std::cout << std::string(110, '-') << "\n";

    // 1. Проверка аргумента (Самое важное!)
    print_row("Arg[0] (Theta)", arg_oc(0), f_arg0, 1e-12);
    
    // 2. Проверка поправок (если аргумент врет, это тоже будет врать)
    print_row("DUT Tide (sec)", diu(0,0), f_dut, 1e-10);
    print_row("DX Tide (arcsec)", diu(1,0), f_dx, 1e-10);
    print_row("DY Tide (arcsec)", diu(2,0), f_dy, 1e-10);

    // 3. Итоговый результат
    print_row("UT1-TAI (Final)", eop_int(0), f_ut1_tai, 1e-11);

    return 0;
}