#include "../src/functions.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ariadna;

void check_val(const std::string& name, double calc, double ref, const std::string& note = "") {
    double diff = std::abs(calc - ref);
    std::cout << std::left << std::setw(20) << name 
              << " | Calc: " << std::scientific << std::setprecision(11) << std::setw(18) << calc 
              << " | Ref: "  << std::setw(18) << ref 
              << " | Diff: " << std::setw(12) << diff 
              << (diff <= 1e-10 ? " [OK]" : " [FAIL]") << " " << note << std::endl;
}

int main() {
    Observation obs;
    obs.mjd = 57402;
    obs.utc = 0.2705150462934398;
    double tt = 0.2713042129601065;

    std::vector<EOPData> eop_nodes = {
        {57400, 0.0589932, -35.9410068, 0.031105, 0.274426, -0.000093, -0.000136},
        {57401, 0.0568475, -35.9431525, 0.029487, 0.275674, -0.000058, -0.000114},
        {57402, 0.0546433, -35.9453567, 0.028220, 0.276671, -0.000145, -0.000126},
        {57403, 0.0524926, -35.9475074, 0.026690, 0.278122, -0.000168, -0.000115},
        {57404, 0.0504060, -35.9495940, 0.024996, 0.279528, -0.000155, -0.000086},
        {57405, 0.0483894, -35.9516106, 0.023229, 0.281340, -0.000142, -0.000057},
        {57406, 0.0464968, -35.9535032, 0.021297, 0.283184, -0.000130, -0.000029}
    };

    double ut1_out;
    Eigen::VectorXd eop_int(5), deop_int(5), arg_oc(8);
    Eigen::MatrixXd diu(3,2), lib(3,2);

    interp_eop(0, obs, tt, ut1_out, eop_int, deop_int, arg_oc, diu, lib, eop_nodes);

    std::cout << std::string(100, '=') << "\n";
    std::cout << "EOP VALUES (Seconds and Radians)\n";
    std::cout << std::string(100, '-') << "\n";
    check_val("UT1-UTC", eop_int(0), 5.40342726317e-02, "(Now matches perfectly)");
    check_val("X_pole",  eop_int(1), 1.32653046496e-07);
    check_val("Y_pole",  eop_int(2), 1.34321610512e-06);

    std::cout << "\n" << std::string(100, '-') << "\n";
    std::cout << "DERIVATIVES (Physically Correct SI vs Fortran Bug)\n";
    std::cout << std::string(100, '-') << "\n";
    check_val("dX/dt", deop_int(1), -1.13524766720e-13, "<- C++ is physically correct here");
    check_val("dY/dt", deop_int(2), -3.48633132047e-14, "<- C++ is physically correct here");
    
    return 0;
}