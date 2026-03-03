#include <iostream>
#include <iomanip>
#include <cmath>
#include "../src/functions.h"

using namespace ariadna;

void check_val(const std::string& name, double calc, double ref, double tol = 1e-12) {
    double diff = std::abs(calc - ref);
    std::cout << std::left << std::setw(25) << name 
              << " | Calc: " << std::scientific << std::setprecision(14) << std::setw(22) << calc 
              << " | Ref: "  << std::setw(22) << ref 
              << " | Diff: " << std::setw(12) << diff 
              << (diff <= tol ? " [OK]" : " [FAIL]") << std::endl;
}

void get_args_from_cent(double cent, Eigen::VectorXd& f, Eigen::VectorXd& fd) {
    double jd_plus_ct = cent * cnst::JUL_CENT + cnst::JD2000;
    double dummy_cent;
    fund_arg(jd_plus_ct, 0.0, dummy_cent, f, fd);
}

void run_test(double cent, const Eigen::MatrixXd& ref_eop) {
    Eigen::VectorXd f(5), fd(5);
    get_args_from_cent(cent, f, fd);

    Eigen::MatrixXd calc_eop(3, 2);
    terms_lib(cent, f, fd, calc_eop);

    std::cout << std::string(95, '-') << "\n";
    std::cout << "TEST FOR CENT = " << std::fixed << std::setprecision(16) << cent << "\n";
    std::cout << std::string(95, '-') << "\n";

    check_val("dEOP_lib(UT1, val)", calc_eop(0,0), ref_eop(0,0));
    check_val("dEOP_lib(X,   val)", calc_eop(1,0), ref_eop(1,0));
    check_val("dEOP_lib(Y,   val)", calc_eop(2,0), ref_eop(2,0));
    
    check_val("dEOP_lib(UT1, der)", calc_eop(0,1), ref_eop(0,1), 1e-18);
    check_val("dEOP_lib(X,   der)", calc_eop(1,1), ref_eop(1,1), 1e-18);
    check_val("dEOP_lib(Y,   der)", calc_eop(2,1), ref_eop(2,1), 1e-18);
}

int main() {
    std::cout << std::string(95, '=') << "\n";
    std::cout << "VERIFICATION: TERMS_LIB (Librations of non-rigid Earth)\n";

    // =========================================================================
    // ТЕСТ 1: Первый дамп (cent = 0.1603770377607929)
    // =========================================================================
    double cent1 = 0.1603770377607929;
    Eigen::MatrixXd ref_eop1(3, 2);
    ref_eop1 <<  0.1348781963026778e-06,  0.6400216191562147e-10,
                -0.3514768458611694e-05,  0.6224462386135502e-09,
                 0.7596889981598940e-05,  0.2548981056098484e-09;
    
    run_test(cent1, ref_eop1);

    // =========================================================================
    // ТЕСТ 2: Второй дамп (cent = 0.1603633667275078)
    // =========================================================================
    double cent2 = 0.1603633667275078;
    Eigen::MatrixXd ref_eop2(3, 2);
    ref_eop2 << -0.2617300740670654e-06, -0.6331574525833121e-10,
                 0.6460839537146565e-05, -0.6486979494179335e-09,
                -0.7993582907186960e-05, -0.4532551128505596e-09;

    run_test(cent2, ref_eop2);

    std::cout << std::string(95, '=') << "\n";
    return 0;
}