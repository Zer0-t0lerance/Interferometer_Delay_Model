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

// Вспомогательная функция: зная cent (столетия), можно найти JD и вызвать fund_arg
void get_args_from_cent(double cent, Eigen::VectorXd& f, Eigen::VectorXd& fd) {
    double jd_plus_ct = cent * cnst::JUL_CENT + cnst::JD2000;
    double dummy_cent;
    fund_arg(jd_plus_ct, 0.0, dummy_cent, f, fd);
}

void run_test(double cent, const Eigen::MatrixXd& ref_eop, const Eigen::VectorXd& ref_arg) {
    Eigen::VectorXd f(5), fd(5);
    get_args_from_cent(cent, f, fd);

    Eigen::MatrixXd calc_eop(3, 2);
    Eigen::VectorXd calc_arg(8);

    terms_71(cent, f, fd, calc_eop, calc_arg);

    std::cout << std::string(95, '-') << "\n";
    std::cout << "TEST FOR CENT = " << std::fixed << std::setprecision(16) << cent << "\n";
    std::cout << std::string(95, '-') << "\n";

    check_val("dEOP_diu(UT1, val)", calc_eop(0,0), ref_eop(0,0));
    check_val("dEOP_diu(X,   val)", calc_eop(1,0), ref_eop(1,0));
    check_val("dEOP_diu(Y,   val)", calc_eop(2,0), ref_eop(2,0));
    
    check_val("dEOP_diu(UT1, der)", calc_eop(0,1), ref_eop(0,1), 1e-18);
    check_val("dEOP_diu(X,   der)", calc_eop(1,1), ref_eop(1,1), 1e-18);
    check_val("dEOP_diu(Y,   der)", calc_eop(2,1), ref_eop(2,1), 1e-18);

    for(int i=0; i<8; ++i) {
        check_val("Arg_oc_tide[" + std::to_string(i) + "]", calc_arg[i], ref_arg[i], 1e-10);
    }
}

int main() {
    std::cout << std::string(95, '=') << "\n";
    std::cout << "VERIFICATION: TERMS_71 (Diurnal/Subdiurnal Tides)\n";
    
    // =========================================================================
    // ТЕСТ 1: Первый дамп (cent = 0.1603770377607929)
    // =========================================================================
    double cent1 = 0.1603770377607929;
    Eigen::MatrixXd ref_eop1(3, 2);
    ref_eop1 << -0.2005352044783782e-04,  0.1291117030440770e-09,
                -0.4757206630941059e-03, -0.7104596731570391e-08,
                 0.2081424338359992e-04, -0.2388812795620709e-07;
    Eigen::VectorXd ref_arg1(8);
    ref_arg1 << 0.6956017380702093, 0.4635503091573969, 2.853621504954258, 0.5554859990429116,
                1.251087737113121,  1.019036308200308,  3.409107504001793, 1.110971998085823;
    
    run_test(cent1, ref_eop1, ref_arg1);

    // =========================================================================
    // ТЕСТ 2: Второй дамп (cent = 0.1603633667275078)
    // =========================================================================
    double cent2 = 0.1603633667275078;
    Eigen::MatrixXd ref_eop2(3, 2);
    ref_eop2 << -0.2622717808331326e-04,  0.1786378767565643e-08,
                -0.3306936145579745e-03,  0.1468487669796942e-07,
                 0.1151092950963298e-03, -0.4205067645305820e-07;
    Eigen::VectorXd ref_arg2(8);
    ref_arg2 << 4.176313186204166, 3.830399938910039, 6.007985623720407, 3.692670223335277,
                1.585798102355233, 1.239884855061106, 3.417470539880720, 1.102155139490967;

    run_test(cent2, ref_eop2, ref_arg2);

    std::cout << std::string(95, '=') << "\n";
    return 0;
}