#include <iostream>
#include <iomanip>
#include <cmath>
#include "../src/functions.h"

using namespace ariadna;

void check_val(const std::string& name, double calc, double ref, double tol = 1e-9) {
    double diff = std::abs(calc - ref);
    std::cout << std::left << std::setw(20) << name 
              << " | Calc: " << std::scientific << std::setprecision(14) << std::setw(22) << calc 
              << " | Ref: "  << std::setw(22) << ref 
              << " | Diff: " << std::setw(12) << diff 
              << (diff <= tol ? " [OK]" : " [FAIL]") << std::endl;
}

int main() {
    std::cout << std::string(90, '=') << "\n";
    std::cout << "VERIFICATION: UT1R_2010 (Zonal Tides, 62 harmonics)\n";
    std::cout << std::string(90, '-') << "\n";

    // =========================================================================
    // ТЕСТ: Интегральный прогон (через fund_arg)
    // Источник: Лог Фортрана для JD = 2457400.5, CT = 0.7891666666666666E-03
    // =========================================================================
    double jd1 = 2457400.5;
    double ct1 = 0.7891666666666666E-03;
    
    double ref_dut1    = -0.01003213994142547;
    double ref_domega1 = -0.2436398437388438E-12;

    double cent1;
    Eigen::VectorXd f1(5), fd1(5);
    fund_arg(jd1, ct1, cent1, f1, fd1);

    double act_dut1, act_lod1, act_domega1;
    ut1r_2010(f1, act_dut1, act_lod1, act_domega1);

    std::cout << "TEST 1: Full Chain (fund_arg -> ut1r_2010) for JD=2457400.5\n";
    check_val("DUT (UT1-UT1R)", act_dut1, ref_dut1);
    // Для domega допуск жестче, так как числа порядка 10^-13
    check_val("DOMEGA", act_domega1, ref_domega1, 1e-18);


    std::cout << std::string(90, '-') << "\n";

/*    // =========================================================================
    // ТЕСТ 2: Чистый прогон (изолированная математика)
    // Источник: Вектор 'f' из дампа 16JAN14XE_N004.cop (cent = 0.1603770377607929)
    // =========================================================================
    Eigen::VectorXd f2(5);
    f2 << 0.1248135956975579e+07, 0.3942055124060065e+05, 0.6760087179270387e+06, 0.2464937860590816e+06, -0.6665271692963659e+06;

    double act_dut2, act_lod2, act_domega2;
    ut1r_2010(f2, act_dut2, act_lod2, act_domega2);

    // ВНИМАНИЕ: Если в дампе 16JAN14 есть вывод UT1R, замени эти числа на эталоны из лога!
    // Сейчас тут стоят расчетные значения C++ для проверки консистентности.
    double ref_dut2    = 1.5224785681955839e-04; 
    double ref_domega2 = 1.3414902129532296e-14;

    std::cout << "TEST 2: Pure Math Check (Hardcoded 'f' from 16JAN14XE)\n";
    check_val("DUT (UT1-UT1R)", act_dut2, ref_dut2);
    check_val("DOMEGA", act_domega2, ref_domega2, 1e-18);*/

    std::cout << std::string(90, '=') << "\n";
    return 0;
}