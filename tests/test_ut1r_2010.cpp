#include <iostream>
#include <iomanip>
#include "..\\src\\functions.h"

int main() {
    // Set output format
    std::cout << std::scientific << std::setprecision(16);
    const Eigen::IOFormat fmt(16, 0, " ", "\t", "", "", "", "");

    // Initialize input data
    Eigen::VectorXd f(5);
    Eigen::VectorXd fd(5);
//    f << 0.1141344609264314E07, 0.3136439828193560E05,  0.5678739386129379E06 , 0.1468483807328939E06 , -0.6660943330937190E06;
//    fd << 0.1717915933443197E10 , 0.1295965808707379E09 , 0.1739527258759306E10 , 0.1602961599166904E10 , -0.6962888146697526E07;

    double JD = 0.2457400500000000E07;
    double ct = 0.7891666666666666E-03;
    double cent;
    
    double dut, dlod, domega;

    ariadna::fund_arg(JD, ct, cent, f, fd);

//    std::cout << "f: " << f.format(fmt) << "\t" << std::endl;

//    f << 0.1141344609264314E07, 0.3136439828193560E05,  0.5678739386129379E06 , 0.1468483807328939E06 , -0.6660943330937190E06;
//    fd << 0.1717915933443197E10 , 0.1295965808707379E09 , 0.1739527258759306E10 , 0.1602961599166904E10 , -0.6962888146697526E07;

//    std::cout << "f: " << f.format(fmt) << "\t" << std::endl;

    ariadna::ut1r_2010(f, dut, dlod, domega);

    // Expected values
    double expected_dut = -0.1003213994142547e-1;
    double expected_domega = -0.2436398437388438e-12;

    // Print results
    std::cout << "e_dut:\t\t" << expected_dut << "\n";
    std::cout << "dut:\t\t" << dut << "\n";
//    std::cout << "e_dlod:\t\t" << expected_dlod << "\n";
    std::cout << "dlod:\t\t" << dlod << "\n";
    std::cout << "e_domega:\t" << expected_domega << "\n";
    std::cout << "domega:\t\t" << domega << "\n";

    // Verify results
    double tol = 1e-10;
    bool pass = true;
    if (std::abs(dut - expected_dut) > tol) {
        std::cout << "Test failed: dut mismatch\n";
        pass = false;
    }
    if (std::abs(domega - expected_domega) > tol) {
        std::cout << "Test failed: domega mismatch\n";
        pass = false;
    }

    if (pass) {
        std::cout << "All tests passed!\n";
    }

    return pass ? 0 : 1;
}