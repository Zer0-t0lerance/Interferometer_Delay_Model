#include <iostream>
#include <iomanip>
#include "..\\src\\functions.h"

int main() {
    // Set output format
    std::cout << std::scientific << std::setprecision(25);
    const Eigen::IOFormat fmt_row(18, 0, " ", " ", "[", "]", "[", "]");

    // Initialize input data
    double jd = 0.2457400500000000E07; // JD from Fortran debug output
//    double jd = 2454465.5; // JD from Fortran debug output
    double ct = 0.7891666666666666E-03; // CT from Fortran debug output
//    double ct = 0.00000000000000000000L; // CT from Fortran debug output
    double cent;
    Eigen::VectorXd f(5), fd(5);

    // Call fund_arg
    ariadna::fund_arg(jd, ct, cent, f, fd);
    // f *= cnst::CARCRAD; // Convert arcseconds to radians

    // Expected values from Fortran debug output
    double expected_cent = 0.07995893223819302;
    
    Eigen::VectorXd expected_f(5);
    expected_f << 0.1141344609264314e7, 0.3136439828193560e5, 0.5678739386129379e6, 0.1468483807328939e6, -0.6660943330937190e6;
    Eigen::VectorXd expected_fd(5);
    expected_fd << 0.1717915933443197e10, 0.1295965808707379e9, 0.1739527258759306e10, 0.1602961599166904e10, -0.6962888146697526e7;

    Eigen::VectorXd expected_f2(5);
    expected_f2 << 2.291187512612069099, 6.212931111003726414, 3.658025792050572989, 4.554139562402433228, -0.5167379217231804489;

    // Print results
//    std::cout << "cent: \t" << cent << "\n";
//    std::cout << "Expected cent: " << expected_cent << "\n";
//    std::cout << "cent diff: " << std::abs(cent - expected_cent) << "\n\n";

    std::cout << "f:  \t\t" << f.transpose().format(fmt_row) << "\n";
    std::cout << "e_f:\t\t" << expected_f.transpose().format(fmt_row) << "\n\n";

    std::cout << "e_fd:\t\t" << expected_fd.transpose().format(fmt_row) << "\n";
    std::cout << "fd:  \t\t" << fd.transpose().format(fmt_row) << "\n\n";

    std::cout << (expected_f - f).transpose().format(fmt_row) << "\n";
//    std::cout << (expected_f2 - f).transpose().format(fmt_row) << "\n";
    std::cout << (expected_fd - fd).transpose().format(fmt_row) << "\n\n";

    // Verify results
    double tol = 1e-17; // Relaxed tolerance
    bool pass = true;
    if (std::abs(cent - expected_cent) > tol) {
        std::cout << "Test failed: cent mismatch\n";
        pass = false;
    }
    if ((f - expected_f).norm() > tol) {
        std::cout << "delta: " << (f - expected_f2).norm() << "\n";
        std::cout << "Test failed: f mismatch\n";
        pass = false;
    }
//    if ((fd - expected_fd).norm() > tol) {
//        std::cout << "Test failed: fd mismatch\n";
//        pass = false;
//    }

    if (pass) {
        std::cout << "All tests passed!\n";
    }

    return pass ? 0 : 1;
}