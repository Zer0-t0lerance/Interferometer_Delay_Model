#include <iostream>
#include <iomanip>
#include "..\\src\\functions.h"

int main() {
    // Формат вывода
    std::cout << std::scientific << std::setprecision(16);
    const Eigen::IOFormat fmt(16, 0, " ", "\n", "", "", "", "");

    // Входные данные
    double cent = 0.1800604837452786;
    Eigen::VectorXd f(5), fd(5);
    f << 70641.4512314200, 1294327.85036814, 1219899.54928458, 694301.835983098, -803580.799113102 ; // F1-F5
    fd << 1717915934.70318, 129596580.848894, 1739527258.25572, 1602961598.91545, -6962887.85145636; // F1-F5 derivatives

    Eigen::MatrixXd dEOP_diu(3, 2);
    Eigen::VectorXd arg_oc_tide(8);

    // Вызов функции
    ariadna::terms_71(cent, f, fd, dEOP_diu, arg_oc_tide);

    // Ожидаемые значения из Fortran
    Eigen::MatrixXd expected_dEOP_diu(3, 2);
    expected_dEOP_diu << 0.2448354017392617e-04, 0.1652554432424468e-08,
                         -0.6833177121466402e-03, 0.2319958527221446e-07,
                          -0.3984824901420720e-03, -0.2660757768984708e-07;

    Eigen::VectorXd expected_arg_oc_tide(8);
    expected_arg_oc_tide << 0.1870410869779353e+01, 0.2212890289881575e+01, 0.2661845560785219e+01,
                            0.2212890289881575e+01, 0.1836856368194230e+01, 0.2179335788296453e+01,
                            0.2628291059200096e+01, 0.6216076304009341e+01;

    // Вывод результатов для dEOP_diu
    std::cout << "dEOP_diu (rows: UT1, X, Y; cols: value, derivative):\n";
    std::cout << "Expected:\n" << expected_dEOP_diu.format(fmt) << "\n";
    std::cout << "Computed:\n" << dEOP_diu.format(fmt) << "\n\n";

    // Вывод результатов для arg_oc_tide
    std::cout << "arg_oc_tide (8 specific tidal arguments):\n";
    std::cout << "Expected:\n" << expected_arg_oc_tide.format(fmt) << "\n";
    std::cout << "Computed:\n" << arg_oc_tide.format(fmt) << "\n\n";

    // Проверка результатов
    double tol = 1e-10;
    bool pass = true;
    if ((dEOP_diu - expected_dEOP_diu).norm() > tol) {
        std::cout << "Test failed: dEOP_diu mismatch\n";
        pass = false;
    }
    if ((arg_oc_tide - expected_arg_oc_tide).norm() > tol) {
        std::cout << "Test failed: arg_oc_tide mismatch\n";
        pass = false;
    }

    if (pass) {
        std::cout << "All tests passed!\n";
    }

    return pass ? 0 : 1;
}