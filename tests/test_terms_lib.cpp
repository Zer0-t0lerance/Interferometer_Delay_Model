#include <iostream>
#include <iomanip>
#include "..\\src\\functions.h"

int main() {
    // Формат вывода
    std::cout << std::scientific << std::setprecision(16);
    const Eigen::IOFormat fmt(16, 0, " ", "\n", "", "", "", "");

    // Входные данные
    double cent = 0.1603770377607929;
    Eigen::VectorXd f(5), fd(5);
    f << 0.1248135956975579e+07, 0.3942055124060065e+05, 0.6760087179270387e+06, 0.2464937860590816e+06, -0.6665271692963659e+06;
    fd << 0.1717915933447163e+10, 0.1295965808706691e+09, 0.1739527258757720e+10, 0.1602961599166112e+10, -0.6962888145768072e+07;

    Eigen::MatrixXd dEOP_lib(3, 2);

    // Вызов функции
    ariadna::terms_lib(cent, f, fd, dEOP_lib);

    // Ожидаемые значения из Fortran
    Eigen::MatrixXd expected_dEOP_lib(3, 2);
    expected_dEOP_lib << 0.1348781963026778e-06, 0.6400216191562147e-10,
                         -0.3514768458611694e-05, 0.6224462386135502e-09,
                          0.7596889981598940e-05, 0.2548981056098484e-09;

    // Вывод результатов для dEOP_lib
    std::cout << "dEOP_lib (rows: UT1, X, Y; cols: value, derivative):\n";
    std::cout << "Expected:\n" << expected_dEOP_lib.format(fmt) << "\n";
    std::cout << "Computed:\n" << dEOP_lib.format(fmt) << "\n\n";

    // Проверка результатов
    double tol = 1e-10;
    bool pass = true;
    if ((dEOP_lib - expected_dEOP_lib).norm() > tol) {
        std::cout << "Test failed: dEOP_lib mismatch\n";
        pass = false;
    }

    if (pass) {
        std::cout << "All tests passed!\n";
    }

    return pass ? 0 : 1;
}