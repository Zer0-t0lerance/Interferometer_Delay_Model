#include <iostream>
#include <iomanip>
#include "../src/functions.h"

int main() {
    std::cout << std::scientific << std::setprecision(16);

    // Входные данные из твоего лога
    double cent = 0.1800604837452786;
    
    Eigen::VectorXd f(5), fd(5);
    f << 70641.4512314200, 1294327.85036814, 1219899.54928458, 694301.835983098, -803580.799113102;
    fd << 1717915934.70318, 129596580.848894, 1739527258.25572, 1602961598.91545, -6962887.85145636;

    // Эталонные EOP из Фортрана
    Eigen::MatrixXd ref_eop(3, 2);
    ref_eop << 0.2448354017392617e-04, 0.1652554432424468e-08,
              -0.6833177121466402e-03, 0.2319958527221446e-07,
              -0.3984824901420720e-03, -0.2660757768984708e-07;

    // Эталонные аргументы из Фортрана
    Eigen::VectorXd ref_arg(8);
    ref_arg << 1.870410869779353, 2.212890289881575, 2.661845560785219, 6.249630805594464,
               1.836856368194230, 2.179335788296453, 2.628291059200096, 6.216076304009341;

    Eigen::MatrixXd calc_eop(3, 2);
    Eigen::VectorXd calc_arg(8);

    // Запуск функции
    ariadna::terms_71(cent, f, fd, calc_eop, calc_arg);

    std::cout << "--- VERIFICATION ARGUMENTS ---" << std::endl;
    for(int i=0; i<8; ++i) {
        double diff = std::abs(calc_arg[i] - ref_arg[i]);
        std::cout << "Arg[" << i << "] Diff: " << diff << (diff < 1e-12 ? " [OK]" : " [FAIL]") << std::endl;
    }

    std::cout << "\n--- VERIFICATION EOP VALUES ---" << std::endl;
    for(int i=0; i<3; ++i) {
        double diff_val = std::abs(calc_eop(i,0) - ref_eop(i,0));
        double diff_der = std::abs(calc_eop(i,1) - ref_eop(i,1));
        std::cout << "Row " << i << " Val Diff: " << diff_val << (diff_val < 1e-14 ? " [OK]" : " [FAIL]") << std::endl;
        std::cout << "Row " << i << " Der Diff: " << diff_der << (diff_der < 1e-18 ? " [OK]" : " [FAIL]") << std::endl;
    }

    return 0;
}