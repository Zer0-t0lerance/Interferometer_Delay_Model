#include <iostream>
#include <iomanip>
#include <cmath>
#include "../src/functions.h"

// Наша стандартная функция для красивого вывода
void check_val(const std::string& name, double calc, double ref, double tol = 1e-12) {
    double diff = std::abs(calc - ref);
    std::cout << std::left << std::setw(15) << name 
              << " | Calc: " << std::scientific << std::setprecision(14) << std::setw(22) << calc 
              << " | Ref: "  << std::setw(22) << ref 
              << " | Diff: " << std::setw(12) << diff 
              << (diff <= tol ? " [OK]" : " [FAIL]") << std::endl;
}

int main() {
    std::cout << std::string(90, '=') << "\n";
    std::cout << "VERIFICATION: SBEND (Atmospheric Refraction)\n";
    std::cout << std::string(90, '-') << "\n";

    // --- Данные для первой станции (строго из дампа Фортрана) ---
    double el_rad1   = 0.674516446350307;
    double temp_k1   = 303.587000000000;
    double humid_f1  = 0.609730000000000;
    double press_hg1 = 755.087095978288;
    double ref_rho1  = 4.534566209457922e-04;

    double calc_rho1 = ariadna::sbend(el_rad1, temp_k1, humid_f1, press_hg1);
    check_val("sbend_1 (rad)", calc_rho1, ref_rho1);

    // --- Данные для второй станции (строго из дампа Фортрана) ---
    double el_rad2   = 0.693637398629522;
    double temp_k2   = 299.044000000000;
    double humid_f2  = 0.433950000000000;
    double press_hg2 = 645.977873180360;
    double ref_rho2  = 3.398313283113267e-04;

    double calc_rho2 = ariadna::sbend(el_rad2, temp_k2, humid_f2, press_hg2);
    check_val("sbend_2 (rad)", calc_rho2, ref_rho2);

    std::cout << std::string(90, '=') << "\n";

    return 0;
}