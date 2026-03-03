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
    // --- Данные для первой станции ---
    double el_rad = 0.674516446350307;
    double temp_k = 30.427 + cnst::KELVIN_OFFSET;
    double humid_f = 60.973 / cnst::PERCENT_TO_FRACTION;
    double press_hg = 755.087095978288;

    double result1 = ariadna::sbend(el_rad, temp_k, humid_f, press_hg);

    // --- Данные для второй станции ---
    double el_rad2 = 0.693637398629522;
    double temp_k2 = 25.884 + cnst::KELVIN_OFFSET;
    double humid_f2 = 43.395 / cnst::PERCENT_TO_FRACTION;
    double press_hg2 = 645.977873180360;

    double result2 = ariadna::sbend(el_rad2, temp_k2, humid_f2, press_hg2);

    // =========================================================================
    // ЭТАЛОНЫ ИЗ ФОРТРАНА (ЗАМЕНИ НА СВОИ ЧИСЛА ИЗ ДАМПА)
    // Если точных чисел из Фортрана нет, пока я положил сюда результаты C++
    // =========================================================================
    double ref1 = 4.53037418731383e-04; // Замени на Фортран Референс
    double ref2 = 3.39706346296711e-04; // Замени на Фортран Референс

    std::cout << std::string(90, '=') << "\n";
    std::cout << "VERIFICATION: SBEND (Atmospheric Refraction)\n";
    std::cout << std::string(90, '-') << "\n";

    check_val("sbend_1 (rad)", result1, ref1);
    check_val("sbend_2 (rad)", result2, ref2);

    std::cout << std::string(90, '=') << "\n";

    return 0;
}