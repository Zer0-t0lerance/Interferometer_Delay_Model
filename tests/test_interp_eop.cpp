#include "..\\src\\functions.h"
#include <iostream>
#include <iomanip>
#include <cmath> // для std::abs

int main() {
    // Настраиваем вывод для научной нотации и высокой точности
    std::cout << std::scientific << std::setprecision(16);
    // const Eigen::IOFormat fmt(16, 0, " ", "\t", "", "", "", ""); // fmt больше не нужен для табличного вывода

    // Входные данные
    int k_int = 0;  // Spline 
    ariadna::Observation obs;
    obs.mjd = 57402;
    obs.utc = 0.2705150462934398;
    // Другие поля obs не используются в interp_eop

    double tt = 0.2713042129601065; // время в формате TAI

    std::vector<ariadna::EOPData> eop_data(7);
    double mjds[7] = {57400, 57401, 57402, 57403, 57404, 57405, 57406};
    double ut1_utcs[7] = {0.0589932, 0.0568475, 0.0546433, 0.0524926, 0.0504060, 0.0483894, 0.0464968};
    double ut1_tais[7] = {-35.9410068, -35.9431525, -35.9453567, -35.9475074, -35.9495940, -35.9516106, -35.9535032};
    double xs[7] = {0.031105, 0.029487, 0.028220, 0.026690, 0.024996, 0.023229, 0.021297};
    double ys[7] = {0.274426, 0.275674, 0.276671, 0.278122, 0.279528, 0.281340, 0.283184};
    double dpsis[7] = {-0.00009, -0.00006, -0.00015, -0.00017, -0.00015, -0.00014, -0.00013};
    double depss[7] = {-0.00014, -0.00011, -0.00013, -0.00012, -0.00009, -0.00006, -0.00003};

    for (int i = 0; i < 7; ++i) {
        eop_data[i].mjd = mjds[i];
        eop_data[i].ut1_utc = ut1_utcs[i];
        eop_data[i].ut1_tai = ut1_tais[i];
        eop_data[i].x = xs[i];
        eop_data[i].y = ys[i];
        eop_data[i].dpsi = dpsis[i];
        eop_data[i].deps = depss[i];
    }

    // Выходные
    double ut1;
    Eigen::VectorXd eop_int(5);
    Eigen::VectorXd deop_int(5);
    Eigen::VectorXd arg_oc_tide(8);
    Eigen::MatrixXd deop_diu(3, 2);
    Eigen::MatrixXd deop_lib(3, 2);

    // Вызов
    ariadna::interp_eop(k_int, obs, tt, ut1, eop_int, deop_int, arg_oc_tide, deop_diu, deop_lib, eop_data);

    // Ожидаемые (из debug: eop_int[0..4] как UT1-UTC, x, y, dpsi, deps)
    double exp_eop_int[5] = {
        0.05403427263165383,
        0.02736165493348263,
        0.2770582096703436,
        -0.0001603168551743361,
        -0.0001263303748918657
    };

    // --------------------------------------------------------------------------
    // НОВАЯ СЕКЦИЯ: Вывод и проверка в табличном формате
    // --------------------------------------------------------------------------
    
    std::cout << "\n------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Idx | EOP Parameter | Calculated Value (eop_int) | Expected Value (exp_eop_int) | Status" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------------" << std::endl;

    const char* names[5] = {"UT1-UTC", "Polar motion X", "Polar motion Y", "DPSI", "DEPS"};
    double tol = 1e-10;
    bool pass = true;
    
    // Используем setw для форматирования в столбцы. 
    // Ширина 25 символов для чисел (с учетом scientific, знака, экспоненты и 16 знаков).
    for (int i = 0; i < 5; ++i) {
        double calculated = eop_int(i);
        double expected = exp_eop_int[i];
        
        bool current_pass = std::abs(calculated - expected) <= tol;
        if (!current_pass) {
            pass = false;
        }

        std::cout << std::setw(3) << i << " | "
                  << std::setw(15) << std::left << names[i] << std::right << " | "
                  << std::setw(26) << calculated << " | "
                  << std::setw(26) << expected << " | "
                  << (current_pass ? "OK" : "FAIL")
                  << std::endl;
    }
    std::cout << "------------------------------------------------------------------------------------------------" << std::endl;

    if (pass) {
        std::cout << "All EOP interpolation tests passed successfully!" << std::endl;
    } else {
        std::cout << "Some EOP interpolation tests failed." << std::endl;
    }

    return pass ? 0 : 1;
}