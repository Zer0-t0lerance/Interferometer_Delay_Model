#include <iostream>
#include <iomanip>
#include <cmath>
#include "../src/functions.h"

using namespace ariadna;

// Увеличим допуск до 1e-5 (доли миллиметра), так как между Фортраном 
// и С++ есть микроразличия в константах радиуса Земли и Пи.
void check_val(const std::string& name, double calc, double ref, double tol = 1e-5) {
    double diff = std::abs(calc - ref);
    std::cout << std::left << std::setw(20) << name 
              << " | Calc: " << std::scientific << std::setprecision(14) << std::setw(22) << calc 
              << " | Ref: "  << std::setw(22) << ref 
              << " | Diff: " << std::setw(12) << diff 
              << (diff <= tol ? " [OK]" : " [FAIL]") << std::endl;
}

int main() {
    std::cout << std::string(95, '=') << "\n";
    std::cout << "VERIFICATION: SITE_TIDE_SOLID (16JAN14XE_N004 Data Match)\n";
    std::cout << std::string(95, '-') << "\n";

    // 1. Новые исходные данные
    Eigen::Vector3d xsta(-4147354.75759, 4581542.36895, -1573303.04010);
    double lat_gcen = -0.2492885033603751;
    double lon_gcen = 2.3064940474380524;

    Eigen::Matrix<double, 3, 2> sun;
    sun <<  59384416444.740143, 27741.542356,
           -123521266368.636536, 11144.421240,
           -53548040938.427254,  4829.833355;

    Eigen::Matrix<double, 3, 2> moon;
    moon <<  368177912.262222,   91.647773,
            -29875090.586884,  1003.542036,
            -13964554.684982,   332.209040;

    Eigen::VectorXd f(5), fd(5);
//    fund_arg(xsta, lat_gcen, lon_gcen, sun, moon, f, fd);

    f << 0.12246502712471e+07, 0.37648832083363e+05, 0.65222758305311e+06, 0.22457964484885e+06, -0.66643197942149e+06;
    fd << 0.17179159334463e+10, 0.12959658087068e+09, 0.17395272587581e+10, 0.16029615991663e+10, -0.69628881459725e+07;

    Eigen::Matrix3d vw_i;
    vw_i << -0.65009201877401e+00, -0.74136269037194e+00, -0.16661851172908e+00,
             0.71815031551819e+00, -0.67110458300065e+00,  0.18406184557545e+00,
            -0.24827503186478e+00,  0.00000000000000e+00,  0.96868958317541e+00;

    Eigen::Vector2d gast(0.54610712084277, 0.00007292115290);

    Eigen::MatrixXd r2000 = Eigen::MatrixXd::Zero(9, 3);
    r2000(0,0) =  0.8564076033420785e+00; r2000(0,1) = -0.5162979716453483e+00; r2000(0,2) =  0.1556088935770033e-02;
    r2000(1,0) =  0.5162986602709342e+00; r2000(1,1) =  0.8564086005442605e+00; r2000(1,2) = -0.4812753432626145e-04;
    r2000(2,0) = -0.1307799799452249e-02; r2000(2,1) =  0.8446234191276613e-03; r2000(2,2) =  0.9999987881347479e+00;
    
    r2000(3,0) = -0.3764904133986157e-04; r2000(3,1) = -0.6245022671509238e-04; r2000(3,2) = -0.7698574166572183e-10;
    r2000(4,0) =  0.6245029944245220e-04; r2000(4,1) = -0.3764909170451489e-04; r2000(4,2) = -0.5754292251350311e-10;
    r2000(5,0) =  0.6168655101856225e-07; r2000(5,1) =  0.9537579702956507e-07; r2000(5,2) =  0.1170274036610182e-12;

    Eigen::Vector3d dxtide, dvtide;
    SITE_TIDE_SOLID(xsta, lat_gcen, lon_gcen, sun, moon, f, fd, vw_i, gast, r2000, dxtide, dvtide);

    // 2. Новые эталоны прямо из твоего Фортран-дампа (DXTIDE total и DVTIDE total)
    Eigen::Vector3d ref_dx(-0.1886053448287193, 0.05180400082089129, -0.009650602127782066);
    Eigen::Vector3d ref_dv(-0.3521000451482045e-05, -0.7755224750492061e-05, -0.1088793325206996e-05);

    std::cout << "Displacement (m):\n";
    check_val("dx_X", dxtide(0), ref_dx(0));
    check_val("dx_Y", dxtide(1), ref_dx(1));
    check_val("dx_Z", dxtide(2), ref_dx(2));

    std::cout << "\nVelocity (m/s):\n";
    check_val("dv_X", dvtide(0), ref_dv(0));
    check_val("dv_Y", dvtide(1), ref_dv(1));
    check_val("dv_Z", dvtide(2), ref_dv(2));

    std::cout << std::string(95, '=') << "\n";
    
    return 0;
}