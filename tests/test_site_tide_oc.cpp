#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include "../src/functions.h"

using namespace ariadna;

void compare_vectors(const std::string& name, const Eigen::Vector3d& calc, const Eigen::Vector3d& ref) {
    Eigen::Vector3d diff = calc - ref;
    double err = diff.norm();
    printf("%-12s | Err: %.3e \n", name.c_str(), err);
    printf("  Calc: % 21.15e % 21.15e % 21.15e\n", calc.x(), calc.y(), calc.z());
    printf("  Ref : % 21.15e % 21.15e % 21.15e\n", ref.x(), ref.y(), ref.z());
}

int main() {
    printf("VERIFICATION: Ocean Tide Loading (Isolated Math Test)\n");
    printf("---------------------------------------------------------------------------\n");

    // Точные данные из твоего main_tst.cpp
    double jd = 2457401.5;
    double ut1_sec = 66630.0551194871;

    // Вшитые данные для станции j1=0 из твоего старого SITE_TIDE_OC.cpp
    const double raw_amp[11][3] = {
        {2.68E-3, 1.62E-3, 1.17E-3}, {8.20E-4, 7.70E-4, 6.30E-4}, {9.00E-4, 3.70E-4, 2.20E-4},
        {3.40E-4, 2.00E-4, 1.70E-4}, {4.15E-3, 2.79E-3, 6.60E-4}, {3.57E-3, 1.65E-3, 4.60E-4},
        {1.44E-3, 8.70E-4, 2.10E-4}, {9.40E-4, 3.50E-4, 1.50E-4}, {6.50E-4, 5.00E-5, 1.50E-4},
        {3.50E-4, 1.00E-5, 9.01E-5}, {3.00E-4, 1.00E-5, 6.999999999999999E-5}
    };
    const double raw_phs[11][3] = {
        {110.70, 80.60, 36.80}, {6.0, 113.6, 94.30}, {87.20, 67.70, 7.8},
        {9.4, 112.0, 93.90}, {-7.3, -169.8, -12.60}, {-8.8, 172.2, -20.20},
        {-6.6, -170.9, -12.90}, {-16.50, 158.2, -29.10}, {-171.7, -48.50, -168.6},
        {-174.8, -92.40, -173.8}, {-178.5, -158.0, -178.4}
    };

    OceanTideData tide_data;
    for (int wave = 0; wave < 11; ++wave) {
        for (int axis = 0; axis < 3; ++axis) {
            tide_data.amplitudes(axis, wave) = raw_amp[wave][axis];
            tide_data.phases(axis, wave)     = raw_phs[wave][axis];
        }
    }

    // Матрицы задаем и СРАЗУ ТРАНСПОНИРУЕМ (имитируем Fortran Column-Major)
    Eigen::Matrix3d vw_i;
    vw_i << 
        -0.650092018774007,  0.718150315518193, -0.248275031864777,
        -0.741362690371940, -0.671104583000652,  0.0,
        -0.166618511729085,  0.184061845575450,  0.968689583175407;
    vw_i.transposeInPlace();

    Eigen::Matrix3d r2000;
    r2000 << 
         0.856407603342079,       0.516298660270934,      -1.307799799452249E-003,
        -0.516297971645348,       0.856408600544260,       8.446234191276613E-004,
         1.556088935770033E-003, -4.812753432626145E-005,  0.999998788134748;
    r2000.transposeInPlace();

    Eigen::Matrix3d dr2000_dt;
    dr2000_dt << 
        -3.764904133986157E-005,  6.245029944245220E-005,  6.168655101856225E-008,
        -6.245022671509238E-005, -3.764909170451489E-005,  9.537579702956507E-008,
        -7.698574166572183E-011, -5.754292251350311E-011,  1.170274036610182E-013;
    dr2000_dt.transposeInPlace();

    Eigen::Vector3d dx_calc, dv_calc;
    SITE_TIDE_OC(jd, ut1_sec, tide_data, vw_i, r2000, dr2000_dt, dx_calc, dv_calc);

    // Эталонные ответы из твоего Фортран дампа
    Eigen::Vector3d dx_ref(7.170258233186649E-004, 3.068970399038139E-003, -8.305815221019871E-004);
    Eigen::Vector3d dv_ref(-6.807624607450670E-007, -5.604918627419783E-008, 6.721200266427836E-008);

    compare_vectors("dx_octide", dx_calc, dx_ref);
    printf("\n");
    compare_vectors("dv_octide", dv_calc, dv_ref);

    return 0;
}