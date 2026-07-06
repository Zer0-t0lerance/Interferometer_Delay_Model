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
    printf("VERIFICATION: Pole Tide Loading (POLE_TIDE)\n");
    printf("---------------------------------------------------------------------------\n");

    // Входные данные из твоего Фортран-дампа
    double cent = 0.160363366727612;
    double lat = -0.250899124329325;
    double lon = 2.30649404743805;
    double x = 2.818958470707882E-002;
    double y = 0.276522249175932;
    double dxdt = -1.439721109240823E-008;
    double dydt = 1.198873834157596E-008;

    // Матрицы транспонируются, т.к. в C++ они по строкам, а в Фортране дамп был по столбцам
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
    
    // Вызов обновленной функции POLE_TIDE
    POLE_TIDE(cent, lat, lon, x, y, dxdt, dydt, vw_i, r2000, dr2000_dt, dx_calc, dv_calc);

    // Эталоны из Фортрана
    Eigen::Vector3d dx_ref(-1.606644517260437E-003, 6.408722967998450E-004, -1.520899371451594E-003);
    Eigen::Vector3d dv_ref(-4.675063934877355E-008, -1.170221638159201E-007, 5.688917305462543E-011);

    compare_vectors("dx_poltide", dx_calc, dx_ref);
    printf("\n");
    compare_vectors("dv_poltide", dv_calc, dv_ref);

    return 0;
}