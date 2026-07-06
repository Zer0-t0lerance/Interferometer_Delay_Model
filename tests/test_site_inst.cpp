#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include "../src/functions.h"

using namespace ariadna;

void compare_vectors(const std::string& name, const Eigen::Vector3d& calc, const Eigen::Vector3d& ref) {
    Eigen::Vector3d diff = calc - ref;
    double err = diff.norm();
    printf("%-15s | Err: %.3e \n", name.c_str(), err);
    printf("  Calc: % 21.15e % 21.15e % 21.15e\n", calc.x(), calc.y(), calc.z());
    printf("  Ref : % 21.15e % 21.15e % 21.15e\n", ref.x(), ref.y(), ref.z());
}

int main() {
    printf("VERIFICATION: SITE_INST (Displacements aggregation to J2000.0)\n");
    printf("---------------------------------------------------------------------------\n");

    // Исходные координаты в ITRF
    Eigen::Vector3d xsta1(-4147354.75759410, 4581542.36895419, -1573303.04010378);
    Eigen::Vector3d xsta2(-2388896.27816111, 5043350.02442442, -3078590.70269967);

    // Матрицы транспонируются, так как Фортран дамп был выведен по столбцам
    Eigen::Matrix3d r2000;
    r2000 << 
         0.856407603342079,       0.516298660270934,      -1.307799799452249E-003,
        -0.516297971645348,       0.856408600544260,       8.446234191276613E-004,
         1.556088935770033E-003, -4.812753432626145E-005,  0.999998788134748;
    r2000.transposeInPlace();

    Eigen::Matrix3d dr2000;
    dr2000 << 
        -3.764904133986157E-005,  6.245029944245220E-005,  6.168655101856225E-008,
        -6.245022671509238E-005, -3.764909170451489E-005,  9.537579702956507E-008,
        -7.698574166572183E-011, -5.754292251350311E-011,  1.170274036610182E-013;
    dr2000.transposeInPlace();

    Eigen::Matrix3d d2r2000;
    d2r2000 << 
        -4.553942307973770E-009, -2.745415038296109E-009,  6.954889614806429E-012,
         2.745411365867157E-009, -4.553947611300586E-009, -4.498109878696921E-012,
         4.290405149448649E-015, -5.725553756373111E-015, -6.961050533090952E-018;
    d2r2000.transposeInPlace();

    // ==========================================
    // СТАНЦИЯ 1 (j1)
    // ==========================================
    Eigen::Vector3d dx_tide1(-0.188605344828719, 5.180400082089129E-002, -9.650602127782066E-003);
    Eigen::Vector3d dx_oc1(7.170258233186649E-004, 3.068970399038139E-003, -8.305815221019871E-004);
    Eigen::Vector3d dx_pol1(-1.606644517260437E-003, 6.408722967998450E-004, -1.520899371451594E-003);
    Eigen::Vector3d dx_atm1(2.584159404070865E-003, 5.296297028081103E-003, -1.532715602586989E-003);
    Eigen::Vector3d dx_temp1(0.0, 0.0, 0.0);

    // В твоем дампе скоростей нет, для теста они не понадобятся (ставим 0)
    Eigen::Vector3d dv_zero(0.0, 0.0, 0.0);
    Eigen::Vector3d xsta_j2000t_1, vsta_j2000t_1, asta_j2000_1;

    SITE_INST(xsta1, r2000, dr2000, d2r2000, dx_tide1, dv_zero, dx_oc1, dv_zero, dx_pol1, dv_zero, dx_atm1, dv_zero, dx_temp1, dv_zero, xsta_j2000t_1, vsta_j2000t_1, asta_j2000_1);

    // Эталон для Станции 1
    Eigen::Vector3d ref_x_j2000t_1(-5919715.56662292, 1782474.36352239, -1564007.55930651);

    // ==========================================
    // СТАНЦИЯ 2 (j2)
    // ==========================================
    Eigen::Vector3d dx_tide2(-9.411056584101973E-002, 3.068191223864809E-002, 6.509665736554310E-003);
    Eigen::Vector3d dx_oc2(-8.893443570289910E-004, -2.931073507656510E-003, 1.898421238644091E-003);
    Eigen::Vector3d dx_pol2(-2.561585420639906E-003, 1.834948799071789E-003, -2.541297039563209E-003);
    Eigen::Vector3d dx_atm2(2.041552619712173E-003, -1.382695330888931E-004, -3.408362951457085E-003);
    Eigen::Vector3d dx_temp2(0.0, 0.0, 0.0);

    Eigen::Vector3d xsta_j2000t_2, vsta_j2000t_2, asta_j2000_2;

    SITE_INST(xsta2, r2000, dr2000, d2r2000, dx_tide2, dv_zero, dx_oc2, dv_zero, dx_pol2, dv_zero, dx_atm2, dv_zero, dx_temp2, dv_zero, xsta_j2000t_2, vsta_j2000t_2, asta_j2000_2);

    // Эталон для Станции 2
    Eigen::Vector3d ref_x_j2000t_2(-4654530.98057079, 3085932.58295861, -3071203.03978917);

    printf("\n--- Station 1 ---\n");
    compare_vectors("X_J2000_T", xsta_j2000t_1, ref_x_j2000t_1);

    printf("\n--- Station 2 ---\n");
    compare_vectors("X_J2000_T", xsta_j2000t_2, ref_x_j2000t_2);

    return 0;
}