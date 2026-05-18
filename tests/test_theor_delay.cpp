#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cstdio>
#include "../src/functions.h"

using namespace ariadna;

int main() {
    printf("VERIFICATION: theor_delay (STRICT FORTRAN MATCH)\n");
    printf("---------------------------------------------------------------------------\n");

    Eigen::Vector3d K_s(0.9340717157118170, 0.08020509321046379, 0.3479614532247550);
    
    Eigen::Matrix<double, 3, 2> base_line;
    base_line.col(0) << 415219.8984698132, 6610745.296371628, -2340691.244571441;
    base_line.col(1) << -482.0557489968165, 30.57309281858454, 0.8338814000745651;

    std::vector<Eigen::Vector3d> xsta(2), vsta(2), asta(2);
    xsta[0] << 4788183.14731, -4190717.19517, -436901.22171;
    xsta[1] << 5203403.04578, 2420028.10120, -2777592.46628;

    Eigen::Matrix3d Earth, Sun, Moon;
    Earth.col(0) << -30331171034.03649, 132858910822.4568, 57575534750.98912;
    Earth.col(1) << -29619.14420142996, -5774.970368033236, -2504.268390466165;
    Earth.col(2) << -0.02643228993058815, -0.01833876281565307, -0.002325872820870766;

    Sun.col(0) << 268184443.84710, 849738216.96790, 348845058.27214;
    Sun.col(1) << -10.14974, 7.86276, 3.67712;
    Sun.col(2).setZero();

    Moon.col(0) << -30457209591.74120, 133170659155.41185, 57696076560.55701;
    Moon.col(1) << -30654.27062, -6147.62886, -2567.30178;
    Moon.col(2).setZero();

    vsta[0] = Eigen::Vector3d(-29313.55154, -5425.75602, -2504.78267) - Earth.col(1);
    vsta[1] = Eigen::Vector3d(-29795.60728676801, -5395.182928854237, -2503.948787164999) - Earth.col(1);
    asta[1] << 0, 0, 0;
    asta[0] = asta[1] - Eigen::Vector3d(-0.002229426718217218, -0.03515211179078134, 0.000002485990340675793);

    Eigen::Matrix2d Datmc_d, Datmc_w, dtau_off, dt_temp;
    Datmc_d << -0.1225039869341168e-07, 0.1023313654522243e-07, 0.9572258037255259e-12, 0.3346159616445429e-12;
    Datmc_w << -0.1348080932573043e-08, 0.7304528775933354e-09, 0.1057065200668694e-12, 0.2396256388890582e-13;
    dtau_off << 0.1658888657089142e-10, -0.3833381489451837e-08, -0.8349692241350816e-15, -0.8724128831589889e-13;
    dt_temp << 0.8286346437956392e-10, 0.1445568241005169e-09, 0, 0;

    double t2_t1, dt2_t1;
    theor_delay(base_line, xsta, vsta, asta, K_s, Earth, Sun, Moon, Datmc_d, Datmc_w, dtau_off, dt_temp, t2_t1, dt2_t1);

    double ref_t2_t1 = -0.3450839807711632e-03;
    double ref_dt2_t1 = 0.1492797110326536e-05;

    double diff_t = std::abs(t2_t1 - ref_t2_t1);
    double diff_dt = std::abs(dt2_t1 - ref_dt2_t1);

    printf("t2_t1  | Calc: % .15le | Ref: % .15le | Diff: % .15le %s\n", 
           t2_t1, ref_t2_t1, diff_t, (diff_t < 1e-13 ? "[OK]" : "[FAIL]"));
           
    printf("dt2_t1 | Calc: % .15le | Ref: % .15le | Diff: % .15le %s\n", 
           dt2_t1, ref_dt2_t1, diff_dt, (diff_dt < 1e-12 ? "[OK]" : "[FAIL]")); // Допуск для производной чуть мягче

    return 0;
}