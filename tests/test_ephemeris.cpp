#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cstdio>
#include "../src/functions.h"

using namespace ariadna;

void compare_vectors(const std::string& name, const Eigen::Vector3d& calc, const Eigen::Vector3d& ref) {
    Eigen::Vector3d diff = calc - ref;
    double err = diff.norm();
    
    printf("%-12s | Err: %.3e \n", name.c_str(), err);
    printf("  Calc: % 18.3f % 18.3f % 18.3f\n", calc.x(), calc.y(), calc.z());
    printf("  Ref : % 18.3f % 18.3f % 18.3f\n", ref.x(), ref.y(), ref.z());
}

int main() {
    printf("VERIFICATION: JPL Ephemeris Adapter (DE440 vs Fortran DE Dump)\n");
    printf("---------------------------------------------------------------------------\n");

    try {
        init_ephemeris("external/dephem-master/linux_p1550p2650.440t");
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }

    // Время сдвинуто на ~17:03:08 UTC, чтобы попасть ровно в момент наблюдения Фортрана
    double jd = 2458120.5;
    double ct = 0.7113110196885; 

    Eigen::Matrix3d Earth, Sun, Moon;
    get_celestial_bodies(jd, ct, Earth, Sun, Moon);

    Eigen::Vector3d ref_Earth_pos(-30331171034.0, 132858910822.4, 57575534750.9);
    Eigen::Vector3d ref_Earth_vel(-29619.144, -5774.970, -2504.268);

    Eigen::Vector3d ref_Sun_pos(268184443.8, 849738216.9, 348845058.2);
    Eigen::Vector3d ref_Sun_vel(-10.149, 7.862, 3.677);

    Eigen::Vector3d ref_Moon_pos(-30457209591.7, 133170659155.4, 57696076560.5);
    Eigen::Vector3d ref_Moon_vel(-30654.270, -6147.628, -2567.301);

    printf("\nNOTE: We calculated the exact fraction of the day (17:03:08) for this observation.\n");
    printf("The position error should drop from 1.8 billion meters to just a few km.\n\n");

    compare_vectors("Earth Pos", Earth.col(0), ref_Earth_pos);
    compare_vectors("Earth Vel", Earth.col(1), ref_Earth_vel);
    printf("\n");
    compare_vectors("Sun Pos", Sun.col(0), ref_Sun_pos);
    compare_vectors("Sun Vel", Sun.col(1), ref_Sun_vel);
    printf("\n");
    compare_vectors("Moon Pos", Moon.col(0), ref_Moon_pos);
    compare_vectors("Moon Vel", Moon.col(1), ref_Moon_vel);

    return 0;
}