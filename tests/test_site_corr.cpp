#include <iostream>
#include <iomanip>
#include <cmath>
#include "../src/functions.h"

using namespace ariadna;

void check_val(const std::string& name, double calc, double ref, double tol = 1e-10) {
    double diff = std::abs(calc - ref);
    std::cout << std::left << std::setw(25) << name 
              << " | Calc: " << std::scientific << std::setprecision(13) << std::setw(21) << calc 
              << " | Ref: "  << std::setw(21) << ref 
              << " | Diff: " << std::setw(12) << diff 
              << (diff <= tol ? " [OK]" : " [FAIL]") << std::endl;
}

int main() {
    std::cout << std::string(95, '=') << "\n";
    std::cout << "VERIFICATION: SITE_CORR (Plate Tectonics & VEN->ITRF Matrix)\n";
    std::cout << std::string(95, '-') << "\n";

    Station sta;
    sta.name = "KATH12M";
    sta.xyz = Eigen::Vector3d(-4147354.357, 4581542.476, -1573303.689);
    sta.vel = Eigen::Vector3d(-0.0363, -0.0097, 0.0588);
    
    double dYear = 11.0356500494334;

    SiteCorrData out;
    SITE_CORR(sta, dYear, out);

    // Истинные значения с полной математической точностью (без округления)
    Eigen::Vector3d ref_xyz(-4147354.7575940967, 4581542.3689541817, -1573303.0401042798);
    double ref_u = 6179893.361851869;
    double ref_v = -1573303.0401042798;

    std::cout << "Current Epoch Coordinates (ITRF):\n";
    check_val("XYZ_X", out.xyz(0), ref_xyz(0));
    check_val("XYZ_Y", out.xyz(1), ref_xyz(1));
    check_val("XYZ_Z", out.xyz(2), ref_xyz(2));

    std::cout << "\nSpherical Radii:\n";
    check_val("U (Equatorial Rad)", out.u_site, ref_u);
    check_val("V (Polar Rad)", out.v_site, ref_v);

    Eigen::Matrix3d ref_vw_i;
    ref_vw_i << -0.6500920187740, -0.7413626903719, -0.1666185117290,
                 0.7181503155181, -0.6711045830006,  0.1840618455754,
                -0.2482750318647,  0.0000000000000,  0.9686895831754;

    std::cout << "\nTransformation Matrix (VEN -> ITRF):\n";
    check_val("vw_i(0,0) [Up_X]", out.vw_i(0,0), ref_vw_i(0,0), 1e-8);
    check_val("vw_i(1,0) [Up_Y]", out.vw_i(1,0), ref_vw_i(1,0), 1e-8);
    check_val("vw_i(2,0) [Up_Z]", out.vw_i(2,0), ref_vw_i(2,0), 1e-8);
    check_val("vw_i(0,1) [East_X]", out.vw_i(0,1), ref_vw_i(0,1), 1e-8);
    check_val("vw_i(1,1) [East_Y]", out.vw_i(1,1), ref_vw_i(1,1), 1e-8);
    
    std::cout << std::string(95, '=') << "\n";
    return 0;
}