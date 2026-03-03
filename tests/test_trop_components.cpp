#include <iostream>
#include <iomanip>
#include <cmath>
#include "../src/functions.h"

using namespace ariadna;

void check_val(const std::string& name, double calc, double ref, double tol = 1e-12) {
    double diff = std::abs(calc - ref);
    std::cout << std::left << std::setw(20) << name 
              << " | Calc: " << std::scientific << std::setprecision(14) << std::setw(22) << calc 
              << " | Ref: "  << std::setw(22) << ref 
              << " | Diff: " << std::setw(12) << diff 
              << (diff <= tol ? " [OK]" : " [FAIL]") << std::endl;
}

int main() {
    std::cout << std::string(90, '=') << "\n";
    std::cout << "VERIFICATION: TROPOSPHERE CORE COMPONENTS\n";
    std::cout << std::string(90, '-') << "\n";

    // Общие данные из дампа
    double epoch = 2458121.209643333;

    // =========================================================================
    // СТАНЦИЯ 1 (FORTLEZA)
    // =========================================================================
    double tc1 = 30.427, pres1 = 1006.700, humid1 = 60.973;
    double lat1 = -0.06768136568262291, h1 = 23.77689785504481, e1 = 0.6745164463503075;

    double z_d1, dot_z_d1, dz_ddh1, z_w1, dot_z_w1;
    Eigen::Vector2d hmf1, wmf1;

    // В trop_delay dpdh вычисляется по формуле, повторим её для теста
    double dxdh = -6.5e-3 / 293.15;
    double dpdh1 = 1013.25 * 5.26 * std::pow(1.0 - 6.5e-3 * h1 / 293.15, 4.26) * dxdh;

    sast_dry(pres1, 0.0, lat1, h1, dpdh1, z_d1, dot_z_d1, dz_ddh1);
    sast_wet(humid1, tc1, 0.0, 0.0, z_w1, dot_z_w1);
    nhmf2(epoch, lat1, h1, e1, hmf1);
    nwmf2(lat1, e1, wmf1);

    std::cout << "--- STATION 1 ---\n";
    check_val("Z_d_1 (m)", z_d1, 2.298126958401461);
    check_val("Z_w_1 (m)", z_w1, 0.2526240304950336);
    check_val("HMF_1", hmf1(0), 1.598074084789667);
    check_val("WMF_1", wmf1(0), 1.599786431904585);

    // =========================================================================
    // СТАНЦИЯ 2 (HART15M)
    // =========================================================================
    double tc2 = 25.884, pres2 = 861.233, humid2 = 43.395;
    double lat2 = -0.4518611056170612, h2 = 1410.125546193682, e2 = 0.6936373986295217;

    double z_d2, dot_z_d2, dz_ddh2, z_w2, dot_z_w2;
    Eigen::Vector2d hmf2, wmf2;

    double dpdh2 = 1013.25 * 5.26 * std::pow(1.0 - 6.5e-3 * h2 / 293.15, 4.26) * dxdh;

    sast_dry(pres2, 0.0, lat2, h2, dpdh2, z_d2, dot_z_d2, dz_ddh2);
    sast_wet(humid2, tc2, 0.0, 0.0, z_w2, dot_z_w2);
    nhmf2(epoch, lat2, h2, e2, hmf2);
    nwmf2(lat2, e2, wmf2);

    std::cout << "\n--- STATION 2 ---\n";
    check_val("Z_d_2 (m)", z_d2, 1.964864699944127);
    check_val("Z_w_2 (m)", z_w2, 0.1401201076263623);
    check_val("HMF_2", hmf2(0), 1.561337611708886);
    check_val("WMF_2", wmf2(0), 1.562832539429761);

    std::cout << std::string(90, '=') << "\n";
    return 0;
}