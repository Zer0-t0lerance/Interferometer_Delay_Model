#include <iostream>
#include <iomanip>
#include <cmath>
#include "../src/functions.h"

using namespace ariadna;

void check_val(const std::string& name, double calc, double ref, double tol = 1e-11) {
    double diff = std::abs(calc - ref);
    std::cout << std::left << std::setw(20) << name 
              << " | Calc: " << std::scientific << std::setprecision(14) << std::setw(22) << calc 
              << " | Ref: "  << std::setw(22) << ref 
              << " | Diff: " << std::setw(12) << diff 
              << (diff <= tol ? " [OK]" : " [FAIL]") << std::endl;
}

int main() {
    Observation obs;
    obs.t1 = 30.427; obs.p1 = 1006.7; obs.e1 = 60.973;
    obs.t2 = 25.884; obs.p2 = 861.233; obs.e2 = 43.395;

    double jd = 2458121.209643333;
    double ct = 0.0;

    Station sta1, sta2;
    sta1.name = "FORTLEZA"; sta1.lat_geod = -0.06768136568262291; sta1.h_geod = 23.77689785504481;
    sta2.name = "HART15M";  sta2.lat_geod = -0.4518611056170612;  sta2.h_geod = 1410.125546193682;

    Eigen::MatrixXd e(2, 2);
    // Значения из дампа + реверс-инжиниринг производных
    e << 0.6745164463503075, -6.524628647339084e-5,
         0.6936373986295217, -6.196985962630866e-5;

    Eigen::MatrixXd az(2, 2);
    // Азимуты (из реверс-инжиниринга DGrad_N и DGrad_E)
    az << 1.04390201532089, -1.499847631644231e-5,
          5.85310838945029, -2.811691949922670e-5;

    Eigen::MatrixXd datmc_d, datmc_w, datmp_hmf, datmp_wmf, dgrad_n, dgrad_e, zen_dry, zen_wet;

    trop_delay(obs, jd, ct, sta1, sta2, e, az, datmc_d, datmc_w, datmp_hmf, datmp_wmf, dgrad_n, dgrad_e, zen_dry, zen_wet);

    std::cout << std::string(100, '=') << "\n";
    std::cout << "VERIFICATION: TROP_DELAY (Matrix Assembly)\n";
    std::cout << std::string(100, '-') << "\n";

    check_val("Datmp_hmf(0,0)", datmp_hmf(0,0), -1.598074084789667);
    check_val("Datmp_hmf(0,1)", datmp_hmf(0,1),  0.1248708542888822e-03);
    check_val("Datmp_hmf(1,0)", datmp_hmf(1,0),  1.561337611708886);
    check_val("Datmp_hmf(1,1)", datmp_hmf(1,1),  0.5105457980404647e-04);

    check_val("Datmp_wmf(0,0)", datmp_wmf(0,0), -1.599786431904585);
    check_val("Datmp_wmf(0,1)", datmp_wmf(0,1),  0.1254434006748068e-03);
    check_val("Datmp_wmf(1,0)", datmp_wmf(1,0),  1.562832539429761);
    check_val("Datmp_wmf(1,1)", datmp_wmf(1,1),  0.5126884392205215e-04);

    check_val("DGrad_N(0,0)", dgrad_n(0,0),  1.000528076216936);
    check_val("DGrad_N(0,1)", dgrad_n(0,1), -0.6137205168079180e-04);
    check_val("DGrad_N(1,0)", dgrad_n(1,0),  1.699607332296508);
    check_val("DGrad_N(1,1)", dgrad_n(1,1),  0.2098137350458779e-04);

    check_val("DGrad_E(0,0)", dgrad_e(0,0),  1.719851167115435);
    check_val("DGrad_E(0,1)", dgrad_e(0,1), -0.2407629516869754e-03);
    check_val("DGrad_E(1,0)", dgrad_e(1,0), -0.7796338866317820);
    check_val("DGrad_E(1,1)", dgrad_e(1,1), -0.1771236331021595e-03);

    check_val("Zen_dry(0,0)", zen_dry(0,0),  0.7665726395296646e-08);
    check_val("Zen_dry(1,0)", zen_dry(1,0),  0.6554083158236512e-08);
    
    check_val("Zen_wet(0,0)", zen_wet(0,0),  0.8426630615738624e-09);
    check_val("Zen_wet(1,0)", zen_wet(1,0),  0.4673903691945521e-09);

    check_val("Datmc_d(0,0)", datmc_d(0,0), -0.1225039869341168e-07);
    check_val("Datmc_d(0,1)", datmc_d(0,1),  0.9572258037255259e-12);
    check_val("Datmc_d(1,0)", datmc_d(1,0),  0.1023313654522243e-07);
    check_val("Datmc_d(1,1)", datmc_d(1,1),  0.3346159616445429e-12);

    check_val("Datmc_w(0,0)", datmc_w(0,0), -0.1348080932573043e-08);
    check_val("Datmc_w(0,1)", datmc_w(0,1),  0.1057065200668694e-12);
    check_val("Datmc_w(1,0)", datmc_w(1,0),  0.7304528775933354e-09);
    check_val("Datmc_w(1,1)", datmc_w(1,1),  0.2396256388890582e-13);

    std::cout << std::string(100, '=') << "\n";
    return 0;
}