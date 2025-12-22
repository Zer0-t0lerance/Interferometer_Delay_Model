// test_trop_delay.cpp
#include <iostream>
#include <iomanip>
#include "..\\src\\functions.h"

int main() {
    // Set output format
    std::cout << std::scientific << std::setprecision(16);
    const Eigen::IOFormat fmt(16, 0, " ", "\n", "", "", "", "");

    // Initialize input data
    ariadna::Observation obs;
    obs.t1 = 30.427; // tc(1) in Celsius
    obs.t2 = 25.884; // tc(2)
    obs.p1 = 1006.7; // pres(1) in mm
    obs.p2 = 861.233; // pres(2)
    obs.e1 = 60.973; // humid(1) in %
    obs.e2 = 43.395; // humid(2)
    obs.kw1 = 0; // Relative humidity
    obs.kw2 = 0;

    double jd = 2458121.209643333; // epoch
    double ct = 0.0; // No coordinate time offset for test

    ariadna::Station sta1, sta2;
    sta1.lat_geod = -0.06768136568262291; // lat_gd(1) in rad
    sta2.lat_geod = -0.4518611056170612; // lat_gd(2)
    sta1.h_geod = 23.77689785504481; // h_g(1) in meters
    sta2.h_geod = 1410.125546193682; // h_g(2)
    sta1.name = "STATION1"; // Ground station
    sta2.name = "STATION2";

    Eigen::MatrixXd e(2, 2);
    e << 0.6745164463503079, -6.524628647339084E-005, // E(1,1), approximate E(1,2)
         0.6936373986295212, -6.196985962630866E-005; // E(2,1), approximate E(2,2)

    Eigen::MatrixXd az(2, 2);
    az << 1.04390201532089, -1.499847631644231E-005,
         5.85310838945029, -2.811691949922670E-005;

    // Initialize output matrices
    Eigen::MatrixXd datmc_d(2, 2), datmc_w(2, 2);
    Eigen::MatrixXd datmp_hmf(2, 2), datmp_wmf(2, 2);
    Eigen::MatrixXd dgrad_n(2, 2), dgrad_e(2, 2);
    Eigen::MatrixXd zen_dry(2, 2), zen_wet(2, 2);

    // Call trop_delay
    ariadna::trop_delay(obs, jd, ct, sta1, sta2, e, az, datmc_d, datmc_w, datmp_hmf, datmp_wmf, dgrad_n, dgrad_e, zen_dry, zen_wet);

    // Expected values from Fortran output
    Eigen::MatrixXd expected_datmp_hmf(2, 2);
    expected_datmp_hmf << -1.598074084789667e+00, 1.248708542888822e-04,
                           1.561337611708886e+00, 5.105457980404647e-05;

    Eigen::MatrixXd expected_datmp_wmf(2, 2);
    expected_datmp_wmf << -1.599786431904585e+00, 1.254434006748068e-04,
                           1.562832539429761e+00, 5.126884392205215e-05;

    Eigen::MatrixXd expected_zen_dry(2, 2);
    expected_zen_dry << 7.665726395296646e-09, 0.0,
                        6.554083158236512e-09, 0.0;

    Eigen::MatrixXd expected_zen_wet(2, 2);
    expected_zen_wet << 8.426630615738624e-10, 0.0,
                        4.673903691945521e-10, 0.0;

    Eigen::MatrixXd expected_datmc_d(2, 2);
    expected_datmc_d << -1.225039869341168e-08, 9.572258037255259e-13,
                         1.023313654522243e-08, 3.346159616445429e-13;

    Eigen::MatrixXd expected_datmc_w(2, 2);
    expected_datmc_w << -1.348080932573043e-09, 1.057065200668694e-13,
                         7.304528775933354e-10, 2.396256388890582e-14;

    // Print results
    std::cout << "datmp_hmf:\n" << datmp_hmf.format(fmt) << "\n\n";
    std::cout << "datmp_wmf:\n" << datmp_wmf.format(fmt) << "\n\n";
    std::cout << "zen_dry:\n" << zen_dry.format(fmt) << "\n\n";
    std::cout << "zen_wet:\n" << zen_wet.format(fmt) << "\n\n";
    std::cout << "datmc_d:\n" << datmc_d.format(fmt) << "\n\n";
    std::cout << "datmc_w:\n" << datmc_w.format(fmt) << "\n\n";

    // Verify results
    double tol = 1e-10;
    bool pass = true;
    if ((datmp_hmf - expected_datmp_hmf).norm() > tol) {
        std::cout << "Test failed: datmp_hmf mismatch\n";
        pass = false;
    }
    if ((datmp_wmf - expected_datmp_wmf).norm() > tol) {
        std::cout << "Test failed: datmp_wmf mismatch\n";
        pass = false;
    }
    if ((zen_dry - expected_zen_dry).norm() > tol) {
        std::cout << "Test failed: zen_dry mismatch\n";
        pass = false;
    }
    if ((zen_wet - expected_zen_wet).norm() > tol) {
        std::cout << "Test failed: zen_wet mismatch\n";
        pass = false;
    }
    if ((datmc_d - expected_datmc_d).norm() > tol) {
        std::cout << "Test failed: datmc_d mismatch\n";
        pass = false;
    }
    if ((datmc_w - expected_datmc_w).norm() > tol) {
        std::cout << "Test failed: datmc_w mismatch\n";
        pass = false;
    }

    if (pass) {
        std::cout << "All tests passed!\n";
    }

    return pass ? 0 : 1;
}