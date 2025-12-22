#include <iostream>
#include "../src/functions.h"

int main() {
    // Initialize observation
    ariadna::Observation obs;
    obs.sta1 = 0;
    obs.sta2 = 1;
    obs.t1 = 30.427;
    obs.t2 = 25.884;
    obs.p1 = 755.087095978288;
    obs.p2 = 645.977873180360;
    obs.e1 = 60.973;
    obs.e2 = 43.395;

    // Initialize stations
    std::vector<ariadna::Station> stations(2);
    stations[0].axsty = "AZEL";
    stations[0].offs = 0.0;
    stations[0].lat_geod = -0.0676813656826229;
    stations[0].h_geod = 23.7768978550448;
    stations[1].axsty = "AZEL";
    stations[1].offs = 0.0;
    stations[1].lat_geod = -0.0451861105617061;
    stations[1].h_geod = 1410.12554619368;

    // Initialize r2000 as 9x3 matrix
    Eigen::MatrixXd r2000(9, 3);
    r2000 << 
        0.9988360617015190, -0.04820310447537553, -0.001727009998570988,
        0.04820309780423414, 0.9988375540050903, -0.00004551047279605037,
        0.001727196188789892, -0.00003778973045624359, 0.9999985076815172,
        0.0, 0.0, 0.0,  // First derivative (zero)
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,  // Second derivative (zero)
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0;

    // Initialize k_star
    std::vector<Eigen::Vector3d> k_star(1);
    k_star[0] << 0.934071715711817, 0.08020509321046379, 0.347961453224755;

    // Initialize vw
    std::vector<Eigen::Matrix3d> vw(2);
    vw[0] << 0.7816194029259159, -0.6200784885114869, -0.06762970541887879,
             0.6215014218927669, 0.7834130344749626, 0.0,
             0.05298199274285165, -0.04203195808002213, 0.9977104905457072;
    vw[1] << 0.7966467840659780, 0.4179699531674714, -0.4366406070059412,
             -0.4645989498200183, 0.8855212113925539, 0.0,
             0.3866545192590811, 0.2028627674637356, 0.8996360265760166;

    // Initialize e and az
    Eigen::MatrixXd e(2, 2);
    e << 0.674516446350307, 0.693637398629522,
         -0.00006524628647339084, -0.00006196985962630866;

    e.transposeInPlace();

    Eigen::MatrixXd az(2, 2);
    az << 1.04390201532089, 5.85310838945029,
          -0.00001499847631644231, -0.00002811691949922670;

    az.transposeInPlace();

    // Output matrices
    Eigen::MatrixXd doff_dl(2, 2);
    Eigen::MatrixXd d_dax(2, 2);
    Eigen::MatrixXd dtau_off(2, 2);

    // Call mount_tel
    ariadna::mount_tel(obs, r2000, stations, k_star, vw, e, az, doff_dl, d_dax, dtau_off);

    // Print results
    std::cout << "doff_dl:\n" << doff_dl << "\n\n";
    std::cout << "d_dax:\n" << d_dax << "\n\n";
    std::cout << "dtau_off:\n" << dtau_off << std::endl;

    return 0;
}