#include <iostream>
#include <iomanip>
#include <cmath>
#include "../src/functions.h"

void check_val(const std::string& name, double calc, double ref, double tol = 1e-12) {
    double diff = std::abs(calc - ref);
    std::cout << std::left << std::setw(15) << name 
              << " | Calc: " << std::scientific << std::setprecision(14) << std::setw(22) << calc 
              << " | Ref: "  << std::setw(22) << ref 
              << " | Diff: " << std::setw(12) << diff 
              << (diff <= tol ? " [OK]" : " [FAIL]") << std::endl;
}

int main() {
    ariadna::Observation obs;
    obs.sta1 = 0;
    obs.sta2 = 1;
    
    // Метеоданные в миллибарах (как в дампе)
    obs.t1 = 30.427;
    obs.t2 = 25.884;
    obs.p1 = 1006.700;
    obs.p2 = 861.233;
    obs.e1 = 60.973;
    obs.e2 = 43.395;

    std::vector<ariadna::Station> stations(2);
    stations[0].axsty = "AZEL";
    stations[0].offs = 0.00637;
    stations[0].lat_geod = -0.0676813656826229; 
    
    stations[1].axsty = "AZEL";
    stations[1].offs = 1.49500;
    stations[1].lat_geod = -0.451861105617061; 

    // r2000 (без производных, так как их нет в дампе)
    Eigen::MatrixXd r2000(9, 3);
    r2000.setZero(); 
    r2000.block<3, 3>(0, 0) << 0.9988360617015190, -0.04820310447537553, -0.001727009998570988,
                               0.04820309780423414, 0.9988375540050903, -0.00004551047279605037,
                               0.001727196188789892, -0.00003778973045624359, 0.9999985076815172;

    std::vector<Eigen::Vector3d> k_star(1);
    k_star[0] << 0.934071715711817, 0.08020509321046379, 0.347961453224755;

    // ВАЖНО: Матрицы vw в Фортране выводились ТРАНСПОНИРОВАННЫМИ.
    std::vector<Eigen::Matrix3d> vw(2);
    Eigen::Matrix3d vw0_t;
    vw0_t << 0.7816194029259159, -0.6200784885114869, -0.06762970541887879,
             0.6215014218927669, 0.7834130344749626, 0.0,
             0.05298199274285165, -0.04203195808002213, 0.9977104905457072;
    vw[0] = vw0_t.transpose();

    Eigen::Matrix3d vw1_t;
    vw1_t << 0.7966467840659780, 0.4179699531674714, -0.4366406070059412,
             -0.4645989498200183, 0.8855212113925539, 0.0,
             0.3866545192590811, 0.2028627674637356, 0.8996360265760166;
    vw[1] = vw1_t.transpose();

    Eigen::MatrixXd e(2, 2);
    e << 0.674516446350307, 0.693637398629522, 0.0, 0.0;
    e.transposeInPlace();

    Eigen::MatrixXd az(2, 2);
    az << 1.04390201532089, 5.85310838945029, 0.0, 0.0;
    az.transposeInPlace();

    Eigen::MatrixXd doff_dl(2, 2), d_dax(2, 2), dtau_off(2, 2);
    ariadna::mount_tel(obs, r2000, stations, k_star, vw, e, az, doff_dl, d_dax, dtau_off);

    // --- ЭТАЛОНЫ ИЗ ТВОЕГО ПОСЛЕДНЕГО ДАМПА ---
    double ref_doff_00 = 0.7807257583313549;
    double ref_doff_10 = 0.7687082670063325;
    double ref_ddax_00 = 0.2604220811756895e-8;
    double ref_ddax_10 = -0.2564134775553068e-8;
    double ref_dtau_00 = 0.1658888657089142e-10;
    double ref_dtau_10 = -0.3833381489451837e-8;

    std::cout << std::string(90, '=') << "\n";
    std::cout << "VERIFICATION: MOUNT_TEL (Instrumental Delay & Refraction)\n";
    std::cout << std::string(90, '-') << "\n";

    check_val("doff_dl(0, 0)", doff_dl(0, 0), ref_doff_00);
    check_val("doff_dl(1, 0)", doff_dl(1, 0), ref_doff_10);
    std::cout << std::string(90, '-') << "\n";
    check_val("d_dax(0, 0)", d_dax(0, 0), ref_ddax_00);
    check_val("d_dax(1, 0)", d_dax(1, 0), ref_ddax_10);
    std::cout << std::string(90, '-') << "\n";
    check_val("dtau_off(0,0)", dtau_off(0, 0), ref_dtau_00);
    check_val("dtau_off(1,0)", dtau_off(1, 0), ref_dtau_10);
    std::cout << std::string(90, '=') << "\n";

    return 0;
}