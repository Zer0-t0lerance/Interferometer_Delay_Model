#include <iostream>
#include <iomanip>
#include <cmath>
#include "../src/functions.h"

// Вспомогательная функция для красивого вывода и сравнения
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
    obs.t1 = 30.427;
    obs.t2 = 25.884;
    obs.p1 = 755.087095978288;
    obs.p2 = 645.977873180360;
    obs.e1 = 60.973;
    obs.e2 = 43.395;

    std::vector<ariadna::Station> stations(2);
    stations[0].axsty = "AZEL";
    stations[0].offs = 0.0;
    stations[0].lat_geod = -0.0676813656826229;
    
    stations[1].axsty = "AZEL";
    stations[1].offs = 0.0;
    stations[1].lat_geod = -0.0451861105617061;

    Eigen::MatrixXd r2000(9, 3);
    r2000.setZero(); 
    r2000.block<3, 3>(0, 0) << 0.9988360617015190, -0.04820310447537553, -0.001727009998570988,
                               0.04820309780423414, 0.9988375540050903, -0.00004551047279605037,
                               0.001727196188789892, -0.00003778973045624359, 0.9999985076815172;

    std::vector<Eigen::Vector3d> k_star(1);
    k_star[0] << 0.934071715711817, 0.08020509321046379, 0.347961453224755;

    std::vector<Eigen::Matrix3d> vw(2);
    vw[0] << 0.7816194029259159, -0.6200784885114869, -0.06762970541887879,
             0.6215014218927669, 0.7834130344749626, 0.0,
             0.05298199274285165, -0.04203195808002213, 0.9977104905457072;
    vw[1] << 0.7966467840659780, 0.4179699531674714, -0.4366406070059412,
             -0.4645989498200183, 0.8855212113925539, 0.0,
             0.3866545192590811, 0.2028627674637356, 0.8996360265760166;

    Eigen::MatrixXd e(2, 2);
    e << 0.674516446350307, 0.693637398629522,
         -0.00006524628647339084, -0.00006196985962630866;
    e.transposeInPlace();

    Eigen::MatrixXd az(2, 2);
    az << 1.04390201532089, 5.85310838945029,
          -0.00001499847631644231, -0.00002811691949922670;
    az.transposeInPlace();

    Eigen::MatrixXd doff_dl(2, 2);
    Eigen::MatrixXd d_dax(2, 2);
    Eigen::MatrixXd dtau_off(2, 2);

    ariadna::mount_tel(obs, r2000, stations, k_star, vw, e, az, doff_dl, d_dax, dtau_off);

    // --- ЭТАЛОНЫ ИЗ ФОРТРАНА ---
    // ЗАМЕНИ ЭТИ ЧИСЛА НА ТЕ, ЧТО ВЫДАЕТ ФОРТРАН В СВОЕМ ЛОГЕ!
    double ref_doff_00 = 7.80775655027398e-01; 
    double ref_doff_01 = 0.00000000000000e+00;
    double ref_doff_10 = 7.68750981776559e-01;
    double ref_doff_11 = 0.00000000000000e+00;

    double ref_ddax_00 =  2.60438724921959e-09;
    double ref_ddax_01 =  0.00000000000000e+00;
    double ref_ddax_10 = -2.56427725668989e-09;
    double ref_ddax_11 =  0.00000000000000e+00;

    std::cout << std::string(90, '=') << "\n";
    std::cout << "VERIFICATION: MOUNT_TEL (Instrumental Delay & Refraction)\n";
    std::cout << std::string(90, '-') << "\n";

    check_val("doff_dl(0, 0)", doff_dl(0, 0), ref_doff_00);
    check_val("doff_dl(0, 1)", doff_dl(0, 1), ref_doff_01);
    check_val("doff_dl(1, 0)", doff_dl(1, 0), ref_doff_10);
    check_val("doff_dl(1, 1)", doff_dl(1, 1), ref_doff_11);

    std::cout << std::string(90, '-') << "\n";

    check_val("d_dax(0, 0)", d_dax(0, 0), ref_ddax_00);
    check_val("d_dax(0, 1)", d_dax(0, 1), ref_ddax_01);
    check_val("d_dax(1, 0)", d_dax(1, 0), ref_ddax_10);
    check_val("d_dax(1, 1)", d_dax(1, 1), ref_ddax_11);

    std::cout << std::string(90, '=') << "\n";

    return 0;
}