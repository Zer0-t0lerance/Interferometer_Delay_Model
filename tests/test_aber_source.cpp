#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "..\\src\\functions.h" 

using namespace std;

void print_result(string name, double act, double exp, string unit = "") {
    double diff = abs(act - exp);
    if (name.find("Azimuth") != string::npos) {
        double d = abs(act - exp);
        double d2pi = abs(d - 6.283185307179586);
        diff = min(d, d2pi);
    }
    cout << left << setw(30) << name 
         << " | Act: " << fixed << setprecision(10) << act 
         << " | Ref: " << exp 
         << " | Diff: " << scientific << setprecision(2) << diff;
    if (diff < 1e-10) cout << " | [ SUCCESS ]" << unit << endl;
    else cout << " | [ FAIL ]" << unit << endl;
}

int main() {
    ariadna::Observation obs;
    obs.sta1 = 0; obs.sta2 = 1;
    Eigen::Matrix2d e, az;

    // 1. Вектор на звезду (K_s)
    Eigen::Vector3d k_s_log(0.2488843299651737, 0.8832588823692126, -0.3973793364200963);
    
    // 2. Матрица R2000 (J2000 -> Crust). 
    // В тесте передаем ТРАНСПОНИРОВАННУЮ, чтобы внутри функции получить оригинал.
    Eigen::Matrix3d r2000_log;
    r2000_log <<  0.856407603342079,  0.516298660270934, -0.001307799799452,
                -0.516297971645348,  0.856408600544260,  0.000844623419127,
                 0.001556088935770, -0.000048127534326,  0.999998788134748;

    Eigen::Matrix3d r2000_der_log;
    r2000_der_log << -3.764904133986157E-05,  6.245029944245220E-05,  6.168655101856225E-08,
                    -6.245022671509238E-05, -3.764909170451489E-05,  9.537579702956507E-08,
                    -7.698574166572183E-11, -5.754292251350311E-11,  1.170274036610182E-13;

    // 3. Матрица VW (для KATH12M)
    Eigen::Matrix3d vw_kath;
    vw_kath << -0.6500920187740070, -0.7413626903719400, -0.1666185117290845,
                0.7181503155181931, -0.6711045830006520,  0.1840618455754499,
               -0.2482750318647766,  0.0000000000000000,  0.9686895831754072;

    // 4. Скорости
    Eigen::Matrix3d earth_vel; earth_vel.setZero();
    earth_vel.col(1) << -27740.31190, -11133.38704, -4825.09396;

    vector<Eigen::Vector3d> vsta_log(2);
    vsta_log[0] << -129.97431, -431.49437, 0.18113;
    vsta_log[1] << -225.01847, -339.06466, 0.33365;

    vector<Eigen::Matrix3d> r2000_v = { r2000_log.transpose(), r2000_der_log.transpose() };
    vector<Eigen::Matrix3d> vw_v = { vw_kath, vw_kath }; 

    ariadna::aber_source(obs, r2000_v, k_s_log, earth_vel, vsta_log, vw_v, e, az);

    cout << "\n================= FINAL VALIDATION (STATION 1) =================" << endl;
    print_result("Elevation (E11)", e(0,0), 0.1142547638, " rad");
    print_result("Azimuth (Az11)",   az(0,0), 4.3191487714, " rad");
    print_result("dE/dt (E12)",      e(0,1), -0.0000652463, " rad/s");
    print_result("dAz/dt (Az12)",    az(0,1), -0.0000149985, " rad/s");

    return 0;
}