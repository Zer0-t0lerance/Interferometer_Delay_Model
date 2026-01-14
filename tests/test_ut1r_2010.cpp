#include <iostream>
#include <iomanip>
#include "../src/functions.h"

using namespace std;

int main() {
    cout << scientific << setprecision(16);

    // 1. Исходные данные из лога
    double jd = 2457400.5;
    double ct = 0.7891666666666666E-03;
    
    // 2. Эталоны из лога
    double exp_dut = -0.01003213994142547; // Это -0.1003...D-01
    double exp_domega = -0.2436398437388438E-12;

    // 3. Вычисляем фундаментальные аргументы
    double cent;
    Eigen::VectorXd f(5), fd(5);
    ariadna::fund_arg(jd, ct, cent, f, fd);

    // 4. Вычисляем UT1R
    double act_dut, act_lod, act_domega;
    ariadna::ut1r_2010(f, act_dut, act_lod, act_domega);

    cout << "================= UT1R FINAL VALIDATION =================\n";
    cout << "DUT (UT1-UT1R):\n";
    cout << "  Act: " << act_dut << endl;
    cout << "  Exp: " << exp_dut << endl;
    double diff_dut = abs(act_dut - exp_dut);
    cout << "  Diff: " << diff_dut << (diff_dut < 1e-8 ? " [ OK ]" : " [ FAIL ]") << endl;

    cout << "\nDOMEGA:\n";
    cout << "  Act: " << act_domega << endl;
    cout << "  Exp: " << exp_domega << endl;
    double diff_omg = abs(act_domega - exp_domega);
    cout << "  Diff: " << diff_omg << (diff_omg < 1e-18 ? " [ OK ]" : " [ FAIL ]") << endl;

    // Дополнительная проверка на f из лога (чтобы убрать влияние fund_arg)
    Eigen::VectorXd f_log(5);
    f_log << 0.1141344609264314D+07, 0.3136439828193560D+05, 0.5678739386129379D+06, 0.1468483807328939D+06, -0.6660943330937190D+06;
    double dut_pure, lod_pure, omg_pure;
    ariadna::ut1r_2010(f_log, dut_pure, lod_pure, omg_pure);
    
    cout << "\nPure UT1R check (using log f vectors):\n";
    cout << "  Act DUT: " << dut_pure << endl;
    cout << "  Diff:    " << abs(dut_pure - exp_dut) << (abs(dut_pure - exp_dut) < 1e-10 ? " [ OK ]" : " [ FAIL ]") << endl;

    return 0;
}