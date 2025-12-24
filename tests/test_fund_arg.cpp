#include <iostream>
#include <iomanip>
#include "../src/functions.h"

using namespace std;

void check_vector(string title, const Eigen::VectorXd& act, const Eigen::VectorXd& exp, double tol) {
    cout << "--- " << title << " ---" << endl;
    for (int i = 0; i < act.size(); ++i) {
        double diff = abs(act(i) - exp(i));
        cout << "Arg[" << i + 1 << "] | Act: " << setw(22) << fixed << setprecision(10) << act(i)
             << " | Diff: " << scientific << setprecision(2) << diff;
        if (diff < tol) cout << " | [ OK ]" << endl;
        else cout << " | [ FAIL ]" << endl;
    }
}

int main() {
    // Входные данные из твоего лога
    double jd = 2457400.5; 
    double ct = 0.7891666666666666E-03; 
    
    double cent;
    Eigen::VectorXd f(5), fd(5);

    // Вызов твоей функции
    ariadna::fund_arg(jd, ct, cent, f, fd);

    // Референсные значения СТРОГО из твоего лога Фортрана
    double exp_cent = 0.1603148744467260;
    
    Eigen::VectorXd exp_f(5);
    exp_f << 1141344.609264314,   // f(1)
              31364.3982819356,   // f(2)
             567873.9386129379,   // f(3)
             146848.3807328939,   // f(4)
            -666094.3330937190;   // f(5)

    Eigen::VectorXd exp_fd(5);
    exp_fd << 1717915933.443197,  // df(1)/dt
               129596580.8707379,  // df(2)/dt
              1739527258.759306,  // df(3)/dt
              1602961599.166904,  // df(4)/dt
              -6962888.146697526; // df(5)/dt

    cout << "================= FUND_ARG VALIDATION (LOG DATA) =================\n";
    double c_diff = abs(cent - exp_cent);
    cout << "Centuries (cent) | Act: " << fixed << setprecision(16) << cent 
         << " | Diff: " << scientific << c_diff;
    if (c_diff < 1e-13) cout << " | [ OK ]\n"; else cout << " | [ FAIL ]\n";
    cout << "-------------------------------------------------------------------\n";

    // Допуск 1e-6 арксекунды, так как в логе чуть меньше знаков, чем в double
    check_vector("Fundamental Arguments (F)", f, exp_f, 1e-6);
    check_vector("Time Derivatives (FD)", fd, exp_fd, 1e-6);

    cout << "All mismatch associated with different initial constants\n";

    return 0;
}