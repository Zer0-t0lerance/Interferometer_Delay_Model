#include <iostream>
#include <iomanip>
#include "..//src//functions.h"

using namespace std;

int main() {
    cout << scientific << setprecision(15);

    // Данные для параметра J=3 из лога Фортрана
    vector<double> x_data = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    vector<double> f_data = {-35.9309746600586, -35.9327674362612, -35.9345265604009, 
                             -35.9362060723846, -35.9378630684868, -35.9395461403042, -35.9412323799315};
    double x_target = 2.27051504629344;

    // Ожидаемые значения из Фортрана
    double exp_f_int = -35.9349878467769;
    double exp_df_int = -1.692677602117368E-03;

    tk::spline s;
    s.set_points(x_data, f_data);
    
    double act_f_int = s(x_target);
    double act_df_int = s.deriv(1, x_target); 

    cout << "--- Spline Check ---" << endl;
    cout << "Value: " << act_f_int << " (Exp: " << exp_f_int << ")" << endl;
    cout << "Diff:  " << abs(act_f_int - exp_f_int) << endl;
    cout << "Deriv: " << act_df_int << " (Exp: " << exp_df_int << ")" << endl;
    cout << "Diff:  " << abs(act_df_int - exp_df_int) << endl;

    return 0;
}