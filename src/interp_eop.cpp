#include "functions.h"
#include <vector>

namespace ariadna {

void interp_eop(int k_int, const Observation& obs, double tt, double& ut1,
                Eigen::VectorXd& eop_int, Eigen::VectorXd& deop_int,
                Eigen::VectorXd& arg_oc_tide, Eigen::MatrixXd& deop_diu,
                Eigen::MatrixXd& deop_lib, const std::vector<EOPData>& eop_data) {
    
    const int N = 7; 
    std::vector<double> x_nodes(N);
    for (int i = 0; i < N; ++i) x_nodes[i] = static_cast<double>(i);

    // Точка интерполяции относительно начала окна данных
    double x_target = (obs.mjd - eop_data[0].mjd) + obs.utc;
    
    double cent, dut, dlod, domega;
    Eigen::VectorXd f(5), fd(5);
    double jd0 = obs.mjd + 2400000.5;
    
    // Вызовы согласно вашим исходникам
    fund_arg(jd0, tt, cent, f, fd);
    ut1r_2010(f, dut, dlod, domega);
    terms_71(cent, f, fd, deop_diu, arg_oc_tide);
    terms_lib(cent, f, fd, deop_lib);

    // Интерполяция 5 параметров
    for (int j = 0; j < 5; ++j) {
        std::vector<double> y_nodes(N);
        for (int i = 0; i < N; ++i) {
            if (j == 0)      y_nodes[i] = eop_data[i].ut1_tai;
            else if (j == 1) y_nodes[i] = eop_data[i].x;
            else if (j == 2) y_nodes[i] = eop_data[i].y;
            else if (j == 3) y_nodes[i] = eop_data[i].dpsi;
            else             y_nodes[i] = eop_data[i].deps;
        }

        tk::spline s;
        s.set_points(x_nodes, y_nodes);
        
        eop_int(j) = s(x_target);
        // Исправленный порядок: (порядок производной, точка)
        deop_int(j) = s.deriv(1, x_target); 
    }

    // Финальные поправки для UT1 (аналог Фортрана)
    eop_int(0) += dut + deop_diu(0, 0) + deop_lib(0, 0);
    deop_int(0) = domega / cnst::CTIMRAD + deop_diu(0, 1) + deop_lib(0, 1) + deop_int(0) / 86400.0;
    
    ut1 = obs.utc + eop_int(0) / 86400.0;
}

} // namespace ariadna