#include "functions.h"

namespace ariadna {

void interp_eop(int k_int, const Observation& obs, double tt, double& ut1,
                Eigen::VectorXd& eop_int, Eigen::VectorXd& deop_int,
                Eigen::VectorXd& arg_oc_tide, Eigen::MatrixXd& deop_diu,
                Eigen::MatrixXd& deop_lib, const std::vector<EOPData>& eop_data) {
    
    const int N = 7;
    double x_target = static_cast<double>(obs.mjd - eop_data[0].mjd) + obs.utc;
    
    std::vector<double> x_nodes(N);
    std::vector<std::vector<double>> y_smooth(5, std::vector<double>(N));

    // TT-UTC для узлов: (TAI-UTC високосные + 32.184)/86400 на эпоху наблюдения.
    // Раньше было жёстко 68.184 (36 сек, только 2016) — исправлено на nsec (эпоха-независимо).
    double idelt; nsec(static_cast<double>(obs.mjd), idelt);
    double tt_node = (idelt + 32.184) / cnst::SECDAY;

    for (int i = 0; i < N; ++i) {
        x_nodes[i] = static_cast<double>(i);
        
        double cent_i, dut_i, dlod_i, domega_i;
        Eigen::VectorXd f_i(5), fd_i(5);
        double jd_node = eop_data[i].mjd + 2400000.5;
        
        fund_arg(jd_node, tt_node, cent_i, f_i, fd_i);
        ut1r_2010(f_i, dut_i, dlod_i, domega_i);

        // Фортран интерполирует (UT1-TAI) - DUT
        y_smooth[0][i] = eop_data[i].ut1_tai - dut_i; 
        y_smooth[1][i] = eop_data[i].x;
        y_smooth[2][i] = eop_data[i].y;
        y_smooth[3][i] = eop_data[i].dpsi;
        y_smooth[4][i] = eop_data[i].deps;
    }

    for (int j = 0; j < 5; ++j) {
        tk::spline s;
        s.set_points(x_nodes, y_smooth[j]);
        eop_int(j) = s(x_target); 
        deop_int(j) = s.deriv(1, x_target); // Пока в "ед/день"
    }

    double cent, dut, dlod, domega;
    Eigen::VectorXd f(5), fd(5);
    double jd0 = static_cast<double>(obs.mjd) + 2400000.5;
    
    fund_arg(jd0, tt, cent, f, fd);
    ut1r_2010(f, dut, dlod, domega);
    terms_71(cent, f, fd, deop_diu, arg_oc_tide);
    terms_lib(cent, f, fd, deop_lib);

    // СБОРКА ЗНАЧЕНИЙ (добавляем приливы и високосные сек idelt для перевода TAI -> UTC)
    eop_int(0) += dut + deop_diu(0, 0) + deop_lib(0, 0) + idelt; // Сразу UT1-UTC
    eop_int(1) += deop_diu(1, 0) + deop_lib(1, 0);       
    eop_int(2) += deop_diu(2, 0) + deop_lib(2, 0);       

    // СБОРКА ПРОИЗВОДНЫХ (Корректная физика в 1/СЕК, исправление бага Фортрана)
    deop_int(0) = (deop_int(0) / cnst::SECDAY) + (domega / cnst::CTIMRAD) + deop_diu(0, 1) + deop_lib(0, 1);
    deop_int(1) = (deop_int(1) / cnst::SECDAY) + deop_diu(1, 1) + deop_lib(1, 1);
    deop_int(2) = (deop_int(2) / cnst::SECDAY) + deop_diu(2, 1) + deop_lib(2, 1);
    deop_int(3) = (deop_int(3) / cnst::SECDAY);
    deop_int(4) = (deop_int(4) / cnst::SECDAY);

    // ПЕРЕВОД УГЛОВ В СИ (Радианы)
    for (int j = 1; j <= 4; ++j) {
        eop_int(j)  *= cnst::CARCRAD; 
        deop_int(j) *= cnst::CARCRAD; 
    }

    for (int i = 1; i <= 2; ++i) {
        deop_diu(i, 0) *= cnst::CARCRAD; 
        deop_diu(i, 1) *= cnst::CARCRAD; 
        deop_lib(i, 0) *= cnst::CARCRAD; 
        deop_lib(i, 1) *= cnst::CARCRAD; 
    }

    ut1 = obs.utc * cnst::SECDAY + eop_int(0);
}
}