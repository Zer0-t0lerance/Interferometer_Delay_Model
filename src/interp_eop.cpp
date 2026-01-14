#include "functions.h"

namespace ariadna {

void interp_eop(int k_int, const Observation& obs, double tt, double& ut1,
                Eigen::VectorXd& eop_int, Eigen::VectorXd& deop_int,
                Eigen::VectorXd& arg_oc_tide, Eigen::MatrixXd& deop_diu,
                Eigen::MatrixXd& deop_lib, const std::vector<EOPData>& eop_data) {
    
    const int N = 7;
    // X_int как в Фортране: разница в целых днях + дробная часть UTC текущего дня
    double x_target = static_cast<double>(obs.mjd - eop_data[0].mjd) + obs.utc;
    
    std::vector<double> x_nodes(N);
    std::vector<std::vector<double>> y_smooth(5, std::vector<double>(N));

    // ПОДГОТОВКА УЗЛОВ (как в Фортране перед CUBSPL_COEF)
    for (int i = 0; i < N; ++i) {
        x_nodes[i] = static_cast<double>(i);
        
        // Считаем приливы строго на ПОЛНОЧЬ каждого дня (UTC = 0, TT = 0)
        double cent_i, dut_i, dlod_i, domega_i;
        Eigen::VectorXd f_i(5), fd_i(5);
        double jd_node = eop_data[i].mjd + 2400000.5;
        
        fund_arg(jd_node, 0.0, cent_i, f_i, fd_i);
        ut1r_2010(f_i, dut_i, dlod_i, domega_i);

        // В Фортране XDATA (FWORK) - это очищенные значения
        y_smooth[0][i] = eop_data[i].ut1_tai - dut_i; // UT1-TAI без приливов
        y_smooth[1][i] = eop_data[i].x;
        y_smooth[2][i] = eop_data[i].y;
        y_smooth[3][i] = eop_data[i].dpsi;
        y_smooth[4][i] = eop_data[i].deps;
    }

    // ИНТЕРПОЛЯЦИЯ ГЛАДКИХ ДАННЫХ
    for (int j = 0; j < 5; ++j) {
        tk::spline s;
        s.set_points(x_nodes, y_smooth[j]);
        eop_int(j) = s(x_target);
        // Производная по дням (в Фортране это ANS1)
        deop_int(j) = s.deriv(1, x_target); 
    }

    // РАСЧЕТ ПОПРАВОК НА МОМЕНТ НАБЛЮДЕНИЯ (obs.mjd + obs.utc + tt)
    double cent, dut, dlod, domega;
    Eigen::VectorXd f(5), fd(5);
    double jd0 = static_cast<double>(obs.mjd) + 2400000.5;
    
    fund_arg(jd0, tt, cent, f, fd);
    ut1r_2010(f, dut, dlod, domega);
    terms_71(cent, f, fd, deop_diu, arg_oc_tide);
    terms_lib(cent, f, fd, deop_lib);

    // СБОРКА ФИНАЛЬНОГО UT1-TAI (F_int)
    // Добавляем зональные приливы, суточные приливы и либрации
    eop_int(0) += dut + deop_diu(0, 0) + deop_lib(0, 0);
    
    // СБОРКА ПРОИЗВОДНОЙ (DF_int)
    // domega (из ut1r) в рад/с переводим в с/с через CTIMRAD
    // Производную сплайна deop_int(0) (сек/день) переводим в сек/сек
    deop_int(0) = (domega / cnst::CTIMRAD) + deop_diu(0, 1) + deop_lib(0, 1) + (deop_int(0) / 86400.0);
    
    // Для остальных параметров (x, y, dpsi, deps) производную просто переводим в сек/сек
    for(int j = 1; j < 5; ++j) {
        deop_int(j) /= 86400.0;
    }

    // Итоговое значение UT1
    ut1 = obs.utc + eop_int(0) / 86400.0;
}
}