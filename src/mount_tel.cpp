#include "functions.h"
#include <stdexcept>
#include <cmath>

namespace ariadna {

void mount_tel(const Observation& obs, const Eigen::MatrixXd& r2000, const std::vector<Station>& stations,
               const std::vector<Eigen::Vector3d>& k_star, const std::vector<Eigen::Matrix3d>& vw,
               const Eigen::MatrixXd& e, const Eigen::MatrixXd& az, Eigen::MatrixXd& doff_dl,
               Eigen::MatrixXd& d_dax, Eigen::MatrixXd& dtau_off) {
    
    // Инициализация выходных матриц нулями
    doff_dl.setZero(2, 2);
    d_dax.setZero(2, 2);
    dtau_off.setZero(2, 2);

    if (r2000.rows() != 9 || r2000.cols() != 3) {
        throw std::runtime_error("r2000 must be a 9x3 matrix");
    }

    // Транспонированные матрицы r2000
    Eigen::Matrix3d rt2000[2];
    rt2000[0] = r2000.block<3, 3>(0, 0).transpose();
    rt2000[1] = r2000.block<3, 3>(3, 0).transpose();

    for (int j = 0; j < 2; ++j) {
        int sta_idx = (j == 0) ? obs.sta1 : obs.sta2;
        const Station& station = stations[sta_idx];

        // Космические станции пропускаем
        if (station.name == "RASTRON") continue;

        const std::string& axtype = station.axsty;
        double offset = station.offs;
        double lat_gd = station.lat_geod;
        double t_c = (j == 0) ? obs.t1 : obs.t2;
        double p_mb = (j == 0) ? obs.p1 : obs.p2;
        double hum = (j == 0) ? obs.e1 : obs.e2;

        // Вектор фиксированной оси монтировки в топоцентрической системе VEN
        Eigen::Vector3d unit_I;
        if (axtype == "AZEL") {
            unit_I << 1.0, 0.0, 0.0;
        } else if (axtype == "EQUA") {
            unit_I << std::sin(lat_gd), 0.0, std::cos(lat_gd);
        } else if (axtype == "X-Y1" || axtype == "X-YN") {
            unit_I << 0.0, 0.0, 1.0;
        } else if (axtype == "X-Y2" || axtype == "X-YE") {
            unit_I << 0.0, 1.0, 0.0;
        } else if (station.name == "RICHMOND") {
            double w1 = cnst::RICHM_LAT * cnst::CDEGRAD;
            double az_rad = cnst::RICHM_AZ * cnst::CDEGRAD;
            unit_I << std::sin(w1), -std::cos(w1) * std::sin(az_rad), std::cos(w1) * std::cos(az_rad);
        } else {
            throw std::runtime_error("Invalid antenna axis type: " + axtype);
        }

        Eigen::Matrix3d vw_tr = vw[sta_idx].transpose();

        // Атмосферная рефракция
        double z = cnst::HALFPI - e(j, 0);
        double az_temp = az(j, 0);
        double temp_k = t_c + cnst::KELVIN_OFFSET;
        double humid_f = hum / cnst::PERCENT_TO_FRACTION;
        double press_hg = p_mb * cnst::MBAR_TO_MMHG;
        
        // Вызов функции sbend (угол изгиба луча)
        double rho = sbend(e(j, 0), temp_k, humid_f, press_hg);

        Eigen::Vector3d app;
        app << std::cos(z - rho), 
               std::sin(z - rho) * std::sin(az_temp), 
               std::sin(z - rho) * std::cos(az_temp);
        Eigen::Vector3d star_unit_app = app.normalized();

        // Расчет векторов задержки
        Eigen::Vector3d work1 = star_unit_app.cross(unit_I);
        Eigen::Vector3d vec_L = unit_I.cross(work1);
        double abs_vec_L = vec_L.norm();
        Eigen::Vector3d unit_vec_L = vec_L / abs_vec_L;

        Eigen::Vector3d unit_cff = vw[sta_idx] * unit_vec_L;
        Eigen::Vector3d unit_ax2000 = r2000.block<3, 3>(0, 0) * unit_cff;

        // Расчет производных
        Eigen::Vector3d work2 = vw[sta_idx] * star_unit_app;
        Eigen::Vector3d star_ab2000 = r2000.block<3, 3>(0, 0) * work2;
        Eigen::Vector3d dstar_ab_cff = rt2000[1] * star_ab2000;
        Eigen::Vector3d dstar_ab_ven = vw_tr * dstar_ab_cff;

        Eigen::Vector3d work3 = dstar_ab_ven.cross(unit_I);
        Eigen::Vector3d dvec_L = unit_I.cross(work3);
        double dabs_vec_L = vec_L.dot(dvec_L) / abs_vec_L;

        Eigen::Vector3d daxis_ven = (dvec_L / abs_vec_L) - (vec_L * dabs_vec_L / (abs_vec_L * abs_vec_L));
        Eigen::Vector3d daxis_cff = vw[sta_idx] * daxis_ven;
        
        Eigen::Vector3d work4 = r2000.block<3, 3>(0, 0) * daxis_cff;
        Eigen::Vector3d work5 = r2000.block<3, 3>(3, 0) * unit_cff;
        Eigen::Vector3d dunit_ax2000 = work4 + work5;

        // Финальные коэффициенты задержки монтировки
        double n_air = cnst::N_AIR_DEFAULT;
        doff_dl(j, 0) = star_ab2000.dot(unit_ax2000);
        doff_dl(j, 1) = star_ab2000.dot(dunit_ax2000);
        
        d_dax(j, 0) = (j == 0 ? 1.0 : -1.0) * (doff_dl(j, 0) / cnst::C) * n_air;
        d_dax(j, 1) = (j == 0 ? 1.0 : -1.0) * (doff_dl(j, 1) / cnst::C) * n_air;
        
        dtau_off(j, 0) = d_dax(j, 0) * offset;
        dtau_off(j, 1) = d_dax(j, 1) * offset;
    }
}

} // namespace ariadna