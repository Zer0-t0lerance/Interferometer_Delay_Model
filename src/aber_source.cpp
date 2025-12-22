// aber_source.cpp
#include "functions.h"
#include <iostream>

namespace ariadna {
// Helper function to compute aberrated source vector
static Eigen::Vector3d compute_aberrated_vector(const Eigen::Vector3d& k_s, const Eigen::Vector3d& v_total, double c) {
    // Total velocity vector norm
    double abs_vec_v = v_total.norm();
    
    // Normalize velocity vector
    Eigen::Vector3d unit_v = v_total / abs_vec_v;
    
    // Scalar product of source vector and velocity
    double v_dot_s = k_s.dot(v_total);
    double n_dot_s = k_s.dot(unit_v);
    
    // Relativistic factor
    double rel_p = 1.0 / std::sqrt(1.0 - std::pow(abs_vec_v / c, 2));
    double w1 = (rel_p - 1.0) / rel_p;
    double w2 = 1.0 / (1.0 + v_dot_s / c);
    
    // Equation (6-66) from the book
    Eigen::Vector3d k_star_aber = w2 * (k_s / rel_p + v_total / c + w1 * n_dot_s * unit_v);
    
    return k_star_aber.normalized(); // Return normalized aberrated vector
}

// Helper function to compute elevation and azimuth angles
static void compute_angles(const Eigen::Vector3d& ven_star, double& elevation, double& azimuth) {
    elevation = std::asin(ven_star(0));
    azimuth = std::atan2(ven_star(1), ven_star(2));
    if (azimuth < 0.0) azimuth += cnst::TWOPI;
}

// Helper function to compute time derivatives of angles
static void compute_angle_derivatives(const Eigen::Vector3d& ven_star, const Eigen::Vector3d& ven_star_der, double elevation, double& elevation_der, double& azimuth_der) {
    elevation_der = ven_star_der(0) / std::cos(elevation);
    
    double ven_ratio = ven_star(1) / ven_star(2);
    double ven_der_ratio = (ven_star_der(1) / ven_star(2)) - (ven_star(1) * ven_star_der(2) / std::pow(ven_star(2),2));
    azimuth_der = ven_der_ratio / (1.0 + std::pow(ven_ratio, 2));
}

void aber_source(const Observation& obs, const Eigen::MatrixXd& r2000, double lat_geod1, double lat_geod2, double h_geod1, double h_geod2, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& vsta_j2000t, const std::vector<Eigen::Matrix3d>& vw, double jd, double ct, Eigen::MatrixXd& e, Eigen::MatrixXd& az) {
    // Индексы станций из наблюдения
    int j1 = obs.sta1;
    int j2 = obs.sta2;

    // Проверка входных данных
    if (vsta_j2000t.size() < 2 || vw.size() < 2 || r2000.rows() < 6 || r2000.cols() < 3) {
        std::cerr << "Invalid input dimensions" << std::endl;
        return;
    }

    // Инициализация выходных массивов
    e.setZero(2, 2);
    az.setZero(2, 2);

    // Цикл по станциям
    for (int j = 0; j < 2; ++j) {
        int sta_idx = (j == 0) ? j1 : j2;

        // Суммарная скорость (Земля + станция)
        Eigen::Vector3d v_total = earth.col(1) + vsta_j2000t[j];
//        std::cout << "Station " << j + 1 << ": v_total = " << v_total.transpose() << std::endl;

        double vdot_s = k_s.dot(v_total);
        double abs_vec_v = v_total.norm();
//        std::cout << "vdot_s = " << vdot_s << ", abs_vec_v = " << abs_vec_v << std::endl;

        // Проверка на нулевую норму
        if (abs_vec_v < 1e-10) {
            std::cerr << "Station " << j + 1 << ": Zero velocity norm" << std::endl;
            return;
        }

        Eigen::Vector3d unit_v = v_total.normalized();
        double ndot_s = k_s.dot(unit_v);
        double rel_p = 1.0 / std::sqrt(1.0 - std::pow(abs_vec_v / cnst::C, 2));
        double w1 = (rel_p - 1.0) / rel_p;
        double w2 = 1.0 / (1.0 + vdot_s / cnst::C);

        Eigen::Vector3d k_star_aber = w2 * (k_s / rel_p + v_total / cnst::C + w1 * ndot_s * unit_v);
//        std::cout << "k_star_aber = " << k_star_aber.transpose() << std::endl;

        // Проверка на нулевую норму аберрированного вектора
        if (k_star_aber.norm() < 1e-10) {
            std::cerr << "Station " << j + 1 << ": Zero aberrated vector norm" << std::endl;
            return;
        }

        Eigen::Vector3d k_unit_aber = k_star_aber.normalized();
//        std::cout << "k_unit_aber = " << k_unit_aber.transpose() << std::endl;

        // Преобразование в земную систему координат
        Eigen::Matrix3d rt2000 = r2000.block(0, 0, 3, 3).transpose();
        Eigen::Vector3d crust_star = rt2000 * k_unit_aber;
//        std::cout << "crust_star = " << crust_star.transpose() << std::endl;

        // Преобразование в топоцентрическую систему
        Eigen::Matrix3d vw_tr = vw[sta_idx].transpose();
        Eigen::Vector3d ven_star = vw_tr * crust_star;
//        std::cout << "ven_star = " << ven_star.transpose() << std::endl;

        // Вычисление угла возвышения
        if (ven_star[0] < -1.0 || ven_star[0] > 1.0) {
            std::cerr << "Station " << j + 1 << ": Invalid ven_star[0] for asin: " << ven_star[0] << std::endl;
            return;
        }
        e(j, 0) = std::asin(ven_star[0]);

        // Вычисление азимута
        if (std::abs(ven_star[2]) < 1e-10 && std::abs(ven_star[1]) < 1e-10) {
            std::cerr << "Station " << j + 1 << ": Both ven_star[1] and ven_star[2] are zero" << std::endl;
            return;
        }
        az(j, 0) = std::atan2(ven_star[1], ven_star[2]);
        if (az(j, 0) < 0.0) {
            az(j, 0) += cnst::TWOPI;
        }

        // Вычисление производных по времени
        Eigen::Matrix3d rt2000_der = r2000.block(3, 0, 3, 3).transpose();
        Eigen::Vector3d crust_star_der = rt2000_der * k_unit_aber;
        Eigen::Vector3d ven_star_der = vw_tr * crust_star_der;
//        std::cout << "ven_star_der = " << ven_star_der.transpose() << std::endl;

        // Производная угла возвышения
        double cos_e = std::cos(e(j, 0));
        if (std::abs(cos_e) < 1e-10) {
            std::cerr << "Station " << j + 1 << ": cos(e) too small: " << cos_e << std::endl;
            return;
        }
        e(j, 1) = ven_star_der[0] / cos_e;

        // Производная азимута
        if (std::abs(ven_star[2]) < 1e-10) {
            std::cerr << "Station " << j + 1 << ": ven_star[2] too small: " << ven_star[2] << std::endl;
            return;
        }
        double ven_ratio = ven_star[1] / ven_star[2];
        double ven_der_ratio = ven_star_der[1] / ven_star[2] - ven_star[1] * ven_star_der[2] / (ven_star[2] * ven_star[2]);
        double denom = 1.0 + std::pow(ven_ratio, 2);
        if (std::abs(denom) < 1e-10) {
            std::cerr << "Station " << j + 1 << ": Denominator too small: " << denom << std::endl;
            return;
        }
        az(j, 1) = ven_der_ratio / denom;

//        std::cout << "Station " << j + 1 << ": e = " << e(j, 0) << ", e_rate = " << e(j, 1)
//                  << ", az = " << az(j, 0) << ", az_rate = " << az(j, 1) << std::endl;
    }
}
}