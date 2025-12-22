#include "functions.h"
#include "constants.h"

void baseline(const Eigen::Matrix3d& r2000, const std::vector<Eigen::Vector3d>& xsta_j2000t, const std::vector<Eigen::Vector3d>& vsta_j2000t, std::vector<Eigen::Vector3d>& base_line, Eigen::Vector3d& b_cfs) {
    base_line[0] = xsta_j2000t[1] - xsta_j2000t[0];
    base_line[1] = vsta_j2000t[1] - vsta_j2000t[0];
    b_cfs.setZero();
}

void theor_delay(const std::vector<Eigen::Vector3d>& base_line, const std::vector<Eigen::Vector3d>& xsta_j2000t, const std::vector<Eigen::Vector3d>& vsta_j2000t, const std::vector<Eigen::Vector3d>& asta_j2000, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const Eigen::MatrixXd& sun, const Eigen::MatrixXd& moon, const Eigen::MatrixXd& datmc_d, const Eigen::MatrixXd& datmc_w, const Eigen::MatrixXd& dtau_off, double& t2_t1, double& dt2_t1) {
    // Упрощенная реализация теоретической задержки
    // Геометрическая задержка
    double geom_delay = base_line[0].dot(k_s) / cnst::C;
    // Релятивистские эффекты (упрощенные)
    double rel_delay = 0.0; // Заглушка для релятивистских эффектов
    // Тропосферная задержка
    double trop_delay = datmc_d(0, 0) + datmc_w(0, 0) + datmc_d(1, 0) + datmc_w(1, 0);
    // Коррекции
    double corr_delay = dtau_off(0, 0) + dtau_off(1, 0);

    t2_t1 = geom_delay + rel_delay + trop_delay + corr_delay;
    dt2_t1 = 0.0; // Упрощение
}
