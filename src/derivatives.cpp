#include "functions.h"

void der_star(double jd, double ct, double dyear, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, Eigen::MatrixXd& dstar, Eigen::MatrixXd& dstar_rate) {
    dstar.setZero();
    dstar_rate.setZero();
}

void der_site(double dyear, const Eigen::MatrixXd& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, std::vector<Eigen::MatrixXd>& dsite, std::vector<Eigen::MatrixXd>& dsite_v) {
    dsite[0].setZero();
    dsite[1].setZero();
    dsite_v[0].setZero();
    dsite_v[1].setZero();
}
