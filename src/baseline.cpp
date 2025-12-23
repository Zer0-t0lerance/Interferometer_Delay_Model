#include "functions.h"

namespace ariadna {
void baseline(const Eigen::Matrix3d& r2000, const Eigen::MatrixXd& xsta_j2000t,
              const Eigen::MatrixXd& vsta_j2000t, const Eigen::MatrixXd& asta_j2000, 
              Eigen::MatrixXd& base_line, Eigen::Vector3d& b_cfs) {
    
    // Вектор базы: Станция 2 - Станция 1
    base_line.col(0) = xsta_j2000t.col(1) - xsta_j2000t.col(0); // Position
    base_line.col(1) = vsta_j2000t.col(1) - vsta_j2000t.col(0); // Velocity

    // b_cfs = R^T * b (перевод из J2000 в Crust-Fixed)
    b_cfs = r2000.transpose() * base_line.col(0);
}
}