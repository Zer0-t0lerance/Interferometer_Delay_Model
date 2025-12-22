#include "functions.h"
#include <stdexcept>

namespace ariadna {
void baseline(const Eigen::Matrix3d& r2000, const std::vector<Eigen::Vector3d>& xsta_j2000t,
              const std::vector<Eigen::Vector3d>& vsta_j2000t, std::vector<Eigen::Vector3d>& base_line,
              Eigen::Vector3d& b_cfs) {
    if (xsta_j2000t.size() != 2 || vsta_j2000t.size() != 2 || base_line.size() != 2) {
        throw std::invalid_argument("Input vectors must have size 2");
    }

    base_line[0] = xsta_j2000t[1] - xsta_j2000t[0]; // Position
    base_line[1] = vsta_j2000t[1] - vsta_j2000t[0]; // Velocity

    if (base_line[0].norm() < cnst::BASELINE_NORM_THRESHOLD) {
        throw std::invalid_argument("Baseline length too small (< " +
                                   std::to_string(cnst::BASELINE_NORM_THRESHOLD) + " m)");
    }

    b_cfs = r2000.transpose() * base_line[0];
}
}