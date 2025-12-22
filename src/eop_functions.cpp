#include "functions.h"
#include "structures.h"
#include "constants.h"
#include <vector>

// Заглушка для interp_eop40
void interp_eop40(int k_int, double mjd, double utc, double tt, double& ut1,
                  Eigen::VectorXd& eop_int, Eigen::VectorXd& deop_int,
                  Eigen::VectorXd& arg_oc_tide, Eigen::MatrixXd& deop_diu,
                  Eigen::MatrixXd& deop_lib, const std::vector<EOPData>& eop_data) {
    // Простая интерполяция - взять ближайшее значение
    // В реальности нужно реализовать сплайн или полином интерполяцию
    if (!eop_data.empty()) {
        const auto& eop = eop_data[0]; // Заглушка
        ut1 = eop.ut1_utc;
        eop_int << eop.x, eop.y, eop.ut1_utc, eop.dpsi, eop.deps, 0.0, 0.0; // 7 элементов
        deop_int << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // Заглушка
        arg_oc_tide.setZero();
        deop_diu.setZero();
        deop_lib.setZero();
    }
}
