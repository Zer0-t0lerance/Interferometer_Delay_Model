#include "functions.h"

namespace ariadna {

void interp_eop(int k_int, const Observation& obs, double tt, double& ut1,
                Eigen::VectorXd& eop_int, Eigen::VectorXd& deop_int,
                Eigen::VectorXd& arg_oc_tide, Eigen::MatrixXd& deop_diu,
                Eigen::MatrixXd& deop_lib, const std::vector<EOPData>& eop_data) {
    // Validate input
    if (eop_data.size() != cnst::EOP_NDATA) {
        throw std::runtime_error("EOP data size must be " + std::to_string(cnst::EOP_NDATA));
    }
    if (k_int != 0 && k_int != 1) {
        throw std::runtime_error("Invalid interpolation method: k_int must be 0 or 1");
    }

    // Initialize output arrays
    eop_int.setZero(5);    // UT1-UTC, x, y, dpsi, deps
    deop_int.setZero(5);   // Derivatives
    arg_oc_tide.setZero(8); // Tidal arguments
    deop_diu.setZero(3, 2); // Diurnal tidal corrections
    deop_lib.setZero(3, 2); // Libration corrections

    // Compute interpolation point
    double x_int = obs.mjd - eop_data[0].mjd + obs.utc;

    // Prepare data for interpolation
    std::vector<double> x_data(cnst::EOP_NDATA);
    std::vector<std::vector<double>> f_data(5, std::vector<double>(cnst::EOP_NDATA));
    for (int i = 0; i < cnst::EOP_NDATA; ++i) {
        x_data[i] = eop_data[i].mjd - eop_data[0].mjd;
        f_data[0][i] = eop_data[i].ut1_tai; // UT1-TAI
        f_data[1][i] = eop_data[i].x;       // x
        f_data[2][i] = eop_data[i].y;       // y
        f_data[3][i] = eop_data[i].dpsi;    // dpsi
        f_data[4][i] = eop_data[i].deps;    // deps
    }

    // Correct UT1-TAI for zonal tides
    Eigen::VectorXd f(5), fd(5);
    double cent, dut, dlod, domega;
    for (int i = 0; i < cnst::EOP_NDATA; ++i) {
        double jd_eop = eop_data[i].mjd + 2400000.5;
        double tai_eop, tt_eop;
        tai_time(eop_data[i].mjd, 0.0, tai_eop, tt_eop); // Using tai_time from functions.h
        fund_arg(jd_eop, tt_eop, cent, f, fd); // Placeholder: need implementation
        ut1r_2010(f, dut, dlod, domega);      // Placeholder: need implementation
        f_data[0][i] -= dut;                   // Remove zonal tide effect
    }

    // Interpolation
    std::vector<double> f_int(5), df_int(5);
    if (k_int == 0) {
        // Cubic spline interpolation using tk::spline
        for (int j = 0; j < 5; ++j) {
            tk::spline spline(x_data, f_data[j], tk::spline::cspline);
            f_int[j] = spline(x_int);
            df_int[j] = spline.deriv(1, x_int);
        }
    } else {
        // Polynomial fit using Eigen
        Eigen::MatrixXd X(cnst::EOP_NDATA, cnst::EOP_NDEG + 1);
        Eigen::VectorXd y(cnst::EOP_NDATA);
        Eigen::VectorXd coeffs(cnst::EOP_NDEG + 1);
        for (int j = 0; j < 5; ++j) {
            for (int i = 0; i < cnst::EOP_NDATA; ++i) {
                y(i) = f_data[j][i];
                for (int k = 0; k <= cnst::EOP_NDEG; ++k) {
                    X(i, k) = std::pow(x_data[i], k);
                }
            }
            coeffs = X.colPivHouseholderQr().solve(y);
            f_int[j] = 0.0;
            df_int[j] = 0.0;
            for (int k = 0; k <= cnst::EOP_NDEG; ++k) {
                f_int[j] += coeffs(k) * std::pow(x_int, k);
                if (k > 0) {
                    df_int[j] += k * coeffs(k) * std::pow(x_int, k - 1);
                }
            }
        }
    }

    // Restore tidal effects for UT1
    double jd0 = obs.mjd + 2400000.5;
    fund_arg(jd0, tt, cent, f, fd); // Placeholder: need implementation
    ut1r_2010(f, dut, dlod, domega); // Placeholder: need implementation
    terms_71(cent, f, fd, deop_diu, arg_oc_tide); // Placeholder: need implementation
    terms_lib(cent, f, fd, deop_lib); // Placeholder: need implementation
    f_int[0] += dut + deop_diu(0, 0) + deop_lib(0, 0);
    df_int[0] = domega / cnst::CTIMRAD + deop_diu(0, 1) + deop_lib(0, 1) + df_int[0] / cnst::SECDAY;

    // Add leap seconds
    double idelt;
    nsec(x_int + eop_data[0].mjd, idelt); // Placeholder: need implementation
    f_int[0] += idelt;

    // Store results
    for (int j = 0; j < 5; ++j) {
        eop_int(j) = f_int[j];
        deop_int(j) = (j == 0) ? df_int[j] : df_int[j] / cnst::SECDAY;
    }

    // Correct polar motion coordinates
    eop_int(1) += deop_diu(1, 0) + deop_lib(1, 0); // x
    eop_int(2) += deop_diu(2, 0) + deop_lib(2, 0); // y
    deop_int(1) += deop_diu(1, 1) + deop_lib(1, 1); // dx/dTAI
    deop_int(2) += deop_diu(2, 1) + deop_lib(2, 1); // dy/dTAI

    // Special case: no interpolation if UTC == 0 (preserved for compatibility)
    if (std::abs(obs.utc) < 1e-15) {
        eop_int(0) = eop_data[2].ut1_utc;
        for (int j = 1; j < 5; ++j) {
            eop_int(j) = f_data[j][2];
        }
    }

    // Compute UT1
    ut1 = eop_int(0) + obs.utc * cnst::SECDAY;
}

} // namespace ariadna