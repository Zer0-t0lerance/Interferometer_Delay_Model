// trop_delay.cpp
#include "functions.h"

namespace ariadna {

void trop_delay(const Observation& obs, double jd, double ct, const Station& sta1, const Station& sta2, const Eigen::MatrixXd& e, const Eigen::MatrixXd& az, Eigen::MatrixXd& datmc_d, Eigen::MatrixXd& datmc_w, Eigen::MatrixXd& datmp_hmf, Eigen::MatrixXd& datmp_wmf, Eigen::MatrixXd& dgrad_n, Eigen::MatrixXd& dgrad_e, Eigen::MatrixXd& zen_dry, Eigen::MatrixXd& zen_wet) {
    // Initialize output matrices (2x2 for two stations, delay and rate)
    datmc_d.setZero(2, 2);
    datmc_w.setZero(2, 2);
    datmp_hmf.setZero(2, 2);
    datmp_wmf.setZero(2, 2);
    dgrad_n.setZero(2, 2);
    dgrad_e.setZero(2, 2);
    zen_dry.setZero(2, 2);
    zen_wet.setZero(2, 2);

    // Station parameters
    const Station* stations[2] = {&sta1, &sta2};
    double pres[2] = {obs.p1, obs.p2};
    double tc[2] = {obs.t1, obs.t2};
    double humid[2] = {obs.e1, obs.e2};

    // Process each station
    for (int j = 0; j < 2; ++j) {
        // Skip space-based stations
        if (stations[j]->name == "CENTER" || stations[j]->name == "RASTRON") {
            continue;
        }

        // Compute mapping functions
        Eigen::Vector2d hmf, wmf;
        nhmf2(jd + ct, stations[j]->lat_geod, stations[j]->h_geod, e(j, 0), hmf);
        nwmf2(stations[j]->lat_geod, e(j, 0), wmf);

        // Set mapping function partials (reverse sign for station 1)
        datmp_hmf(j, 0) = (j == 0) ? -hmf(0) : hmf(0);
        datmp_hmf(j, 1) = (j == 0) ? -hmf(1) * e(j, 1) : hmf(1) * e(j, 1);
        datmp_wmf(j, 0) = (j == 0) ? -wmf(0) : wmf(0);
        datmp_wmf(j, 1) = (j == 0) ? -wmf(1) * e(j, 1) : wmf(1) * e(j, 1);

        // Compute atmosphere gradient partials (Chen and Herring)
        double denom = std::sin(e(j, 0)) * std::tan(e(j, 0)) + cnst::TROP_GRADIENT_COEFF;
        dgrad_n(j, 0) = std::cos(az(j, 0)) / denom;
        dgrad_e(j, 0) = std::sin(az(j, 0)) / denom;
        dgrad_n(j, 1) = -std::sin(az(j, 0)) / denom * az(j, 1) - std::cos(az(j, 0)) / (denom * denom) * e(j, 1) * std::tan(e(j, 0)) * (2.0 - std::tan(e(j, 0)));
        dgrad_e(j, 1) = std::cos(az(j, 0)) / denom * az(j, 1) - std::sin(az(j, 0)) / (denom * denom) * e(j, 1) * std::tan(e(j, 0)) * (2.0 - std::tan(e(j, 0)));

        // Compute zenith delays
        double z_d, dot_z_d, dz_ddh, z_w, dot_z_w;
        double dot_pres = 0.0, dot_hum = 0.0, dot_tc = 0.0; // Rates set to zero as in original
        double dxdh = cnst::PRES_DXDH_COEFF;
        double dpdh = cnst::PRES_REF * cnst::PRES_EXP * std::pow(1.0 - 6.5e-3 * stations[j]->h_geod / cnst::TEMP_REF, cnst::PRES_EXP - 1.0) * dxdh;

        sast_dry(pres[j], dot_pres, stations[j]->lat_geod, stations[j]->h_geod, dpdh, z_d, dot_z_d, dz_ddh);
        sast_wet(humid[j], tc[j], dot_hum, dot_tc, z_w, dot_z_w);

        // Store zenith delays
        zen_dry(j, 0) = z_d / cnst::C;
        zen_dry(j, 1) = dot_z_d / cnst::C;
        zen_wet(j, 0) = z_w / cnst::C;
        zen_wet(j, 1) = dot_z_w / cnst::C;

        // Compute tropospheric contributions
        datmc_d(j, 0) = datmp_hmf(j, 0) * zen_dry(j, 0);
        datmc_d(j, 1) = datmp_hmf(j, 1) * zen_dry(j, 0) + datmp_hmf(j, 0) * zen_dry(j, 1);
        datmc_w(j, 0) = datmp_wmf(j, 0) * zen_wet(j, 0);
        datmc_w(j, 1) = datmp_wmf(j, 1) * zen_wet(j, 0) + datmp_wmf(j, 0) * zen_wet(j, 1);
    }
}

} // namespace ariadna