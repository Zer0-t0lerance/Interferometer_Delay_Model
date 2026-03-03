#include "functions.h"
#include <cmath>

namespace ariadna {

void trop_delay(const Observation& obs, double jd, double ct, const Station& sta1, const Station& sta2, const Eigen::MatrixXd& e, const Eigen::MatrixXd& az, Eigen::MatrixXd& datmc_d, Eigen::MatrixXd& datmc_w, Eigen::MatrixXd& datmp_hmf, Eigen::MatrixXd& datmp_wmf, Eigen::MatrixXd& dgrad_n, Eigen::MatrixXd& dgrad_e, Eigen::MatrixXd& zen_dry, Eigen::MatrixXd& zen_wet) {
    
    datmc_d.setZero(2, 2); datmc_w.setZero(2, 2);
    datmp_hmf.setZero(2, 2); datmp_wmf.setZero(2, 2);
    dgrad_n.setZero(2, 2); dgrad_e.setZero(2, 2);
    zen_dry.setZero(2, 2); zen_wet.setZero(2, 2);

    const Station* stations[2] = {&sta1, &sta2};
    double pres[2]  = {obs.p1, obs.p2};
    double tc[2]    = {obs.t1, obs.t2};
    double humid[2] = {obs.e1, obs.e2};

    double epoch = jd + ct;

    for (int j = 0; j < 2; ++j) {
        if (stations[j]->name == "CENTER" || stations[j]->name == "RASTRON") continue;

        Eigen::Vector2d hmf, wmf;
        nhmf2(epoch, stations[j]->lat_geod, stations[j]->h_geod, e(j, 0), hmf);
        nwmf2(stations[j]->lat_geod, e(j, 0), wmf);

        // Для первой станции (j=0) картирующие функции берутся со знаком минус
        double sign = (j == 0) ? -1.0 : 1.0;

        datmp_hmf(j, 0) = sign * hmf(0);
        datmp_hmf(j, 1) = sign * hmf(1) * e(j, 1);
        datmp_wmf(j, 0) = sign * wmf(0);
        datmp_wmf(j, 1) = sign * wmf(1) * e(j, 1);

        // Градиенты (знак минус для первой станции НЕ применяется, судя по дампу)
        double denom = std::sin(e(j, 0)) * std::tan(e(j, 0)) + 0.0032;
        dgrad_n(j, 0) = std::cos(az(j, 0)) / denom;
        dgrad_e(j, 0) = std::sin(az(j, 0)) / denom;
        
        dgrad_n(j, 1) = -std::sin(az(j, 0)) / denom * az(j, 1) - std::cos(az(j, 0)) / (denom * denom) * e(j, 1) * std::tan(e(j, 0)) * (2.0 - std::tan(e(j, 0)));
        dgrad_e(j, 1) =  std::cos(az(j, 0)) / denom * az(j, 1) - std::sin(az(j, 0)) / (denom * denom) * e(j, 1) * std::tan(e(j, 0)) * (2.0 - std::tan(e(j, 0)));

        // Производные метеопараметров обнулены, так как они считаются в отдельной функции
        double dot_pres = 0.0, dot_hum = 0.0, dot_tc = 0.0;
        
        double x_val = 1.0 - (6.5e-3) * stations[j]->h_geod / 293.15;
        double dxdh = -6.5e-3 / 293.15;
        double dpdh = 1013.25 * 5.26 * std::pow(x_val, 4.26) * dxdh;

        double z_d, dot_z_d, dz_ddh, z_w, dot_z_w;
        sast_dry(pres[j], dot_pres, stations[j]->lat_geod, stations[j]->h_geod, dpdh, z_d, dot_z_d, dz_ddh);
        sast_wet(humid[j], tc[j], dot_hum, dot_tc, z_w, dot_z_w);

        // Перевод зенитных задержек в секунды (деление на скорость света)
        zen_dry(j, 0) = z_d / cnst::C;
        zen_dry(j, 1) = dot_z_d / cnst::C;
        zen_wet(j, 0) = z_w / cnst::C;
        zen_wet(j, 1) = dot_z_w / cnst::C;

        // Финальные задержки для МНК (умножаем на n_air, как в дампе)
        double n_air = 1.0 + (pres[j] / 1013.25) * 0.000283; // Приближение N_air из Фортрана
        
        datmc_d(j, 0) = datmp_hmf(j, 0) * zen_dry(j, 0) * n_air;
        datmc_d(j, 1) = (datmp_hmf(j, 1) * zen_dry(j, 0) + datmp_hmf(j, 0) * zen_dry(j, 1)) * n_air;
        
        datmc_w(j, 0) = datmp_wmf(j, 0) * zen_wet(j, 0) * n_air;
        datmc_w(j, 1) = (datmp_wmf(j, 1) * zen_wet(j, 0) + datmp_wmf(j, 0) * zen_wet(j, 1)) * n_air;
    }
}
} // namespace ariadna