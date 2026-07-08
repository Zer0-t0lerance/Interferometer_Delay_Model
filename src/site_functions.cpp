#include "functions.h"

namespace ariadna {

// GEOD (декартовы -> геодезические, Borkowski 1989) реализован в GEOD.cpp,
// объявлен в functions.h. Здесь дубликат удалён — линковать с GEOD.cpp.

// Функция для вращения вокруг оси (аналог R_123)
void rotation_matrix(int axis, double angle, Eigen::Matrix3d& R) {
    R.setIdentity();
    double c = std::cos(angle);
    double s = std::sin(angle);
    if (axis == 1) { // X
        R(1, 1) = c; R(1, 2) = -s;
        R(2, 1) = s; R(2, 2) = c;
    } else if (axis == 2) { // Y
        R(0, 0) = c; R(0, 2) = s;
        R(2, 0) = -s; R(2, 2) = c;
    } else if (axis == 3) { // Z
        R(0, 0) = c; R(0, 1) = -s;
        R(1, 0) = s; R(1, 1) = c;
    }
}

// Функция для умножения матриц (ARRAY_2)
Eigen::Matrix3d matrix_multiply(const Eigen::Matrix3d& A, const Eigen::Matrix3d& B) {
    return A * B;
}

// Основная функция site
void site(const std::vector<Station>& stations, const std::vector<SpaceStation>& space_stations, int n_stations, double dyear,
          std::vector<Eigen::Vector3d>& site_xyz, std::vector<Eigen::Vector3d>& site_vel, std::vector<double>& lat_geod,
          std::vector<double>& h_geod, std::vector<double>& lat_gcen,
          std::vector<double>& lon_gcen, std::vector<double>& sph_rad,
          std::vector<double>& u_site, std::vector<double>& v_site,
          std::vector<Eigen::Matrix3d>& vw) {
    site_xyz.resize(n_stations);
    site_vel.resize(n_stations);
    lat_geod.resize(n_stations);
    h_geod.resize(n_stations);
    lat_gcen.resize(n_stations);
    lon_gcen.resize(n_stations);
    sph_rad.resize(n_stations);
    u_site.resize(n_stations);
    v_site.resize(n_stations);
    vw.resize(n_stations);

    for (int i = 0; i < n_stations; ++i) {
        const auto& sta = stations[i];

        // Пропустить центр Земли
        if (sta.name == "CENTER") continue;

        // Координаты и скорости на эпоху наблюдения
        site_xyz[i] = sta.xyz + sta.vel * dyear;
        site_vel[i] = sta.vel;

        // Для космического телескопа использовать данные из space_stations
        if (sta.name == "RASTRON" && !space_stations.empty()) {
            // Предполагаем, что space_stations[0] - текущие данные
            const auto& sp = space_stations[0];
            site_xyz[i] << sp.xyz(0) * 1000.0, sp.xyz(1) * 1000.0, sp.xyz(2) * 1000.0; // km to m
            site_vel[i] << sp.vel(0) * 1000.0, sp.vel(1) * 1000.0, sp.vel(2) * 1000.0; // km/s to m/s
        }

        // Сферический радиус
        sph_rad[i] = site_xyz[i].norm();

        // Геодентрическая широта и долгота
        lat_gcen[i] = std::asin(site_xyz[i](2) / sph_rad[i]);
        lon_gcen[i] = std::atan2(site_xyz[i](1), site_xyz[i](0));
        if (lon_gcen[i] < 0.0) lon_gcen[i] += cnst::TWOPI;

        // Расстояние от оси вращения и экваториальной плоскости
        u_site[i] = std::sqrt(site_xyz[i](0) * site_xyz[i](0) + site_xyz[i](1) * site_xyz[i](1)) * 1e-3; // km
        v_site[i] = site_xyz[i](2) * 1e-3; // km

        // Геодезические координаты (только для наземных станций)
        if (sta.name != "RASTRON") {
            double req = std::sqrt(site_xyz[i](0) * site_xyz[i](0) + site_xyz[i](1) * site_xyz[i](1));
            GEOD(req, site_xyz[i](2), lat_geod[i], h_geod[i]);

            // Матрица преобразования VEN (Up, East, North) -> ITRF.
            // Строим напрямую (прежняя композиция rotation_matrix давала неверные знаки Up).
            double cla = std::cos(lat_geod[i]), sla = std::sin(lat_geod[i]);
            double clo = std::cos(lon_gcen[i]), slo = std::sin(lon_gcen[i]);
            vw[i].col(0) << cla * clo, cla * slo, sla;        // Up (Vertical)
            vw[i].col(1) << -slo, clo, 0.0;                   // East
            vw[i].col(2) << -sla * clo, -sla * slo, cla;      // North
        } else {
            // Для космического телескопа заглушки
            lat_geod[i] = 0.0;
            h_geod[i] = 0.0;
            vw[i].setIdentity();
        }
    }
}

// Функция для расчета векторов источников
void source_vec(const std::vector<Source>& sources, double t_mean, std::vector<Eigen::Vector3d>& k_star) {
    k_star.resize(sources.size());
    // delta_t - разница времени от J2000.0 в годах
    double delta_t = (t_mean - (cnst::JD2000 - 2400000.5)) / 365.25;

    for (size_t i = 0; i < sources.size(); ++i) {
        const auto& sou = sources[i];
        double ra0 = sou.ra;
        double dec0 = sou.dec;
        double dot_RA = sou.ra_rate;
        double dot_DEC = sou.dec_rate;

        // RA и DEC на время t_mean
        double ra = ra0 + dot_RA * 1e-6 * delta_t * cnst::CARCRAD;
        double dec = dec0 + dot_DEC * 1e-6 * delta_t * cnst::CARCRAD;

        // Единичный вектор
        k_star[i](0) = std::cos(dec) * std::cos(ra);
        k_star[i](1) = std::cos(dec) * std::sin(ra);
        k_star[i](2) = std::sin(dec);
    }
}

} // namespace ariadna