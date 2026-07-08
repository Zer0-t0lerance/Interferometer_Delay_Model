// process_ariadna.cpp
//
// Верхний слой оркестрации ARIADNA: пакетный расчёт теоретической задержки для
// всех наблюдений сессии. Соответствует основному циклу ARIADNA4_5corr (без слоя
// оценивания der_*/create_matr — он вне области переноса) и переиспользует уже
// реализованные и сверенные с ARIADNA блоки:
//
//   ОДИН РАЗ (на сессию):
//     * site()        — геодезия станций + матрицы VEN->ITRF + дрейф на эпоху;
//     * source_vec()  — единичные векторы направлений на источники в J2000.
//
//   НА КАЖДОЕ НАБЛЮДЕНИЕ:
//     * tai_time()/t_eph40()      — шкалы времени (TAI, TT, TDB);
//     * (окно 7 узлов EOP)+interp_iers() — x, y, UT1-UTC на момент;
//     * get_celestial_bodies()/jpl_eph() — эфемериды SSB и геоцентрические;
//     * get_r2000_matrices()/fund_arg()/gast_iau2006() — ориентация Земли;
//     * compute_delay_obs()        — весь конвейер задержки (site_pair->...->theor_delay).
//
// Результат (mjd, utc, станции, источник, tau, dtau) пишется в текстовый файл
// output_path. Космический VLBI (space_stations/orbit_data) и сегментация вывода
// (n_segm/k_ch_*/delta_sec) здесь не задействованы — заглушены как параметры.

#include "functions.h"
#include "constants.h"
#include <fstream>
#include <cstdio>
#include <algorithm>

namespace ariadna {
namespace {

// Выбор 7 узлов EOP вокруг момента наблюдения (mjd+utc) из полной таблицы EOPData,
// отсортированной по MJD. Возвращает окно, гарантированно охватывающее момент
// (для 4-точечного Лагранжа внутри interp_iers).
void window_eop(const std::vector<EOPData>& eop, double rjd, int n,
                std::vector<EOPData>& nodes) {
    const int m = static_cast<int>(eop.size());
    // Первый узел с mjd >= rjd (нижняя граница окна — на 3 узла левее).
    int idx = 0;
    while (idx < m && eop[idx].mjd < rjd) ++idx;
    int start = idx - n / 2 - 1;            // центрируем окно на момент
    if (start < 0) start = 0;
    if (start > m - n) start = m - n;
    nodes.assign(eop.begin() + start, eop.begin() + start + n);
}

} // anonymous namespace

void process_ariadna(const std::vector<Station>& stations, const std::vector<Source>& sources,
                     const std::vector<Observation>& observations,
                     const std::vector<SpaceStation>& space_stations,
                     const std::vector<OrbitData>& orbit_data,
                     int n_segm, int k_ch_c, int k_ch_z, double delta_sec,
                     const std::string& output_path,
                     const std::vector<EOPData>& eop_data,
                     double mjd_mean, double utc_mean, double t_mean) {
    (void)orbit_data; (void)n_segm; (void)k_ch_c; (void)k_ch_z; (void)delta_sec;

    const int n_stations = static_cast<int>(stations.size());

    // --- ОДИН РАЗ: геодезия станций на среднюю эпоху сессии ---
    const double dyear = (mjd_mean + utc_mean - cnst::MJD_J2000) / cnst::DAYS_PER_YEAR;
    std::vector<Eigen::Vector3d> site_xyz, site_vel;
    std::vector<double> lat_geod, h_geod, lat_gcen, lon_gcen, sph_rad, u_site, v_site;
    std::vector<Eigen::Matrix3d> vw;
    site(stations, space_stations, n_stations, dyear, site_xyz, site_vel, lat_geod, h_geod,
         lat_gcen, lon_gcen, sph_rad, u_site, v_site, vw);

    // Подготовка «постоянной» части SitePrep для каждой станции (метео — на наблюдение).
    std::vector<SitePrep> prep(n_stations);
    for (int i = 0; i < n_stations; ++i) {
        SitePrep& s = prep[i];
        s.name      = stations[i].name;
        s.xsta_itrf = site_xyz[i];
        s.lat_gcen  = lat_gcen[i];
        s.lon_gcen  = lon_gcen[i];
        s.lat_geod  = lat_geod[i];
        s.h_geod    = h_geod[i];
        s.vw_i      = vw[i];
        s.tide_data = stations[i].tide_data;
        s.atm_load  = stations[i].atm_load;
        s.axsty     = stations[i].axsty;
        s.offs      = stations[i].offs;
    }

    // --- ОДИН РАЗ: направления на источники ---
    std::vector<Eigen::Vector3d> k_star;
    source_vec(sources, t_mean, k_star);

    // --- Вывод ---
    std::ofstream out(output_path);
    if (!out) {
        std::fprintf(stderr, "process_ariadna: не удалось открыть %s для записи\n", output_path.c_str());
        return;
    }
    out.setf(std::ios::scientific);
    out.precision(15);
    out << "# ARIADNA theoretical delays\n";
    out << "# mjd utc sta1 sta2 sou tau[s] dtau[s/s]\n";

    // --- НА КАЖДОЕ НАБЛЮДЕНИЕ ---
    for (const Observation& obs : observations) {
        const int mjd = obs.mjd;
        const double utc = obs.utc;
        const double jd0 = cnst::MJD_OFFSET + static_cast<double>(mjd); // JD в 0h суток

        // Шкалы времени: TAI, TT (доли суток).
        double TAI, TT;
        tai_time(static_cast<double>(mjd), utc, TAI, TT);

        // EOP на момент: окно 7 узлов -> Лагранж + приливные поправки.
        std::vector<EOPData> nodes;
        window_eop(eop_data, static_cast<double>(mjd) + utc, cnst::EOP_NDATA, nodes);
        Eigen::VectorXd eop_int;
        interp_iers(obs, nodes, eop_int);
        const double ut1_utc = eop_int(0);                 // сек
        const double xp = eop_int(1) * cnst::CARCRAD;      // угл.сек -> рад
        const double yp = eop_int(2) * cnst::CARCRAD;

        // Производные величины времени.
        const double ut1_frac = utc + ut1_utc / cnst::SECDAY;
        const double ct = t_eph40(jd0, TT, ut1_frac, 0.0, 0.0, 0.0); // TDB, геоцентр
        const double jd_tt  = jd0 + TT;
        const double jd_ut1 = jd0 + ut1_frac;
        const double cent   = (jd_tt - cnst::JD2000) / cnst::JUL_CENT;
        const double ut1_sec = utc * cnst::SECDAY + ut1_utc;

        // Эфемериды: SSB (для theor_delay) и геоцентрические (для приливов).
        Eigen::Matrix3d Earth, Sun, Moon;
        get_celestial_bodies(jd0, ct, Earth, Sun, Moon);
        Eigen::Matrix3d Eg; Eigen::MatrixXd sunM, moonM;
        jpl_eph(jd0, ct, Eg, sunM, moonM);
        Eigen::Matrix<double, 3, 2> sun_geo = sunM, moon_geo = moonM;

        // Ориентация Земли: ITRF->J2000 + производные, фунд.аргументы, GAST.
        Eigen::Matrix3d R, dR, d2R;
        get_r2000_matrices(jd_tt, jd_ut1, xp, yp, R, dR, d2R);
        Eigen::VectorXd f(5), fd(5); double c2;
        fund_arg(jd0, ct, c2, f, fd);
        Eigen::Vector2d gast = gast_iau2006(jd_tt, jd_ut1);

        // Станции этого наблюдения: постоянная часть + метео из записи.
        SitePrep s1 = prep[obs.sta1];
        SitePrep s2 = prep[obs.sta2];
        s1.pres = obs.p1; s1.dPdt = 0.0;
        s2.pres = obs.p2; s2.dPdt = 0.0;

        // Полный конвейер задержки.
        double tau, dtau;
        compute_delay_obs(s1, s2, k_star[obs.sou], obs,
                          mjd, utc, jd0, ct, cent, ut1_sec, f, fd, gast,
                          Earth, Sun, Moon, sun_geo, moon_geo,
                          xp, yp, 0.0, 0.0, R, dR, d2R, tau, dtau);

        out << mjd << ' ' << utc << ' ' << obs.sta1 << ' ' << obs.sta2 << ' '
            << obs.sou << ' ' << tau << ' ' << dtau << '\n';
    }
}

} // namespace ariadna
