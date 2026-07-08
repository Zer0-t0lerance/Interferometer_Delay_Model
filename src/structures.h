#pragma once
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#include "constants.h"
#include "..\\external\\eigen\\Eigen\\Dense"
#include "..\\external\\spline.h"

namespace ariadna {
// Описание: Данные океанического прилива для одной станции (11 волн, 3 компоненты (Up, North, East))
struct OceanTideData {
    // AMP_OCEAN[3, 11] - Амплитуды [m] в (Up, North, East) для 11 волн
    Eigen::Matrix<double, 3, cnst::NUM_TIDES> amplitudes;
    // PHA_OCEAN[3, 11] - Фазы [rad] в (Up, North, East) для 11 волн
    Eigen::Matrix<double, 3, cnst::NUM_TIDES> phases;
};

// Описание: Данные атмосферной нагрузки для одной станции (каталог VLBI_atmload4_12.cat).
// Модель: смещение = b0 + Sum_{i=1..3}[A_i*cos(omega_i*t) + B_i*sin(omega_i*t)]
//         + b1*(P - p_0), где omega = 2Pi/{365.2422, 182.12, 1.0} сут, t = MJD - 45700.
// Амплитуды хранятся В МИЛЛИМЕТРАХ (как в каталоге); перевод в метры — в SITE_ATM40.
struct AtmLoadData {
    // coef(c, k): c = 0(Vertical),1(East),2(North);
    //             k = 0..5 -> A1,B1,A2,B2,A3,B3;  6 -> b0 (const);  7 -> b1 (регрессия по давлению)
    Eigen::Matrix<double, 3, 8> coef;
    double p_0;      // Опорное (среднее) давление станции [мбар] — из каталога ANTENNA_INFO
    bool has_data;   // true, если станция найдена в каталоге атмнагрузки
    AtmLoadData() : p_0(0.0), has_data(false) { coef.setZero(); }
};

// Параметры антенны из каталога ANTENNA_INFO (Nothnagel 2009) для термодеформации
// и опорных метеоусловий. Читаются форматом READ_CAT42corr FORMAT 155.
// Значения по умолчанию — для мобильных антенн (DATA hf/1/, hp/5/, gamma_f/1e-5/, gamma_a/1.2e-5/).
struct DefPar {
    std::string name;
    double t_0 = 20.0;        // Опорная температура [°C]
    double p_0 = 0.0;         // Опорное давление [мбар]
    double ant_diam = 0.0;    // Диаметр зеркала [м]
    double hf = 1.0, gamma_hf = 1.0e-5;   // Высота фундамента [м] + коэфф. расширения [1/°C]
    double hp = 5.0, gamma_hp = 1.2e-5;   // Высота опоры/колонны [м] + коэфф. расширения
    double AO = 0.0, gamma_AO = 1.2e-5;   // Осевое смещение [м] + коэфф. (в THERM_DEF40 не используется)
    double hv = 0.0, gamma_hv = 1.2e-5;   // Высота вершины зеркала [м] (не используется)
    double hs = 0.0, gamma_hs = 1.2e-5;   // Высота субрефлектора [м] (не используется)
};

struct Station {
    std::string name; // Station name (e.g., "RASTRON" for space telescope)
    Eigen::Vector3d xyz; // Coordinates in ITRF (m)
    Eigen::Vector3d vel; // Velocities (m/year)
    std::string axsty; // Mount type (e.g., "AZEL")
    double offs; // Axis offset (m)
    double lat_geod; // Geodetic latitude (rad)
    double lon_geod; // Geodetic longitude (rad)
    double h_geod; // Geodetic height (m)
    std::string domes; // DOMES number
    std::string descr; // Description
    OceanTideData tide_data;
    AtmLoadData atm_load;
    DefPar def_par; // Параметры антенны (термодеформация, опорные метео) из ANTENNA_INFO
};

struct Source {
    std::string name; // Source name
    std::string icrf_name; // ICRF name
    double ra; // Right ascension (rad)
    double dec; // Declination (rad)
    double ra_rate; // RA rate (rad/s)
    double dec_rate; // Dec rate (rad/s)
};

struct Observation {
    int mjd; // Modified Julian Date
    double utc; // UTC time (fraction of day)
    int sta1; // Index of station 1
    int sta2; // Index of station 2
    int sou; // Index of source
    double tau; // Delay (s)
    double dtau; // Delay error (s)
    double f; // Fringe frequency (Hz)
    double df; // Fringe frequency error (Hz)
    double cab1; // Cable correction for station 1 (s)
    double cab2; // Cable correction for station 2 (s)
    double t1; // Temperature at station 1 (C)
    double p1; // Pressure at station 1 (mb)
    double e1; // Humidity parameter at station 1
    double t2; // Temperature at station 2 (C)
    double p2; // Pressure at station 2 (mb)
    double e2; // Humidity parameter at station 2
    int kw1; // Humidity flag for station 1 (0=rel.hum, 1=dew point, 2=wet bulb)
    int kw2; // Humidity flag for station 2
};

struct SpaceStation {
    int mjd; // Modified Julian Date
    double utc; // UTC time (fraction of day)
    Eigen::Vector3d xyz; // Position in GCRS (km)
    Eigen::Vector3d vel; // Velocity in GCRS (km/s)
    Eigen::Vector3d acc; // Acceleration in GCRS (km/s^2)
};

struct OrbitData {
    int mjd; // Modified Julian Date
    double utc; // UTC time (fraction of day)
    double fres; // Doppler frequency (Hz)
    double phi; // Doppler phase (rad)
    double cor; // Corrected phase (rad)
};

struct AberrationResult {
    Eigen::MatrixXd elevation; // Elevation angles for two stations (rad, rad/s)
    Eigen::MatrixXd azimuth;   // Azimuth angles for two stations (rad, rad/s)
};

struct EOPData {
    double mjd;           // Modified Julian Date
    double ut1_utc;       // UT1-UTC (seconds)
    double ut1_tai;       // UT1-TAI (seconds)
    double x;             // Polar motion x (arcseconds)
    double y;             // Polar motion y (arcseconds)
    double dpsi;          // Nutation angle dpsi (arcseconds)
    double deps;          // Nutation angle deps (arcseconds)
};

// Подготовленные данные одной станции для посуточной сборки координат (site_pair):
// координаты ITRF на эпоху + геодезия/матрица VEN->ITRF + каталоги нагрузок + метео.
struct SitePrep {
    std::string name;
    Eigen::Vector3d xsta_itrf = Eigen::Vector3d::Zero(); // координаты ITRF на эпоху [м]
    double lat_gcen = 0.0, lon_gcen = 0.0, lat_geod = 0.0, h_geod = 0.0;
    Eigen::Matrix3d vw_i = Eigen::Matrix3d::Identity();   // VEN -> ITRF
    OceanTideData tide_data;    // океаническая нагрузка
    AtmLoadData atm_load;       // атмосферная нагрузка
    double pres = 0.0, dPdt = 0.0; // давление станции и его производная
    double tC = 0.0, dTdt = 0.0;   // температура станции [°C] и её производная [°C/с] (для термодеформации)
    DefPar def_par;             // параметры антенны (термодеформация) из ANTENNA_INFO
    std::string axsty = "AZEL"; // тип монтировки (для mount_tel)
    double offs = 0.0;          // смещение оси монтировки [м]
};

struct SiteCorrData {
    Eigen::Vector3d xyz;     // Координаты на эпоху наблюдения ITRF [м]
    Eigen::Vector3d vel;     // Скорости станции ITRF [м/год]
    double lat_geod;         // Геодезическая широта [рад]
    double lon_geod;         // Геодезическая долгота [рад]
    double h_geod;           // Геодезическая высота над эллипсоидом [м]
    double u_site;           // Сферический радиус экваториальный (U) [м]
    double v_site;           // Сферический радиус полярный (V) [м]
    Eigen::Matrix3d vw_i;    // Матрица перехода VEN (Vertical, East, North) -> ITRF
};
} // namespace ariadna