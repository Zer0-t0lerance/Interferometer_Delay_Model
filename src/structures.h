#pragma once
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#include "constants.h"
#include "..\\external\\eigen\\Eigen\\Dense"
#include "..\\external\\spline.h"

namespace ariadna {
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

// Описание: Данные океанического прилива для одной станции (11 волн, 3 компоненты (Up, North, East))
struct OceanTideData {
    // AMP_OCEAN[3, 11] - Амплитуды [m] в (Up, North, East) для 11 волн
    Eigen::Matrix<double, 3, cnst::NUM_TIDES> amplitudes;
    // PHA_OCEAN[3, 11] - Фазы [rad] в (Up, North, East) для 11 волн
    Eigen::Matrix<double, 3, cnst::NUM_TIDES> phases;
};
} // namespace ariadna