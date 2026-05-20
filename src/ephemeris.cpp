#include "functions.h"
#include <stdexcept>
#include <iostream>

#include "../external/dephem-master/include/dephem/EphemerisRelease.h"
#include "../external/dephem-master/include/dephem/help.h"

namespace ariadna {

static dph::EphemerisRelease* eph_release = nullptr;

void init_ephemeris(const std::string& eph_file) {
    if (eph_release != nullptr) {
        delete eph_release;
    }
    
    eph_release = new dph::EphemerisRelease(eph_file);
    
    if (!eph_release->isReady()) {
        throw std::runtime_error("Ephemeris load error: " + eph_file);
    }
    
    std::cout << "JPL Ephemeris loaded successfully!\n" 
              << eph_release->releaseLabel() << "\n";
}

// Численное дифференцирование скорости для получения ускорения
static Eigen::Vector3d get_acceleration(unsigned target, unsigned center, double JED) {
    double dt_days = 1.0 / 86400.0; // Шаг ровно 1 секунда
    double state_plus[6], state_minus[6];

    eph_release->calculateBody(dph::Calculate::STATE, target, center, JED + dt_days, state_plus);
    eph_release->calculateBody(dph::Calculate::STATE, target, center, JED - dt_days, state_minus);

    Eigen::Vector3d a;
    for (int i = 0; i < 3; ++i) {
        // DEPHEM возвращает скорость в км/сек. Разница по времени = 2 секунды.
        double a_km_s2 = (state_plus[i + 3] - state_minus[i + 3]) / 2.0;
        // Перевод км/сек^2 -> м/сек^2
        a(i) = a_km_s2 * 1000.0;
    }
    return a;
}

void get_celestial_bodies(double jd, double ct, Eigen::Matrix3d& Earth, Eigen::Matrix3d& Sun, Eigen::Matrix3d& Moon) {
    if (eph_release == nullptr || !eph_release->isReady()) {
        throw std::runtime_error("Ephemeris not initialized. Call init_ephemeris first.");
    }

    double JED = jd + ct; 
    unsigned CENTER = 12; // 12 - Барицентр Солнечной системы (SSB)
    double state_arr[6];

    auto fill_matrix = [&](unsigned body_id, Eigen::Matrix3d& matrix) {
        eph_release->calculateBody(dph::Calculate::STATE, body_id, CENTER, JED, state_arr);
        
        for(int i = 0; i < 3; ++i) {
            matrix(i, 0) = state_arr[i] * 1000.0;     // Позиция: км -> м
            matrix(i, 1) = state_arr[i + 3] * 1000.0; // Скорость: км/сек -> м/сек
        }
        matrix.col(2) = get_acceleration(body_id, CENTER, JED);
    };

    fill_matrix(dph::Body::EARTH, Earth);
    fill_matrix(dph::Body::SUN, Sun);
    fill_matrix(dph::Body::MOON, Moon);
}

} // namespace ariadna