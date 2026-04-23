#include "functions.h"
#include <cmath>

namespace ariadna {

void SITE_CORR(const Station& station, double dYear, SiteCorrData& out) {
    
    // 1. Проверка на космический телескоп
    // Космическим аппаратам не нужна тектоника плит и земная геодезия
    if (station.name == "RADIOASTRON" || station.name == "SPACE_TEL") { // Можно заменить на флаг station.is_space
        out.xyz = station.xyz;
        out.vel = station.vel;
        out.lat_geod = 0.0;
        out.lon_geod = 0.0;
        out.h_geod   = 0.0;
        out.u_site   = 0.0;
        out.v_site   = 0.0;
        out.vw_i.setIdentity();
        return;
    }

    // 2. Учет движения тектонических плит (пересчет на эпоху наблюдения)
    out.xyz = station.xyz + station.vel * dYear;
    out.vel = station.vel;

    // 3. Экваториальный и полярный радиусы
    out.u_site = std::sqrt(out.xyz(0) * out.xyz(0) + out.xyz(1) * out.xyz(1));
    out.v_site = out.xyz(2);

    // 4. Вычисление геодезической широты и высоты (через уже имеющуюся функцию GEOD)
    GEOD(out.u_site, out.v_site, out.lat_geod, out.h_geod);

    // Геодезическая долгота (азимут в плоскости XY)
    out.lon_geod = std::atan2(out.xyz(1), out.xyz(0));

    // 5. Расчет матрицы перехода из топоцентрической системы (VEN: Vertical, East, North) в ITRF
    double cos_lat = std::cos(out.lat_geod);
    double sin_lat = std::sin(out.lat_geod);
    double cos_lon = std::cos(out.lon_geod);
    double sin_lon = std::sin(out.lon_geod);

    // Столбец 0: Up (Вертикаль)
    out.vw_i(0, 0) = cos_lat * cos_lon;
    out.vw_i(1, 0) = cos_lat * sin_lon;
    out.vw_i(2, 0) = sin_lat;

    // Столбец 1: East (Восток)
    out.vw_i(0, 1) = -sin_lon;
    out.vw_i(1, 1) = cos_lon;
    out.vw_i(2, 1) = 0.0;

    // Столбец 2: North (Север)
    out.vw_i(0, 2) = -sin_lat * cos_lon;
    out.vw_i(1, 2) = -sin_lat * sin_lon;
    out.vw_i(2, 2) = cos_lat;
}

} // namespace ariadna