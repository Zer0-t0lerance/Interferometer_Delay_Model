// scf_reader.cpp
//
// Чтение орбиты космического телескопа (RadioAstron) из файла .scf — стандарт
// CCSDS OEM (Orbit Ephemeris Message): геоцентрические координаты в системе
// EME2000 (=J2000), время UTC. Строки данных после META_STOP:
//   YYYY-MM-DDThh:mm:ss.sss  X Y Z VX VY VZ   ([км], [км/с])
//
// Возвращает вектор точек SpaceStation (mjd, utc, xyz [км], vel [км/с]); ускорение
// не задаётся (его даёт orbit_interp через производную сплайна). Кадр EME2000 ≈ J2000
// берётся напрямую (без r2000) — так же, как в orbit_interp.

#include "functions.h"
#include <fstream>
#include <sstream>
#include <cstdio>

namespace ariadna {

// Григорианская дата -> MJD на 0h (алгоритм Fliegel-Van Flandern).
static int ymd_to_mjd(int y, int m, int d) {
    long a = (14 - m) / 12;
    long yy = y + 4800 - a;
    long mm = m + 12 * a - 3;
    long jdn = d + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045;
    return static_cast<int>(jdn - 2400001); // MJD на 0h суток
}

bool read_scf_orbit(const std::string& path, std::vector<SpaceStation>& orbit) {
    orbit.clear();
    std::ifstream f(path);
    if (!f) { std::fprintf(stderr, "read_scf_orbit: не открыть %s\n", path.c_str()); return false; }

    std::string line;
    bool in_data = false;
    while (std::getline(f, line)) {
        if (line.find("META_STOP") != std::string::npos) { in_data = true; continue; }
        if (!in_data) continue;
        if (line.find("META_START") != std::string::npos) { in_data = false; continue; } // на случай нескольких сегментов
        if (line.empty()) continue;
        // Строка данных начинается с года (цифра).
        if (!std::isdigit(static_cast<unsigned char>(line[0]))) continue;

        // Разбор: дата в формате YYYY-MM-DDThh:mm:ss.sss, затем 6 чисел.
        std::string date;
        std::istringstream ss(line);
        ss >> date;
        int Y, M, D, h, mi; double s;
        if (std::sscanf(date.c_str(), "%d-%d-%dT%d:%d:%lf", &Y, &M, &D, &h, &mi, &s) != 6) continue;

        SpaceStation p;
        p.mjd = ymd_to_mjd(Y, M, D);
        p.utc = (h * 3600.0 + mi * 60.0 + s) / cnst::SECDAY;
        ss >> p.xyz.x() >> p.xyz.y() >> p.xyz.z() >> p.vel.x() >> p.vel.y() >> p.vel.z();
        p.acc.setZero();
        orbit.push_back(p);
    }
    return !orbit.empty();
}

} // namespace ariadna
