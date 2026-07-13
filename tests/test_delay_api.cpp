// test_delay_api.cpp
//
// Проверка адаптера compute_task_polys (полиномы всех станций В ПАМЯТИ, без файлов):
// считаем полиномы задержки и координат (u,v,w) для примера и сверяем P0 каждого блока с
// эталонами example/*.TXT / *_uvw.txt по совпадающему ВРЕМЕНИ старта (посканая модель ->
// у станции может быть меньше блоков, чем в непрерывном эталоне). Ожидаемо те же невязки,
// что у файлового пути (process_task): наземные 1-4 м, RASTRON задержка 1.2 м, uvw <1 м.

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif
#include "../src/delay_api.h"
#include "../src/functions.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

using namespace ariadna;

static std::string find_eph() {
    const char* c[] = {"external/dephem-master/linux_p1550p2650.440t", "../external/dephem-master/linux_p1550p2650.440t"};
    for (auto p : c) { std::ifstream f(p, std::ios::binary); if (f.good()) return p; } return c[0];
}
static long ymd2mjd(int y, int m, int d) {
    long a = (14 - m) / 12, yy = y + 4800 - a, mm = m + 12 * a - 3;
    long jdn = d + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045; return jdn - 2400001;
}
// Ключ момента: секунды от MJD-эпохи (целые), для сопоставления блоков по времени старта.
static long ref_key(const std::string& start) {
    int D, M, Y, h, mi, s;
    if (std::sscanf(start.c_str(), "start = %d/%d/%d %dh%dm%ds", &D, &M, &Y, &h, &mi, &s) != 6) return -1;
    return ymd2mjd(Y, M, D) * 86400L + h * 3600 + mi * 60 + s;
}
static long mem_key(int mjd, double utc) { return (long)mjd * 86400L + std::lround(utc * 86400.0); }
// Чтение эталона: ключ времени -> вектор коэффициентов P0 (1 для задержки, 3 для uvw).
static std::map<long, std::vector<double>> read_ref(const std::string& path) {
    std::map<long, std::vector<double>> m; std::ifstream f(path); std::string line; long key = -1;
    while (std::getline(f, line)) {
        if (line.rfind("start", 0) == 0) key = ref_key(line);
        else if (line.rfind("P0", 0) == 0 && key >= 0) {
            std::vector<double> v; std::stringstream ss(line.substr(line.find('=') + 1)); std::string tok;
            while (std::getline(ss, tok, ',')) v.push_back(std::stod(tok));
            m[key] = v;
        }
    }
    return m;
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    printf("=====================================================================\n");
    printf("  АДАПТЕР В ПАМЯТИ: compute_task_polys(cfx, орбита) -> полиномы всех станций\n");
    printf("---------------------------------------------------------------------\n");
    try { init_ephemeris(find_eph()); } catch (const std::exception& e) { printf("SKIP: %s\n", e.what()); return 2; }

    TaskPolys tp;
    bool ok_call = compute_task_polys(
        "example/GVLBI_RAKS01A0_L_20140423T130000_ASC_V1.cfx",
        "example/RA140423_1200_v02.scf",
        "external/catalogs/EOPC04_14_IAU2000_62-now.cat",
        60.0, 5, 6.0, /*with_tropo=*/false, tp);
    if (!ok_call || tp.delay.empty()) { printf("  FAIL: compute_task_polys не вернул результат\n"); return 1; }
    printf("  станций получено: %zu (delay+uvw в памяти, файлы не писались)\n\n", tp.delay.size());

    // Имя станции -> базовое имя эталонного файла.
    auto refbase = [](const std::string& tel) -> std::string {
        if (tel == "RASTRON") return "RA_L"; if (tel == "KALYAZIN") return "KALYAZIN_L";
        if (tel == "HARTRAO") return "HARTRAO_L"; return tel; // BADARY -> BADARY
    };
    const double C = 2.99792458e8;
    bool ok = true; double worst_d = 0, worst_u = 0;

    printf("  Сверка задержки (P0) и uvw (P0: u,v,w) с эталонами:\n");
    for (size_t i = 0; i < tp.delay.size(); ++i) {
        const std::string& tel = tp.delay[i].telescope; std::string base = refbase(tel);
        auto rd = read_ref("example/" + base + ".TXT");
        auto ru = read_ref("example/" + base + "_uvw.txt");
        double ed = 0, eu = 0; int nd = 0, nu = 0;
        for (const auto& b : tp.delay[i].blocks) { long k = mem_key(b.mjd, b.utc_start);
            auto it = rd.find(k); if (it != rd.end()) { ed = std::max(ed, std::fabs(b.coef[0] - it->second[0])); ++nd; } }
        for (const auto& b : tp.uvw[i].blocks) { long k = mem_key(b.mjd, b.utc_start);
            auto it = ru.find(k); if (it != ru.end() && it->second.size() == 3) {
                eu = std::max(eu, std::fabs(b.u[0] - it->second[0]));
                eu = std::max(eu, std::fabs(b.v[0] - it->second[1]));
                eu = std::max(eu, std::fabs(b.w[0] - it->second[2])); ++nu; } }
        printf("    %-9s задержка %.2f м (%d блоков), uvw %.2f м (%d блоков)\n",
               tel.c_str(), ed * C, nd, eu * C, nu);
        if (nd == 0 || nu == 0) ok = false;                 // не сопоставилось -> провал
        if (ed > 2e-8 || eu > 2e-8) ok = false;             // > ~6 м -> провал
        worst_d = std::max(worst_d, ed); worst_u = std::max(worst_u, eu);
    }
    printf("  ХУДШЕЕ: задержка %.2f м, uvw %.2f м\n", worst_d * C, worst_u * C);
    printf("---------------------------------------------------------------------\n");
    printf("  РЕЗУЛЬТАТ: %s (адаптер в памяти == сверенный файловый путь)\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
