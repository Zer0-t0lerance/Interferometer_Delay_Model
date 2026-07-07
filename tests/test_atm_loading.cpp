// test_atm_loading.cpp
//
// Unit-тест для чтения каталога атмосферной нагрузки (ReadATM) и адаптера
// map_atm_loading_to_stations.
//
// Сборка (из корня репозитория):
//   g++ -std=c++17 -Isrc -Iexternal\eigen ^
//       src\READ_CAT.cpp src\catalog_bridge.cpp tests\test_atm_loading.cpp ^
//       -o tests\test_atm_loading.exe
// Запуск: из корня репозитория (путь к каталогу относительный).
//
// Эталон — значения из самого каталога VLBI_atmload4_12.cat (первые станции).

#include "functions.h"
#include "catalog_bridge.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

using namespace ariadna;

// Находит файл вне зависимости от рабочей папки (корень репо или tests/).
static std::string resolve_path(const char* argv_override, std::vector<std::string> candidates) {
    if (argv_override) candidates.insert(candidates.begin(), argv_override);
    for (const auto& c : candidates) {
        std::ifstream f(c);
        if (f.good()) return c;
    }
    return candidates.back();
}

static int g_fail = 0;

static void check(const char* what, double got, double exp) {
    double d = std::fabs(got - exp);
    bool ok = d < 1e-9;
    if (!ok) ++g_fail;
    printf("  %-28s got=% .4f exp=% .4f %s\n", what, got, exp, ok ? "OK" : "FAIL");
}

int main(int argc, char** argv) {
    printf("=====================================================================\n");
    printf("  UNIT: ReadATM + map_atm_loading_to_stations\n");
    printf("---------------------------------------------------------------------\n");

    std::string path = resolve_path(argc > 1 ? argv[1] : nullptr,
        {"external/catalogs/VLBI_atmload4_12.cat",
         "../external/catalogs/VLBI_atmload4_12.cat"});
    printf("  Каталог: %s\n", path.c_str());

    std::vector<atm_record> atm;
    char fn[256];
    std::snprintf(fn, sizeof(fn), "%s", path.c_str());
    int rc = ReadATM(atm, fn);
    printf("  ReadATM rc=%d, stations=%zu (ожид. rc=1, >=198)\n", rc, atm.size());
    if (rc != 1) { printf("  FAIL: каталог не прочитан\n"); return 1; }
    if (atm.size() < 198) { ++g_fail; printf("  FAIL: мало станций\n"); }

    // --- Проверка сырого чтения: первая станция AIRA ---
    printf("\n  [ReadATM] станция[0] name='%s' (ожид. AIRA)\n", atm[0].telescope);
    if (std::strcmp(atm[0].telescope, "AIRA") != 0) { ++g_fail; printf("  FAIL: имя не AIRA\n"); }
    // Vertical row: A1 B1 A2 B2 A3 B3 b0 b1
    const double aira_V[8] = {-3.918, -0.328, 0.282, 0.695, -0.588, 0.513, 0.055, -0.502};
    for (int k = 0; k < 8; ++k) {
        char lbl[32]; std::snprintf(lbl, sizeof(lbl), "AIRA V coeff[%d]", k);
        check(lbl, atm[0].coeff[k], aira_V[k]);
    }

    // --- Проверка адаптера: маппинг на станции по имени ---
    std::vector<Station> st(3);
    st[0].name = "ALGOPARK";
    st[1].name = "AIRA";
    st[2].name = "NO_SUCH_STATION";
    map_atm_loading_to_stations(atm, st);

    printf("\n  [adapter] AIRA has_data=%d (ожид.1), NO_SUCH has_data=%d (ожид.0)\n",
           (int)st[1].atm_load.has_data, (int)st[2].atm_load.has_data);
    if (!st[1].atm_load.has_data) ++g_fail;
    if (st[2].atm_load.has_data) ++g_fail;

    // AIRA через адаптер должен совпасть с сырыми V-коэффициентами
    for (int k = 0; k < 8; ++k) {
        char lbl[40]; std::snprintf(lbl, sizeof(lbl), "adapter AIRA coef(V,%d)", k);
        check(lbl, st[1].atm_load.coef(0, k), aira_V[k]);
    }
    // ALGOPARK Vertical row
    const double algo_V[8] = {-0.089, -0.834, -0.212, 0.017, 0.587, -0.012, -0.047, -0.430};
    for (int k = 0; k < 8; ++k) {
        char lbl[44]; std::snprintf(lbl, sizeof(lbl), "adapter ALGOPARK coef(V,%d)", k);
        check(lbl, st[0].atm_load.coef(0, k), algo_V[k]);
    }
    // Несуществующая станция — нули
    check("NO_SUCH coef(V,0)==0", st[2].atm_load.coef(0, 0), 0.0);

    printf("=====================================================================\n");
    printf("  РЕЗУЛЬТАТ: %s (провалов: %d)\n", g_fail == 0 ? "PASS" : "FAIL", g_fail);
    return g_fail == 0 ? 0 : 1;
}
