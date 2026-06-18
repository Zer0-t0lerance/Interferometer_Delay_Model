#include <iostream>
#include <vector>
#include <string>
#include "../src/catalog_bridge.h"

using namespace std;
using namespace ariadna;

int main() {
    cout << "VERIFICATION: Ocean Tide Catalog Parsing and Mapping" << endl;
    cout << "--------------------------------------------------------" << endl;

    // 1. Читаем сырой каталог
    vector<oc_record> raw_oc_data;
    char filename[] = "external/catalogs/VLBI_ocload_40.cat";

    int res = ReadOC(raw_oc_data, filename);
    if (res != 1) {
        cerr << "ERROR: Cannot open or parse file: " << filename << endl;
        cerr << "Make sure the file exists in 'external/catalogs/'" << endl;
        return -1;
    }

    cout << "[OK] Successfully read " << raw_oc_data.size() << " records from catalog." << endl;
    if (raw_oc_data.empty()) return 0;

    // 2. Создаем вектор станций и берем первые две станции из каталога для теста
    vector<Station> stations;
    for (size_t i = 0; i < min(size_t(2), raw_oc_data.size()); ++i) {
        Station st;
        st.name = string(raw_oc_data[i].telescope);
        stations.push_back(st);
    }

    cout << "[OK] Created test stations: " << stations[0].name << ", " << stations[1].name << endl;

    // 3. Вызываем нашу функцию-прокладку!
    map_ocean_tides_to_stations(raw_oc_data, stations);

    // 4. Проверяем результат (выводим матрицы Eigen для первой станции)
    cout << "\n========================================================" << endl;
    cout << "Mapped Data for Station: " << stations[0].name << endl;
    cout << "========================================================" << endl;
    
    cout << "--- AMPLITUDES (3 axes x 11 waves) ---" << endl;
    cout << stations[0].tide_data.amplitudes << endl;
    
    cout << "\n--- PHASES (3 axes x 11 waves) ---" << endl;
    cout << stations[0].tide_data.phases << endl;
    cout << "========================================================" << endl;

    return 0;
}