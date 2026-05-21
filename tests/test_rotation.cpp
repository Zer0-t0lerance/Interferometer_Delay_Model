#include <iostream>
#include <cstdio>
#include "../src/functions.h"

using namespace ariadna;

int main() {
    printf("VERIFICATION: Earth Rotation Matrices (ITRF -> J2000, IAU 1980)\n");
    printf("---------------------------------------------------------------------------\n");

    // Исходная эпоха теста (примерно 1 января 2018 года)
    double jd_tdb = 2458120.500000;
    double jd_ut1 = 2458120.501234; // С учетом поправки UT1-UTC
    
    // Смещение полюса в радианах (типичные значения порядка сотых долей секунд дуги)
    double xp = 0.1234 * (cnst::PI / (180.0 * 3600.0)); 
    double yp = 0.4567 * (cnst::PI / (180.0 * 3600.0));

    Eigen::Matrix3d R2000, dR2000_dt;
    get_r2000_matrices(jd_tdb, jd_ut1, xp, yp, R2000, dR2000_dt);

    // Вывод матрицы координат
    printf("R2000 Matrix (Coordinates rotation):\n");
    for (int i = 0; i < 3; ++i) {
        printf("  [% .15le  % .15le  % .15le]\n", R2000(i,0), R2000(i,1), R2000(i,2));
    }

    // Вывод матрицы скоростей
    printf("\ndR2000_dt Matrix (Velocities derivative):\n");
    for (int i = 0; i < 3; ++i) {
        printf("  [% .15le  % .15le  % .15le]\n", dR2000_dt(i,0), dR2000_dt(i,1), dR2000_dt(i,2));
    }

    // Проверка ортогональности матрицы поворота (Det должна быть равна 1.0)
    double det = R2000.determinant();
    printf("\nMatrix Orthogonality Check:\n");
    printf("  Determinant of R2000: %.15f %s\n", det, (std::abs(det - 1.0) < 1e-12 ? "[OK]" : "[FAIL]"));

    return 0;
}