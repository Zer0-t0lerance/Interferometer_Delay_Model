// test_interp_iers.cpp
//
// Сборка (из корня репозитория):
//   g++ -std=c++17 -I.\external\eigen -I.\external ^
//       src\interp_iers.cpp tests\test_interp_iers.cpp -o tests\test_interp_iers.exe
//
// Входные данные взяты из старого бенчмарка:
//   D:\prog\_Arriadna\INTERP_IERS\INTERP_TEST\interp_test\interp_test\interp_test.cpp
// 7 узлов EOP (MJD 57421..57427), точка интерполяции mjd=57422, utc=0.7711805555555555.
// Выход: eop_int[0]=UT1-UTC [сек], eop_int[1..2]=x,y [угл.сек].

#include "../src/functions.h"
#include <cstdio>

using namespace ariadna;

int main() {
    // MJD,   UT1-UTC,    UT1-TAI,     x,        y,        dpsi,      deps
    const double nodes[7][7] = {
        {57421.0, 0.0237813, -35.9762187, -0.005710, 0.302684, -0.00010, -0.00008},
        {57422.0, 0.0225960, -35.9774040, -0.006449, 0.304508, -0.00008, -0.00008},
        {57423.0, 0.0213968, -35.9786032, -0.007507, 0.306500, -0.00006, -0.00007},
        {57424.0, 0.0201401, -35.9798599, -0.008681, 0.308137, -0.00007, -0.00007},
        {57425.0, 0.0187401, -35.9812599, -0.009586, 0.309827, -0.00009, -0.00008},
        {57426.0, 0.0171296, -35.9828704, -0.009637, 0.311463, -0.00011, -0.00009},
        {57427.0, 0.0153150, -35.9846850, -0.009694, 0.313165, -0.00013, -0.00010},
    };

    std::vector<EOPData> eop(7);
    for (int i = 0; i < 7; ++i) {
        eop[i].mjd     = nodes[i][0];
        eop[i].ut1_utc = nodes[i][1];
        eop[i].ut1_tai = nodes[i][2];
        eop[i].x       = nodes[i][3];
        eop[i].y       = nodes[i][4];
        eop[i].dpsi    = nodes[i][5];
        eop[i].deps    = nodes[i][6];
    }

    Observation obs{};
    obs.mjd = 57422;
    obs.utc = 0.7711805555555555;

    Eigen::VectorXd eop_int(5);
    interp_iers(obs, eop, eop_int);

    printf("=====================================================================\n");
    printf("  INTERP_IERS benchmark (mjd=%d, utc=%.16f)\n", obs.mjd, obs.utc);
    printf("---------------------------------------------------------------------\n");
    printf("  eop_int[0]  UT1-UTC = % .15f  [sec]\n",    eop_int(0));
    printf("  eop_int[1]  x       = % .15f  [arcsec]\n", eop_int(1));
    printf("  eop_int[2]  y       = % .15f  [arcsec]\n", eop_int(2));
    printf("  eop_int[3]  dpsi    = % .15f  [arcsec]\n", eop_int(3));
    printf("  eop_int[4]  deps    = % .15f  [arcsec]\n", eop_int(4));
    printf("=====================================================================\n");
    printf("  Сверить с прогоном старой ARIADNA INTERP_IERS на этих же узлах.\n");
    return 0;
}
