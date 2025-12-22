// test_aber_source.cpp
#include "../src/functions.h"
#include <iostream>

int main() {
    // Initialize observation
    ariadna::Observation obs;
    obs.sta1 = 0;
    obs.sta2 = 1;
    obs.mjd = 2451545; // J2000 epoch for simplicity
    obs.utc = 0.5;     // Half day
    obs.sou = 0;

    // Initialize r2000 as 9x3 matrix
    Eigen::MatrixXd r2000(9, 3);
    double r2000_init[3][3][3] = {
        {{0.856407603342079, 0.516298660270934, -1.307799799452249E-003},
         {-0.516297971645348, 0.856408600544260, 8.446234191276613E-004},
         {1.556088935770033E-003, -4.812753432626145E-005, 0.999998788134748}},
        {{-3.764904133986157E-005, 6.245029944245220E-005, 6.168655101856225E-008},
         {-6.245022671509238E-005, -3.764909170451489E-005, 9.537579702956507E-008},
         {-7.698574166572183E-011, -5.754292251350311E-011, 1.170274036610182E-013}},
        {{-4.553942307973770E-009, -2.745415038296109E-009, 6.954889614806429E-012},
         {2.745411365867157E-009, -4.553947611300586E-009, -4.498109878696921E-012},
         {4.290405149448649E-015, -5.725553756373111E-015, -6.961050533090952E-018}}
    };
    for (int k = 0; k < 3; ++k) {
        Eigen::Matrix3d temp;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                temp(i, j) = r2000_init[k][i][j];
        r2000.block<3, 3>(k * 3, 0) = temp.transpose();
    }

    // Initialize station parameters
    double lat_gd1 = -0.250899124329325;
    double lat_gd2 = -0.506968279695139;
    double h_geod1 = 189.973700276110;
    double h_geod2 = 248.925031604711;

    // Initialize vw
    std::vector<Eigen::Matrix3d> vw(2);
    double vw_init[2][3][3] = {
        {{-0.650092018774007, 0.718150315518193, -0.248275031864777},
         {-0.741362690371940, -0.671104583000652, 0.000000000000000E+000},
         {-0.166618511729085, 0.184061845575450, 0.968689583175407}},
        {{-0.374234314756092, 0.790069731247788, -0.485529090194486},
         {-0.903741928876677, -0.428077710223580, 0.000000000000000E+000},
         {-0.207844181177394, 0.438792996498103, 0.874220511412833}}
    };
    for (int i = 0; i < 2; ++i) {
        Eigen::Matrix3d temp;
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                temp(j, k) = vw_init[i][j][k];
        vw[i] = temp.transpose();
    }

    // Initialize Earth (corrected to use only velocity column)
    Eigen::Matrix3d earth;
    double earth_t[3][3] = {
        {-58823615826.18116, -27740.31190, 0.00251},  // Position, velocity, acceleration
        {123736441534.45139, -11133.38704, -0.00515},
        {53614977263.53795, -4825.09396, -0.00223}
    };
    earth.transposeInPlace();
    // Use only velocity column (col 1) and set others to zero
    earth.setZero();
    earth.col(1) << earth_t[0][1], earth_t[1][1], earth_t[2][1];

    // Initialize vsta_j2000t
    std::vector<Eigen::Vector3d> vsta_j2000t(2);
    double vsta_init[2][3] = {
        {-129.97431, -431.49437, 0.18113},
        {-225.01847, -339.06466, 0.33365}
    };
    for (int i = 0; i < 2; ++i)
        vsta_j2000t[i] << vsta_init[i][0], vsta_init[i][1], vsta_init[i][2];

    // Initialize k_s
    Eigen::Vector3d k_s;
    k_s << 0.2488843299651737, 0.8832588823692126, -0.3973793364200963;

    // Initialize output matrices
    Eigen::MatrixXd e(2, 2);
    Eigen::MatrixXd az(2, 2);

    // Call aber_source
    ariadna::aber_source(obs, r2000, lat_gd1, lat_gd2, h_geod1, h_geod2, k_s, earth, vsta_j2000t, vw, 2451545.5, 0.5, e, az);

    // Print results
    std::cout << "Elevation (rad, rad/s):\n" << e << "\n\n";
    std::cout << "Azimuth (rad, rad/s):\n" << az << std::endl;

    return 0;
}