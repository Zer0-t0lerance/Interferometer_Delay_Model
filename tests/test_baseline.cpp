// test_baseline.cpp
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include "../src/functions.h"

int main() {
    // Настройка формата вывода
    std::cout << std::scientific << std::setprecision(16);
    // Формат для Eigen
    const Eigen::IOFormat fmt(16, 0, " ", "\n", "", "", "", "");

    // Initialize r2000 as 3x3 matrix
    Eigen::Matrix3d r2000;
    double r2000_init[3][3] = {
        {0.9988360617015190, 0.4820309780423414E-01, 0.1727196188789892E-02},
        {-0.4820310447537553E-01, 0.9988375540050903, -0.3778973045624359E-04},
        {-0.1727009998570988E-02, -0.4551047279605037E-04, 0.9999985076815172}
    };
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            r2000(i, j) = r2000_init[i][j];

    // Вывод r2000
    std::cout << "r2000 (transformation matrix):\n" << r2000.format(fmt) << "\n\n";

    // Initialize xsta_j2000t
    std::vector<Eigen::Vector3d> xsta_j2000t(2);
    double xsta_init[2][3] = {
        {0.4788183147308413E+07, -0.4190717195172735E+07, -0.4369012217094282E+06},
        {0.5203403045778226E+07, 0.2420028101198892E+07, -0.2777592466280870E+07}
    };
    for (int i = 0; i < 2; ++i)
        xsta_j2000t[i] << xsta_init[i][0], xsta_init[i][1], xsta_init[i][2];

    // Вывод xsta_j2000t
    std::cout << "xsta_j2000t[0] (station 1 position, m):\n" << xsta_j2000t[0].format(fmt) << "\n\n";
    std::cout << "xsta_j2000t[1] (station 2 position, m):\n" << xsta_j2000t[1].format(fmt) << "\n\n";

    // Initialize vsta_j2000t
    std::vector<Eigen::Vector3d> vsta_j2000t(2);
    double vsta_init[2][3] = {
        {0.3055926636587654E+03, 0.3492143463604144E+03, -0.5142780989077815E+00},
        {-0.1764630853380511E+03, 0.3797874391789990E+03, 0.3196033011667836E+00}
    };
    for (int i = 0; i < 2; ++i)
        vsta_j2000t[i] << vsta_init[i][0], vsta_init[i][1], vsta_init[i][2];

    // Вывод vsta_j2000t
    std::cout << "vsta_j2000t[0] (station 1 velocity, m/s):\n" << vsta_j2000t[0].format(fmt) << "\n\n";
    std::cout << "vsta_j2000t[1] (station 2 velocity, m/s):\n" << vsta_j2000t[1].format(fmt) << "\n\n";

    // Initialize output vectors
    std::vector<Eigen::Vector3d> base_line(2);
    Eigen::Vector3d b_cfs;

    // Call baseline
    ariadna::baseline(r2000, xsta_j2000t, vsta_j2000t, base_line, b_cfs);

    // Print results
    std::cout << "base_line[0] (position, m):\n" << base_line[0].format(fmt) << "\n\n";
    std::cout << "base_line[1] (velocity, m/s):\n" << base_line[1].format(fmt) << "\n\n";
    std::cout << "b_cfs (crust-fixed baseline, m):\n" << b_cfs.format(fmt) << std::endl;

    return 0;
}