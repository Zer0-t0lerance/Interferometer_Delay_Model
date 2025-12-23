#include <iostream>
#include <iomanip>
#include "..\\src\\functions.h"

using namespace std;

void check(string name, const Eigen::Vector3d& act, double x, double y, double z) {
    Eigen::Vector3d exp(x, y, z);
    double diff = (act - exp).norm();
    cout << left << setw(25) << name 
         << " | Diff: " << scientific << setprecision(2) << diff;
    if (diff < 1e-8) cout << " | [ OK ]" << endl;
    else cout << " | [ FAIL ]" << endl;
}

int main() {
    // 1. Исходные данные (Debug output 1)
    Eigen::MatrixXd xsta(3, 2);
    xsta.col(0) << 0.4788183147308413E+07, -0.4190717195172735E+07, -0.4369012217094282E+06;
    xsta.col(1) << 0.5203403045778226E+07,  0.2420028101198892E+07, -0.2777592466280870E+07;

    Eigen::MatrixXd vsta(3, 2);
    vsta.col(0) <<  0.3055926636587654E+03,  0.3492143463604144E+03, -0.5142780989077815E+00;
    vsta.col(1) << -0.1764630853380511E+03,  0.3797874391789990E+03,  0.3196033011667836E+00;

    Eigen::MatrixXd asta = Eigen::MatrixXd::Zero(3, 2);

    // Матрица r2000 из лога (Debug output 2)
    Eigen::Matrix3d r2000;
    r2000 << 0.9988360617015190E+00,  0.4820309780423414E-01,  0.1727196188789892E-02,
            -0.4820310447537553E-01,  0.9988375540050903E+00, -0.3778973045624359E-04,
            -0.1727009998570988E-02, -0.4551047279605037E-04,  0.9999985076815172E+00;

    // 2. Вызов функции
    Eigen::MatrixXd base_line(3, 2);
    Eigen::Vector3d b_cfs;
    ariadna::baseline(r2000, xsta, vsta, asta, base_line, b_cfs);

    // 3. Валидация (сравнение с Debug output 2)
    cout << "================= BASELINE VALIDATION =================" << endl;
    
    // Проверка позиции базы (base_line col 1 в фортране)
    check("Base Position (J2000)", base_line.col(0), 
          0.4152198984698132E+06, 0.6610745296371628E+07, -0.2340691244571441E+07);

    // Проверка скорости базы (base_line col 2 в фортране)
    check("Base Velocity (J2000)", base_line.col(1), 
          -0.4820557489968165E+03, 0.3057309281858454E+02, 0.8338814000745651E+00);

    // Проверка b_cfs (Crust-fixed)
    // b_cfs = rT2000 * base_pos_j2000
    // Считаем вручную эталон b_cfs, так как в твоем дропе самого значения b_cfs нет, 
    // но есть матрицы для его получения.
    Eigen::Matrix3d rT2000;
    rT2000 << 0.9988360617015190E+00, -0.4820310447537553E-01, -0.1727009998570988E-02,
              0.4820309780423414E-01,  0.9988375540050903E+00, -0.4551047279605037E-04,
              0.1727196188789892E-02, -0.3778973045624359E-04,  0.9999985076815172E+00;
    
    Eigen::Vector3d base_pos_j2000(0.4152198984698132E+06, 0.6610745296371628E+07, -0.2340691244571441E+07);
    Eigen::Vector3d b_cfs_ref = rT2000 * base_pos_j2000;

    check("Crust-Fixed (b_cfs)", b_cfs, b_cfs_ref(0), b_cfs_ref(1), b_cfs_ref(2));

    return 0;
}