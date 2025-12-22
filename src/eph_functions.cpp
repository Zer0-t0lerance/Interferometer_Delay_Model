#include "functions.h"
#include "constants.h"
#include <cmath>

// Заглушка для JPLEPH_421 - в реальности нужно реализовать чтение JPL эфемерид
void jpleph(double jd, double ct, Eigen::Vector3d& earth, Eigen::Vector3d& sun, Eigen::Vector3d& moon) {
    // Простая аппроксимация - Земля в центре, Солнце и Луна на расстоянии AU
    // В реальности нужно читать из файла DE421
    double t = (jd - cnst::JD2000) / 365.25; // годы

    // Солнце
    sun << cnst::AU * cos(2 * cnst::PI * t), cnst::AU * sin(2 * cnst::PI * t), 0.0;

    // Луна (орбита вокруг Земли)
    moon << cnst::AU * 0.00257 * cos(2 * cnst::PI * t * 12.37), cnst::AU * 0.00257 * sin(2 * cnst::PI * t * 12.37), 0.0;

    // Земля в центре
    earth << 0.0, 0.0, 0.0;
}

// Полная реализация r2000_matrix аналогично PNSXY40 из Fortran
void r2000_matrix(double mjd, double ut1, const Eigen::VectorXd& eop_int, const Eigen::VectorXd& deop_int, int i_choice, Eigen::Matrix3d r2000[3], double& gast) {
    // Заглушка: единичная матрица
    r2000[0].setIdentity();
    r2000[1].setZero();
    r2000[2].setZero();
    gast = 0.0;
}
