// jpleph.cpp
//
// Порт интерфейса JPLEPH_421 (FORTRAN_SOURCE/JPLEPH_421.FOR.TXT).
// Соглашение оригинала:
//   EARTH(3,3) — барицентрические (SSB) позиция, скорость, ускорение [м, м/с, м/с^2];
//   SUN(3,2)   — ГЕОЦЕНТРИЧЕСКИЕ позиция и скорость Солнца [м, м/с] (Солнце_SSB - Земля_SSB);
//   MOON(3,2)  — ГЕОЦЕНТРИЧЕСКИЕ позиция и скорость Луны  [м, м/с] (Луна_SSB - Земля_SSB).
//
// Считаем всё через уже протестированную get_celestial_bodies() (SSB для всех тел),
// затем переводим Солнце и Луну в геоцентрическую систему.

#include "functions.h"

namespace ariadna {

// Матричный формат: earth(3x3) SSB pos/vel/acc; sun,moon(3x2) геоцентрические pos/vel.
void jpl_eph(double jd, double ct, Eigen::Matrix3d& earth,
             Eigen::MatrixXd& sun, Eigen::MatrixXd& moon) {
    Eigen::Matrix3d E, S, M;
    get_celestial_bodies(jd, ct, E, S, M);

    earth = E; // SSB: col0=pos, col1=vel, col2=acc

    sun.resize(3, 2);
    moon.resize(3, 2);
    // Геоцентрические позиция и скорость = тело_SSB - Земля_SSB
    sun.col(0)  = S.col(0) - E.col(0);
    sun.col(1)  = S.col(1) - E.col(1);
    moon.col(0) = M.col(0) - E.col(0);
    moon.col(1) = M.col(1) - E.col(1);
}

// Векторный формат: только позиции. earth — SSB; sun, moon — геоцентрические.
void jpleph(double jd, double ct, Eigen::Vector3d& earth,
            Eigen::Vector3d& sun, Eigen::Vector3d& moon) {
    Eigen::Matrix3d E, S, M;
    get_celestial_bodies(jd, ct, E, S, M);

    earth = E.col(0);            // SSB позиция Земли
    sun   = S.col(0) - E.col(0); // геоцентрическая позиция Солнца
    moon  = M.col(0) - E.col(0); // геоцентрическая позиция Луны
}

} // namespace ariadna
