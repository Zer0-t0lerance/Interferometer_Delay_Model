#include "functions.h"
#include "constants.h"
#include <cmath>

namespace ariadna {

void SITE_TIDE_SOLID(const Eigen::Vector3d& xsta_itrf, double lat_gcen, double lon_gcen,
                     const Eigen::Matrix<double, 3, 2>& sun, const Eigen::Matrix<double, 3, 2>& moon,
                     const Eigen::VectorXd& f, const Eigen::VectorXd& fd,
                     const Eigen::Matrix3d& vw_i, const Eigen::Vector2d& gast,
                     const Eigen::MatrixXd& r2000,
                     Eigen::Vector3d& dxtide, Eigen::Vector3d& dvtide) {
    
    const double DAS2R = 4.84813681109535993589e-06; 

    Eigen::Matrix3d R0 = r2000.block<3, 3>(0, 0);
    Eigen::Matrix3d R1 = r2000.block<3, 3>(3, 0);

    Eigen::Vector3d xsta_j2000 = R0 * xsta_itrf;
    Eigen::Vector3d vsta_j2000 = R1 * xsta_itrf;

    double rsta = xsta_j2000.norm();
    Eigen::Vector3d xsta_n = xsta_j2000 / rsta;
    Eigen::Vector3d dxsta_hat = vsta_j2000 / rsta;

    Eigen::Vector3d xsun = sun.col(0), vsun = sun.col(1);
    Eigen::Vector3d xmon = moon.col(0), vmon = moon.col(1);

    double rsun = xsun.norm(), rmon = xmon.norm();
    Eigen::Vector3d xsun_n = xsun / rsun, xmon_n = xmon / rmon;

    double pl[4] = {std::cos(lat_gcen), std::sin(lat_gcen), std::cos(lon_gcen), std::sin(lon_gcen)};
    
    double sc_sun = xsta_n.dot(xsun_n);
    double sc_mon = xsta_n.dot(xmon_n);

    // ========================================================================
    // STEP 1: Базовые приливы (Degree 2, 3)
    // ========================================================================
    double w1 = 1.5 * pl[1] * pl[1] - 0.5;
    double h2 = cnst::H02 + cnst::H22 * w1;
    double l2 = cnst::L02 + cnst::L22 * w1;
    double h2_2 = h2 / 2.0;

    double p2_sun = 3.0 * (h2_2 - l2) * sc_sun * sc_sun - h2_2;
    double p2_mon = 3.0 * (h2_2 - l2) * sc_mon * sc_mon - h2_2;

    double p3_sun = 2.5 * (cnst::H3 - 3.0 * cnst::L3) * sc_sun * sc_sun * sc_sun + 1.5 * (cnst::L3 - cnst::H3) * sc_sun;
    double p3_mon = 2.5 * (cnst::H3 - 3.0 * cnst::L3) * sc_mon * sc_mon * sc_mon + 1.5 * (cnst::L3 - cnst::H3) * sc_mon;

    double x2_sun = 3.0 * l2 * sc_sun, x2_mon = 3.0 * l2 * sc_mon;
    double x3_sun = 1.5 * cnst::L3 * (5.0 * sc_sun * sc_sun - 1.0);
    double x3_mon = 1.5 * cnst::L3 * (5.0 * sc_mon * sc_mon - 1.0);

    double mass_ratio_sun = cnst::GSUN / cnst::GEARTH;

    double fgm[2][2]; 
    fgm[0][0] = cnst::AE * std::pow(cnst::AE / rsun, 3) * mass_ratio_sun;
    fgm[1][0] = (cnst::AE / rsun) * fgm[0][0];
    fgm[0][1] = cnst::AE * std::pow(cnst::AE / rmon, 3) * cnst::MU; 
    fgm[1][1] = (cnst::AE / rmon) * fgm[0][1];

    Eigen::Vector3d dx1 = fgm[0][0] * (x2_sun * xsun_n + p2_sun * xsta_n) +
                          fgm[0][1] * (x2_mon * xmon_n + p2_mon * xsta_n) +
                          fgm[1][0] * (x3_sun * xsun_n + p3_sun * xsta_n) +
                          fgm[1][1] * (x3_mon * xmon_n + p3_mon * xsta_n);

    double a1 = 6.0 * (h2_2 - l2) * sc_sun, b1 = 6.0 * (h2_2 - l2) * sc_mon, c1 = 3.0 * l2;
    double a2 = 7.5 * (cnst::H3 - 3.0 * cnst::L3) * sc_sun * sc_sun + 1.5 * (cnst::L3 - cnst::H3);
    double b2 = 7.5 * (cnst::H3 - 3.0 * cnst::L3) * sc_mon * sc_mon + 1.5 * (cnst::L3 - cnst::H3);
    double c2 = 15.0 * cnst::L3;

    double dr_sun = xsun.dot(vsun) / rsun;
    double dr_mon = xmon.dot(vmon) / rmon;

    Eigen::Vector3d drsun_hat = vsun / rsun - xsun * dr_sun / (rsun * rsun);
    Eigen::Vector3d drmon_hat = vmon / rmon - xmon * dr_mon / (rmon * rmon);

    double drsun_r_hat = drsun_hat.dot(xsta_n) + xsun_n.dot(dxsta_hat);
    double drmon_r_hat = drmon_hat.dot(xsta_n) + xmon_n.dot(dxsta_hat);

    Eigen::Vector3d dv1 = fgm[0][0] * (-3.0 * dr_sun / rsun) * (x2_sun * xsun_n + p2_sun * xsta_n) +
                          fgm[0][1] * (-3.0 * dr_mon / rmon) * (x2_mon * xmon_n + p2_mon * xsta_n) +
                          fgm[0][0] * (x2_sun * drsun_hat + p2_sun * dxsta_hat + (a1 * xsta_n + c1 * xsun_n) * drsun_r_hat) +
                          fgm[0][1] * (x2_mon * drmon_hat + p2_mon * dxsta_hat + (b1 * xsta_n + c1 * xmon_n) * drmon_r_hat) +
                          fgm[1][0] * (-4.0 * dr_sun / rsun) * (x3_sun * xsun_n + p3_sun * xsta_n) +
                          fgm[1][1] * (-4.0 * dr_mon / rmon) * (x3_mon * xmon_n + p3_mon * xsta_n) +
                          fgm[1][0] * (x3_sun * drsun_hat + p3_sun * dxsta_hat + (a2 * xsta_n + c2 * xsun_n) * drsun_r_hat) +
                          fgm[1][1] * (x3_mon * drmon_hat + p3_mon * dxsta_hat + (b2 * xsta_n + c2 * xmon_n) * drmon_r_hat);

    // ========================================================================
    // Crust-Fixed Векторы для дальнейших коррекций
    // ========================================================================
    Eigen::Vector3d xsun_cf = R0.transpose() * xsun;
    Eigen::Vector3d xmon_cf = R0.transpose() * xmon;
    Eigen::Vector3d vsun_cf = R0.transpose() * vsun + R1.transpose() * xsun;
    Eigen::Vector3d vmon_cf = R0.transpose() * vmon + R1.transpose() * xmon;

    double rsun_cf = xsun_cf.norm(), rmon_cf = xmon_cf.norm();
    Eigen::Vector3d xsun_ncf = xsun_cf / rsun_cf, xmon_ncf = xmon_cf / rmon_cf;
    double drsun_cf = xsun_cf.dot(vsun_cf) / rsun_cf;
    double drmon_cf = xmon_cf.dot(vmon_cf) / rmon_cf;

    Eigen::Vector3d vsun_ncf = vsun_cf / rsun_cf - xsun_ncf / (rsun_cf * rsun_cf * drsun_cf);
    Eigen::Vector3d vmon_ncf = vmon_cf / rmon_cf - xmon_ncf / (rmon_cf * rmon_cf * drmon_cf);

    double p2[5][2], dp2_m[5][2]; 
    
    // Sun
    p2[0][0] = 1.5 * xsun_ncf(2) * xsun_ncf(2) - 0.5;
    p2[1][0] = 3.0 * xsun_ncf(0) * xsun_ncf(2);
    p2[2][0] = 3.0 * xsun_ncf(1) * xsun_ncf(2);
    p2[3][0] = 3.0 * (xsun_ncf(0) * xsun_ncf(0) - xsun_ncf(1) * xsun_ncf(1));
    p2[4][0] = 6.0 * xsun_ncf(0) * xsun_ncf(1);
    
    dp2_m[0][0] = 3.0 * xsun_ncf(2) * vsun_ncf(2);
    dp2_m[1][0] = 3.0 * (vsun_ncf(0) * xsun_ncf(2) + xsun_ncf(0) * vsun_ncf(2));
    dp2_m[2][0] = 3.0 * (vsun_ncf(1) * xsun_ncf(2) + xsun_ncf(1) * vsun_ncf(2));
    dp2_m[3][0] = 6.0 * (vsun_ncf(0) * xsun_ncf(0) - xsun_ncf(1) * vsun_ncf(1));
    dp2_m[4][0] = 6.0 * (vsun_ncf(0) * xsun_ncf(1) + xsun_ncf(0) * vsun_ncf(1));

    // Moon
    p2[0][1] = 1.5 * xmon_ncf(2) * xmon_ncf(2) - 0.5;
    p2[1][1] = 3.0 * xmon_ncf(0) * xmon_ncf(2);
    p2[2][1] = 3.0 * xmon_ncf(1) * xmon_ncf(2);
    p2[3][1] = 3.0 * (xmon_ncf(0) * xmon_ncf(0) - xmon_ncf(1) * xmon_ncf(1));
    p2[4][1] = 6.0 * xmon_ncf(0) * xmon_ncf(1);

    dp2_m[0][1] = 3.0 * xmon_ncf(2) * vmon_ncf(2);
    dp2_m[1][1] = 3.0 * (vmon_ncf(0) * xmon_ncf(2) + xmon_ncf(0) * vmon_ncf(2));
    dp2_m[2][1] = 3.0 * (vmon_ncf(1) * xmon_ncf(2) + xmon_ncf(1) * vmon_ncf(2));
    dp2_m[3][1] = 6.0 * (vmon_ncf(0) * xmon_ncf(0) - xmon_ncf(1) * vmon_ncf(1));
    dp2_m[4][1] = 6.0 * (vmon_ncf(0) * xmon_ncf(1) + xmon_ncf(0) * vmon_ncf(1));

    double s2phi = 2.0 * pl[0] * pl[1];
    double c2phi = pl[0] * pl[0] - pl[1] * pl[1];
    double s2lam = 2.0 * pl[2] * pl[3];
    double c2lam = pl[2] * pl[2] - pl[3] * pl[3];
    double coef[2] = {-3.0 * drsun_cf / rsun_cf, -3.0 * drmon_cf / rmon_cf};

    double dr_s = 0, dn_s = 0, de_s = 0, dvr = 0, dvn = 0, dve = 0;

    auto rotate_to_j2000 = [&](double dr, double de, double dn, double dvr_in, double dve_in, double dvn_in, Eigen::Vector3d& dx_out, Eigen::Vector3d& dv_out) {
        Eigen::Vector3d dr_ven(dr, de, dn);
        Eigen::Vector3d dv_ven(dvr_in, dve_in, dvn_in);
        Eigen::Vector3d work = vw_i * dr_ven;
        dx_out = R0 * work;
        Eigen::Vector3d dwork = vw_i * dv_ven;
        dv_out = R0 * dwork + R1 * work;
    };

    // ========================================================================
    // CORR_H2L2
    // ========================================================================
    for (int j = 0; j < 2; j++) {
        double dr1_0 = -0.5 * cnst::HI_1 * fgm[0][j] * s2phi * (p2[1][j] * pl[3] - p2[2][j] * pl[2]);
        double dn1_0 = -cnst::LI_1 * fgm[0][j] * c2phi * (p2[1][j] * pl[3] - p2[2][j] * pl[2]);
        double de1_0 = -cnst::LI_1 * fgm[0][j] * pl[1] * (p2[1][j] * pl[2] + p2[2][j] * pl[3]);

        double dr1_1 = -0.5 * cnst::HI_1 * fgm[0][j] * s2phi * (coef[j] * (p2[1][j] * pl[3] - p2[2][j] * pl[2]) + dp2_m[1][j] * pl[3] - dp2_m[2][j] * pl[2]);
        double dn1_1 = -cnst::LI_1 * fgm[0][j] * c2phi * (coef[j] * (p2[1][j] * pl[3] - p2[2][j] * pl[2]) + dp2_m[1][j] * pl[3] - dp2_m[2][j] * pl[2]);
        double de1_1 = -cnst::LI_1 * fgm[0][j] * pl[1] * (coef[j] * (p2[1][j] * pl[2] + p2[2][j] * pl[3]) + dp2_m[1][j] * pl[2] + dp2_m[2][j] * pl[3]);

        double dr2_0 = -0.25 * cnst::HI_2 * fgm[0][j] * pl[0] * pl[0] * (p2[3][j] * s2lam - p2[4][j] * c2lam);
        double dn2_0 =  0.25 * cnst::LI_2 * fgm[0][j] * s2phi * (p2[3][j] * s2lam - p2[4][j] * c2lam);
        double de2_0 = -0.50 * cnst::LI_2 * fgm[0][j] * pl[0] * (p2[3][j] * c2lam + p2[4][j] * s2lam);

        double dr2_1 = -0.25 * cnst::HI_2 * fgm[0][j] * pl[1] * pl[1] * (coef[j] * (p2[4][j] * s2lam - p2[4][j] * c2lam) + dp2_m[4][j] * s2lam - dp2_m[4][j] * c2lam);
        double dn2_1 =  0.25 * cnst::LI_2 * fgm[0][j] * s2phi * (coef[j] * (p2[4][j] * s2lam - p2[4][j] * c2lam) + dp2_m[4][j] * s2lam - dp2_m[4][j] * c2lam);
        double de2_1 = -0.50 * cnst::LI_2 * fgm[0][j] * pl[1] * (coef[j] * (p2[4][j] * c2lam + p2[4][j] * s2lam) + dp2_m[4][j] * c2lam + dp2_m[4][j] * s2lam);

        dr_s += dr1_0 + dr2_0; dn_s += dn1_0 + dn2_0; de_s += de1_0 + de2_0;
        dvr  += dr1_1 + dr2_1; dvn  += dn1_1 + dn2_1; dve  += de1_1 + de2_1;
    }
    Eigen::Vector3d dxyz_h2l2, dvxyz_h2l2;
    rotate_to_j2000(dr_s, de_s, dn_s, dvr, dve, dvn, dxyz_h2l2, dvxyz_h2l2);
    Eigen::Vector3d dx2 = dx1 + dxyz_h2l2;
    Eigen::Vector3d dv2 = dv1 + dvxyz_h2l2;

    // ========================================================================
    // CORR_L1
    // ========================================================================
    dr_s = 0; dn_s = 0; de_s = 0; dvr = 0; dvn = 0; dve = 0;
    for (int j = 0; j < 2; j++) {
        double dn1_0 = -cnst::L1_1 * fgm[0][j] * pl[1] * pl[1] * (p2[1][j] * pl[2] + p2[2][j] * pl[3]);
        double de1_0 =  cnst::L1_1 * fgm[0][j] * pl[1] * c2phi * (p2[1][j] * pl[3] - p2[2][j] * pl[2]);
        double dn1_1 = -cnst::L1_1 * fgm[0][j] * pl[1] * pl[1] * (coef[j] * (p2[1][j] * pl[2] + p2[2][j] * pl[3]) + (dp2_m[1][j] * pl[2] + dp2_m[2][j] * pl[3]));
        double de1_1 =  cnst::L1_1 * fgm[0][j] * pl[1] * c2phi * (coef[j] * (p2[1][j] * pl[3] - p2[2][j] * pl[2]) + (dp2_m[1][j] * pl[3] - dp2_m[2][j] * pl[2]));

        double dn2_0 = -0.25 * cnst::L1_2 * fgm[0][j] * s2phi * (p2[3][j] * c2lam + p2[4][j] * s2lam);
        double de2_0 = -0.25 * cnst::L1_2 * fgm[0][j] * s2phi * pl[1] * (p2[3][j] * s2lam - p2[4][j] * c2lam);
        double dn2_1 = -0.25 * cnst::L1_2 * fgm[0][j] * s2phi * (coef[j] * (p2[3][j] * c2lam + p2[4][j] * s2lam) + (dp2_m[3][j] * c2lam + dp2_m[4][j] * s2lam));
        double de2_1 = -0.25 * cnst::L1_2 * fgm[0][j] * s2phi * pl[2] * (coef[j] * (p2[3][j] * s2lam - p2[4][j] * c2lam) + (dp2_m[3][j] * s2lam - dp2_m[4][j] * c2lam));

        dn_s += dn1_0 + dn2_0; de_s += de1_0 + de2_0;
        dvn  += dn1_1 + dn2_1; dve  += de1_1 + de2_1;
    }
    Eigen::Vector3d dxyz_l1, dvxyz_l1;
    rotate_to_j2000(dr_s, de_s, dn_s, dvr, dve, dvn, dxyz_l1, dvxyz_l1);
    Eigen::Vector3d dx3 = dx2 + dxyz_l1;
    Eigen::Vector3d dv3 = dv2 + dvxyz_l1;

    // ========================================================================
    // CORR_DIU
    // ========================================================================
    dr_s = 0; dn_s = 0; de_s = 0; dvr = 0; dvn = 0; dve = 0;
    double lam = std::atan2(pl[3], pl[2]);
    double cos2phi = pl[0]*pl[0] - pl[1]*pl[1];

    // АВТОМАТИЧЕСКОЕ ОПРЕДЕЛЕНИЕ РАЗМЕРА МАССИВА
    constexpr int num_diu = sizeof(cnst::AMP_DIU) / sizeof(cnst::AMP_DIU[0]);

    for (int i = 0; i < num_diu; i++) {
        double theta_f = 0.0, dtheta_f = 0.0;
        for (int j = 0; j < 5; j++) {
            theta_f += cnst::AMP_DIU[i][j] * f(j);
            dtheta_f += cnst::AMP_DIU[i][j] * fd(j);
        }
        theta_f = std::fmod(-theta_f * DAS2R + (gast(0) + cnst::PI), cnst::TWOPI);
        dtheta_f = -dtheta_f * DAS2R / 3.15576e9 + gast(1);

        double ctheta = std::cos(theta_f + lam), stheta = std::sin(theta_f + lam);
        double dctheta = -stheta * dtheta_f, dstheta = ctheta * dtheta_f;

        dr_s += (cnst::AMP_DIU[i][5] * stheta + cnst::AMP_DIU[i][6] * ctheta) * s2phi;
        dn_s += (cnst::AMP_DIU[i][7] * stheta + cnst::AMP_DIU[i][8] * ctheta) * cos2phi;
        de_s += (cnst::AMP_DIU[i][7] * ctheta - cnst::AMP_DIU[i][8] * stheta) * pl[1];

        dvr += (cnst::AMP_DIU[i][5] * dstheta + cnst::AMP_DIU[i][6] * dctheta) * s2phi;
        dvn += (cnst::AMP_DIU[i][7] * dstheta + cnst::AMP_DIU[i][8] * dctheta) * cos2phi;
        dve += (cnst::AMP_DIU[i][7] * dctheta - cnst::AMP_DIU[i][8] * dstheta) * pl[1];
    }
    Eigen::Vector3d dxyz_diu, dvxyz_diu;
    rotate_to_j2000(dr_s, de_s, dn_s, dvr, dve, dvn, dxyz_diu, dvxyz_diu);
    Eigen::Vector3d dx4 = dx3 + dxyz_diu * 1e-3; 
    Eigen::Vector3d dv4 = dv3 + dvxyz_diu * 1e-3;

    // ========================================================================
    // CORR_LON
    // ========================================================================
    dr_s = 0; dn_s = 0; de_s = 0; dvr = 0; dvn = 0; dve = 0;
    double coef_lon = (3.0 * pl[1] * pl[1] - 1.0) / 2.0;

    // АВТОМАТИЧЕСКОЕ ОПРЕДЕЛЕНИЕ РАЗМЕРА МАССИВА
    constexpr int num_lon = sizeof(cnst::AMP_LON) / sizeof(cnst::AMP_LON[0]);

    for (int i = 0; i < num_lon; i++) {
        double theta_f = 0.0, dtheta_f = 0.0;
        for (int j = 0; j < 5; j++) {
            theta_f += cnst::AMP_LON[i][j] * f(j);
            dtheta_f += cnst::AMP_LON[i][j] * fd(j);
        }
        theta_f = -theta_f * DAS2R;
        dtheta_f = -dtheta_f * DAS2R / 3.15576e9;

        double ctheta = std::cos(theta_f), stheta = std::sin(theta_f);
        double dctheta = -stheta * dtheta_f, dstheta = ctheta * dtheta_f;

        dr_s += (cnst::AMP_LON[i][5] * ctheta + cnst::AMP_LON[i][7] * stheta) * coef_lon;
        dn_s += (cnst::AMP_LON[i][7] * ctheta + cnst::AMP_LON[i][8] * stheta) * s2phi;

        dvr += (cnst::AMP_LON[i][5] * dctheta + cnst::AMP_LON[i][6] * dstheta) * coef_lon;
        dvn += (cnst::AMP_LON[i][7] * dctheta + cnst::AMP_LON[i][8] * dstheta) * s2phi;
    }
    Eigen::Vector3d dxyz_lon, dvxyz_lon;
    rotate_to_j2000(dr_s, de_s, dn_s, dvr, dve, dvn, dxyz_lon, dvxyz_lon);
    
    // Итоговый вектор
    dxtide = dx4 + dxyz_lon * 1e-3;
    dvtide = dv4 + dvxyz_lon * 1e-3;
}

} // namespace ariadna