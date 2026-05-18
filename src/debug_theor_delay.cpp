#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>

namespace cnst {
    constexpr double C = 2.99792458e8;
    constexpr double C2 = C * C;
    constexpr double C3 = C2 * C;
    constexpr double gamma_PPN = 1.0;
    
    // Выведено из твоего дампа (U_Sun * C2 и U_Earth * C2)
    constexpr double GSUN    = 1.327124400409440e20; 
    constexpr double GEARTH  = 3.98600436233e14; 
    constexpr double GMOON   = 4.902801e12; // Стандарт IERS, влияет на 16-й знак
}

void check_val(const std::string& name, double calc, double ref, double tol = 1e-12) {
    double diff = std::abs(calc - ref);
    std::cout << std::left << std::setw(25) << name 
              << " | Calc: " << std::scientific << std::setprecision(15) << std::setw(24) << calc 
              << " | Ref: "  << std::setw(24) << ref 
              << " | Diff: " << std::setw(12) << diff 
              << (diff <= tol ? " [OK]" : " [FAIL]") << std::endl;
}

int main() {
    std::cout << std::string(100, '=') << "\n";
    std::cout << "STRICT FORTRAN TRANSLATION: THEOR_DELAY\n";
    std::cout << std::string(100, '-') << "\n";

    double C = cnst::C;
    double C2 = cnst::C2;
    double C3 = cnst::C3;

    // 1. Исходные векторы из дампа (18JAN02XA_N004.ngs)
    Eigen::Vector3d K_s(0.9340717157118170, 0.08020509321046379, 0.3479614532247550);
    
    Eigen::Matrix<double, 3, 2> base_line;
    base_line.col(0) << 415219.89847, 6610745.29637, -2340691.24457;
    
    Eigen::Vector3d X_1(-30326382850.88918, 132854720105.26163, 57575097849.76741);
    Eigen::Vector3d X_2(-30325967630.99071, 132861330850.55800, 57572757158.52284);
    Eigen::Vector3d dotX_1(-29313.55154, -5425.75602, -2504.78267);
    Eigen::Vector3d dotX_2(-29795.60728676801, -5395.182928854237, -2503.948787164999);

    Eigen::Matrix3d Earth;
    Earth.col(0) << -30331171034.03649, 132858910822.4568, 57575534750.98912;
    Earth.col(1) << -29619.14420142996, -5774.970368033236, -2504.268390466165;
    
    Eigen::Vector3d X_Sun(268184443.84710, 849738216.96790, 348845058.27214);
    Eigen::Vector3d V_Sun(-10.14974, 7.86276, 3.67712);

    Eigen::Vector3d X_Moon(-30457209591.74120, 133170659155.41185, 57696076560.55701);
    Eigen::Vector3d V_Moon(-30654.27062, -6147.62886, -2567.30178);

    // Восстанавливаем геоцентрические векторы станций
    Eigen::Vector3d xsta_1 = X_1 - Earth.col(0);
    Eigen::Vector3d xsta_2 = X_2 - Earth.col(0);
    Eigen::Vector3d vsta_1 = dotX_1 - Earth.col(1);
    Eigen::Vector3d vsta_2 = dotX_2 - Earth.col(1);

    // 2. Скалярные произведения
    double K_starB = K_s.dot(base_line.col(0));
    
    // Фортран: R_geocen(j) = -Sun(j,1)
    Eigen::Vector3d Sun_geo = X_Sun - Earth.col(0);
    Eigen::Vector3d R_geocen = -Sun_geo;
    double R_geo = R_geocen.norm();

    double U_Sun = cnst::GSUN / C2;
    double U_Earth = cnst::GEARTH / C2;
    double U = U_Sun / R_geo; // Фортран v.1.0 U = U_Sun/R_geo !!!!

    // =========================================================================
    // ДОСЛОВНЫЙ ПЕРЕВОД БЛОКА ШАПИРО (term1) ИЗ ФОРТРАНА
    // =========================================================================
    
    // -- Солнце --
    Eigen::Vector3d X1_Sun = X_Sun - X_1;
    double t1_Sun = K_s.dot(X1_Sun) / C;
    if (t1_Sun < 0.0) t1_Sun = 0.0;

    Eigen::Vector3d X1_Sun_t1 = X_Sun - V_Sun * t1_Sun;
    Eigen::Vector3d R1_Sun_t1 = X_1 - X1_Sun_t1;
    Eigen::Vector3d R2_Sun_t1 = X_2 - Earth.col(1) * (K_starB / C) - X1_Sun_t1;

    double C_Sun = (1.0 + cnst::gamma_PPN) * cnst::GSUN / C3;
    double w1 = R1_Sun_t1.norm();
    double numer1 = w1 + K_s.dot(R1_Sun_t1);
    double w2 = R2_Sun_t1.norm();
    double denom1 = w2 + K_s.dot(R2_Sun_t1);
    
    double delta_t_grav_Sun = C_Sun * std::log(numer1 / denom1);

    // Добавочный член от Солнца (Eq. 14)
    double C_Sun1 = C_Sun * C_Sun * C;
    Eigen::Vector3d N_hat = R1_Sun_t1 / w1;
    Eigen::Vector3d NplusK = N_hat + K_s;
    double V1 = base_line.col(0).dot(NplusK);
    double V2 = numer1 * numer1;
    double add_grav_Sun1 = C_Sun1 * V1 / V2;

    // -- Луна --
    Eigen::Vector3d X1_Moon = X_Moon - X_1;
    double t1_Moon = K_s.dot(X1_Moon) / C;
    if (t1_Moon < 0.0) t1_Moon = 0.0;

    Eigen::Vector3d X1_Moon_t1 = X_Moon - V_Moon * t1_Moon;
    Eigen::Vector3d R1_Moon_t1 = X_1 - X1_Moon_t1;
    Eigen::Vector3d R2_Moon_t1 = X_2 - Earth.col(1) * (K_starB / C) - X1_Moon_t1;

    double C_Moon = (1.0 + cnst::gamma_PPN) * cnst::GMOON / C3;
    double w1_m = R1_Moon_t1.norm();
    double numer2 = w1_m + K_s.dot(R1_Moon_t1);
    double w2_m = R2_Moon_t1.norm();
    double denom2 = w2_m + K_s.dot(R2_Moon_t1);
    
    double delta_t_grav_Moon = C_Moon * std::log(numer2 / denom2);

    // -- Земля --
    double C_Earth_grav = (1.0 + cnst::gamma_PPN) * cnst::GEARTH / C3;
    double w1_e = xsta_1.norm();
    double numer4 = w1_e + K_s.dot(xsta_1);
    double w2_e = xsta_2.norm();
    double denom4 = w2_e + K_s.dot(xsta_2);
    
    double delta_t_grav_Earth = C_Earth_grav * std::log(numer4 / denom4);

    // -- Итог гравитации --
    double delta_t_grav_Pl = 0.0; // Планеты в дампе отключены
    double delta_t_grav = delta_t_grav_Sun + delta_t_grav_Moon + delta_t_grav_Pl + delta_t_grav_Earth;
    
    double term1 = delta_t_grav + add_grav_Sun1;

    // =========================================================================
    // ОСТАЛЬНАЯ ВАКУУМНАЯ ЗАДЕРЖКА
    // =========================================================================
    double term2a = K_starB / C;
    double term2b = 1.0 - ((1.0 + cnst::gamma_PPN) * U);
    double term2c = Earth.col(1).squaredNorm() / (2.0 * C2);
    double term2d = Earth.col(1).dot(vsta_2) / C2;

    double term2bcd = term2b - term2c - term2d;
    double term2 = term2a * term2bcd;

    double VdotB = Earth.col(1).dot(base_line.col(0));
    double KdotV = K_s.dot(Earth.col(1));

    double term3a = VdotB / C2;
    double term3b = 1.0 + KdotV / (2.0 * C);
    double term3 = term3a * term3b;

    double numer_Eq9 = term1 - term2 - term3;
    
    Eigen::Vector3d vec_sum = Earth.col(1) + vsta_2;
    double den_Eq9 = 1.0 + K_s.dot(vec_sum) / C;
    
    double tv2_tv1 = numer_Eq9 / den_Eq9;

    // =========================================================================
    // ТРОПОСФЕРА И АППАРАТУРА (Debug 7, 8, 9)
    // =========================================================================
    Eigen::Vector3d w2_w1 = vsta_2 - vsta_1;
    double K_dotw2w1 = K_s.dot(w2_w1);

    double Datmc_d1 = -0.1225039869341168e-07;
    double Datmc_w1 = -0.1348080932573043e-08;
    double Datmc11 = Datmc_d1 + Datmc_w1; // Datmc(1,1)
    
    double tg2_tg1 = tv2_tv1 + Datmc11 * K_dotw2w1 / C;

    double Datmc_d2 =  0.1023313654522243e-07;
    double Datmc_w2 =  0.7304528775933354e-09;
    double Datmc21 = Datmc_d2 + Datmc_w2; // Datmc(2,1)

    double t2_t1_a = tg2_tg1 + Datmc11 + Datmc21;

    double Dtau_off1 = 0.1658888657089142e-10;
    double Dtau_off2 = -0.3833381489451837e-08;
    double Dt_temp1 =  0.8286346437956392e-10;
    double Dt_temp2 =  0.1445568241005169e-09;

    double t2_t1 = t2_t1_a + Dtau_off1 + Dtau_off2 + Dt_temp1 - Dt_temp2;

    // =========================================================================
    // ВЫВОД РЕЗУЛЬТАТОВ И СВЕРКА
    // =========================================================================
    std::cout << "\n--- REPLICATING FORTRAN DEBUG (6) ---\n";
    check_val("U_Sun", U_Sun, 1476.625038505688);
    check_val("U_Earth", U_Earth, 0.004435027977176442);
    check_val("R_geo", R_geo, 147097369863.7544, 1e-4); 
    check_val("U", U, 1.003841904089365e-08);
    check_val("term1", term1, -3.334437597304269e-10);
    check_val("term2a", term2a, 0.3455399756453089e-03);
    check_val("term2b", term2b, 0.9999999799231619);
    check_val("term2c", term2c, 0.5101029556441379e-08);
    check_val("term2d", term2d, 0.3374249253578403e-10);
    check_val("term2bcd", term2bcd, 0.9999999747883899);
    check_val("term2", term2, 0.3455399669336898e-03);
    check_val("term3", term3, -4.963692039367605e-07);
    check_val("numer_Eq9", numer_Eq9, -0.3450439311735128e-03);
    check_val("den_Eq9", den_Eq9, 0.9999028153242351);
    check_val("tv2_tv1", tv2_tv1, -0.3450774674152973e-03);

    std::cout << "\n--- REPLICATING FORTRAN DEBUG (8 & 9) ---\n";
    check_val("tg2_tg1", tg2_tg1, -0.3450774673949974e-03);
    check_val("t2_t1_a", t2_t1_a, -0.3450801022852005e-03);
    check_val("t2_t1 (FINAL)", t2_t1, -0.3450839807711632e-03);

    std::cout << std::string(100, '=') << "\n";
    return 0;
}