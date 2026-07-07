// interp_iers.cpp
//
// Порт подпрограммы INTERP_IERS (Bizouard/Zharov) из FORTRAN_SOURCE/INTERP_IERS.for.TXT.
//
// Интерполирует UT1-UTC, координаты полюса (x, y) на момент наблюдения по
// 7 опорным узлам методом 4-точечного Лагранжа (LAGINT) и добавляет
// высокочастотные (суточные/полусуточные) поправки:
//   * PMUT1_OCEANS — океанические приливные члены (71 гармоника) для x, y, UT1, LOD;
//   * PM_GRAVI     — лунно-солнечная либрация (10 гармоник) для x, y.
//
// ВНИМАНИЕ по единицам: как и в оригинале, выход НЕ переводится в радианы —
//   eop_int(0) = UT1-UTC [сек],  eop_int(1..2) = x, y [угл. сек].
// (Ср. interp_eop(), который отдаёт радианы.) Значения dpsi/deps в исходном
// Фортране здесь не вычисляются; для полноты интерполируем их тем же LAGINT.
//
// Все табличные константы вынесены в constants.h (cnst::IERS_*).

#include "functions.h"

namespace ariadna {
namespace {

// 4-точечная интерполяция Лагранжа (порт SUBROUTINE LAGINT), окно из 4 узлов.
double lagint(const double* X, const double* Y, int n, double xint) {
    int K = 0;
    for (int I = 1; I <= n - 1; ++I) {
        if (xint >= X[I - 1] && xint < X[I]) K = I;
    }
    if (K < 2) K = 2;
    if (K > n - 2) K = n - 2;
    double yout = 0.0;
    for (int M = K - 1; M <= K + 2; ++M) {
        double term = Y[M - 1];
        for (int J = K - 1; J <= K + 2; ++J) {
            if (M != J) term *= (xint - X[J - 1]) / (X[M - 1] - X[J - 1]);
        }
        yout += term;
    }
    return yout;
}

// Фундаментальные аргументы chi=GMST+pi, l, l', F, D, Omega (Simon и др.) [рад]
// и их производные [рад/сут]. arg6_lin — линейный коэффициент ARG(6):
//   cnst::IERS_ARG6_LIN_OCEANS в PMUT1_OCEANS, cnst::IERS_ARG6_LIN_GRAVI в PM_GRAVI.
void compute_args(double rjd, double arg6_lin, double ARG[6], double DARG[6]) {
    const double SR = cnst::CARCRAD; // arcsec -> rad
    double T  = (rjd - 51544.5) / 36525.0;
    double T2 = T * T, T3 = T2 * T, T4 = T3 * T;

    ARG[0] = (67310.54841
              + (876600.0 * 3600.0 + 8640184.812866) * T
              + 0.093104 * T2 - 6.2e-6 * T3) * 15.0 + 648000.0;
    ARG[0] = std::fmod(ARG[0], 1296000.0) * SR;
    DARG[0] = (876600.0 * 3600.0 + 8640184.812866
               + 2.0 * 0.093104 * T - 3.0 * 6.2e-6 * T2) * 15.0;
    DARG[0] = DARG[0] * SR / 36525.0;

    ARG[1] = -0.00024470 * T4 + 0.051635 * T3 + 31.8792 * T2
             + 1717915923.2178 * T + 485868.249036;
    ARG[1] = std::fmod(ARG[1], 1296000.0) * SR;
    DARG[1] = -4.0 * 0.00024470 * T3 + 3.0 * 0.051635 * T2
              + 2.0 * 31.8792 * T + 1717915923.2178;
    DARG[1] = DARG[1] * SR / 36525.0;

    ARG[2] = -0.00001149 * T4 - 0.000136 * T3 - 0.5532 * T2
             + 129596581.0481 * T + 1287104.79305;
    ARG[2] = std::fmod(ARG[2], 1296000.0) * SR;
    DARG[2] = -4.0 * 0.00001149 * T3 - 3.0 * 0.000136 * T2
              - 2.0 * 0.5532 * T + 129596581.0481;
    DARG[2] = DARG[2] * SR / 36525.0;

    ARG[3] = 0.00000417 * T4 - 0.001037 * T3 - 12.7512 * T2
             + 1739527262.8478 * T + 335779.526232;
    ARG[3] = std::fmod(ARG[3], 1296000.0) * SR;
    DARG[3] = 4.0 * 0.00000417 * T3 - 3.0 * 0.001037 * T2
              - 2.0 * 12.7512 * T + 1739527262.8478;
    DARG[3] = DARG[3] * SR / 36525.0;

    ARG[4] = -0.00003169 * T4 + 0.006593 * T3 - 6.3706 * T2
             + 1602961601.2090 * T + 1072260.70369;
    ARG[4] = std::fmod(ARG[4], 1296000.0) * SR;
    DARG[4] = -4.0 * 0.00003169 * T3 + 3.0 * 0.006593 * T2
              - 2.0 * 6.3706 * T + 1602961601.2090;
    DARG[4] = DARG[4] * SR / 36525.0;

    ARG[5] = -0.00005939 * T4 + 0.007702 * T3 + 7.4722 * T2
             - arg6_lin * T + 450160.398036;
    ARG[5] = std::fmod(ARG[5], 1296000.0) * SR;
    DARG[5] = -4.0 * 0.00005939 * T3 + 3.0 * 0.007702 * T2
              + 2.0 * 7.4722 * T - arg6_lin;
    DARG[5] = DARG[5] * SR / 36525.0;
}

// Океанические суточные/полусуточные приливы -> поправки x, y ["], UT1, LOD [сек].
void pmut1_oceans(double rjd, double& cor_x, double& cor_y,
                  double& cor_ut1, double& cor_lod) {
    double ARG[6], DARG[6];
    compute_args(rjd, cnst::IERS_ARG6_LIN_OCEANS, ARG, DARG);

    cor_x = cor_y = cor_ut1 = cor_lod = 0.0;
    for (int j = 0; j < cnst::IERS_PMUT1_OCEANS_N; ++j) {
        const double* row = cnst::IERS_PMUT1_OCEANS[j];
        double ag = 0.0, dag = 0.0;
        for (int i = 0; i < 6; ++i) {
            ag  += row[i] * ARG[i];
            dag += row[i] * DARG[i];
        }
        ag = std::fmod(ag, cnst::TWOPI);
        double c = std::cos(ag), s = std::sin(ag);
        const double XSIN = row[6],  XCOS = row[7];
        const double YSIN = row[8],  YCOS = row[9];
        const double UTSIN = row[10], UTCOS = row[11];
        cor_x   += XCOS * c + XSIN * s;
        cor_y   += YCOS * c + YSIN * s;
        cor_ut1 += UTCOS * c + UTSIN * s;
        cor_lod -= (-UTCOS * s + UTSIN * c) * dag;
    }
    cor_x   *= 1.0e-6; // arcsec
    cor_y   *= 1.0e-6; // arcsec
    cor_ut1 *= 1.0e-6; // sec
    cor_lod *= 1.0e-6; // sec
}

// Лунно-солнечная либрация -> поправки x, y ["].
void pm_gravi(double rjd, double& cor_x, double& cor_y) {
    double ARG[6], DARG[6];
    compute_args(rjd, cnst::IERS_ARG6_LIN_GRAVI, ARG, DARG);

    cor_x = cor_y = 0.0;
    for (int j = 0; j < cnst::IERS_PM_GRAVI_N; ++j) {
        const double* row = cnst::IERS_PM_GRAVI[j];
        double ag = 0.0;
        for (int i = 0; i < 6; ++i) ag += row[i] * ARG[i];
        ag = std::fmod(ag, cnst::TWOPI);
        double c = std::cos(ag), s = std::sin(ag);
        const double XSIN = row[6], XCOS = row[7];
        const double YSIN = row[8], YCOS = row[9];
        cor_x += XCOS * c + XSIN * s;
        cor_y += YCOS * c + YSIN * s;
    }
    cor_x *= 1.0e-6; // arcsec
    cor_y *= 1.0e-6; // arcsec
}

} // anonymous namespace

void interp_iers(const Observation& obs, const std::vector<EOPData>& eop_data,
                 Eigen::VectorXd& eop_int) {
    const int N = cnst::EOP_NDATA; // 7 узлов
    if (eop_int.size() < 5) eop_int.resize(5);

    // Момент интерполяции в шкале UTC (mjd + доля суток)
    double rjd_int = static_cast<double>(obs.mjd) + obs.utc;

    std::vector<double> RJD(N), X(N), Y(N), UT1(N), DPSI(N), DEPS(N);
    for (int i = 0; i < N; ++i) {
        RJD[i]  = eop_data[i].mjd;
        X[i]    = eop_data[i].x;        // arcsec
        Y[i]    = eop_data[i].y;        // arcsec
        UT1[i]  = eop_data[i].ut1_utc;  // sec
        DPSI[i] = eop_data[i].dpsi;     // arcsec
        DEPS[i] = eop_data[i].deps;     // arcsec
    }

    // lagint(узлы, значения, N, точка): RJD — независимая переменная.
    double x_int   = lagint(RJD.data(), X.data(),   N, rjd_int);
    double y_int   = lagint(RJD.data(), Y.data(),   N, rjd_int);
    double ut1_int = lagint(RJD.data(), UT1.data(), N, rjd_int);

    // Океанический эффект (x, y, UT1)
    double cor_x, cor_y, cor_ut1, cor_lod;
    pmut1_oceans(rjd_int, cor_x, cor_y, cor_ut1, cor_lod);
    x_int   += cor_x;
    y_int   += cor_y;
    ut1_int += cor_ut1;

    // Лунно-солнечный эффект (x, y)
    pm_gravi(rjd_int, cor_x, cor_y);
    x_int += cor_x;
    y_int += cor_y;

    eop_int(0) = ut1_int; // UT1-UTC [сек]
    eop_int(1) = x_int;   // x [угл. сек]
    eop_int(2) = y_int;   // y [угл. сек]
    // В оригинале dpsi/deps здесь не заполняются; интерполируем для полноты.
    eop_int(3) = lagint(RJD.data(), DPSI.data(), N, rjd_int);
    eop_int(4) = lagint(RJD.data(), DEPS.data(), N, rjd_int);
}

} // namespace ariadna
