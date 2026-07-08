// dmeteo.cpp
//
// Порт DMETEO1_DT / DMETEO2_DT — аппроксимация метеопараметров станции (температура,
// давление, влажность) полиномом по времени и вычисление их производных на момент
// наблюдения. В оригинале ARIADNA полином подгоняется по методу наименьших квадратов
// (POLCOEF) степени NDEG; производные dTdt/dPdt/dHumdt идут в SITE_ATM40 (dPdt),
// THERM_DEF40 (dTdt) и тропосферу — как СКОРОСТИ (влияют на скорость задержки dtau,
// не на саму задержку).
//
//   dmeteo1_dt: для каждой станции собирает метео-точки всех её наблюдений и строит
//               МНК-полином степени ndeg в аргументе (t - t_mean) [сутки];
//   dmeteo2_dt: возвращает производную полинома на момент (mjd+utc) — ед./СУТКИ.
//
// Единицы производных — на СУТКИ (как ждут SITE_ATM40 dPdt и THERM_DEF40 dTdt).
// Оригинального DMETEO*.for в архиве нет; реализация по согласованной модели
// (МНК-полином степени NDEG), степень задаётся вызывающим (обычно 2).

#include "functions.h"

namespace ariadna {
namespace {

// МНК-подгонка полинома степени deg: p(x) = c0 + c1 x + ... + c_deg x^deg.
// Возвращает коэффициенты по возрастанию (deg+1). При нехватке точек степень снижается.
Eigen::VectorXd polyfit(const std::vector<double>& x, const std::vector<double>& y, int deg) {
    const int npts = static_cast<int>(x.size());
    int d = deg;
    if (d > npts - 1) d = npts - 1;    // не больше, чем позволяют точки
    if (d < 0) return Eigen::VectorXd::Zero(deg + 1);

    Eigen::MatrixXd A(npts, d + 1);
    Eigen::VectorXd b(npts);
    for (int i = 0; i < npts; ++i) {
        double p = 1.0;
        for (int k = 0; k <= d; ++k) { A(i, k) = p; p *= x[i]; }
        b(i) = y[i];
    }
    Eigen::VectorXd c = A.householderQr().solve(b);

    // Дополнить нулями до deg+1 (старшие коэффициенты, если степень снижена).
    Eigen::VectorXd out = Eigen::VectorXd::Zero(deg + 1);
    out.head(d + 1) = c;
    return out;
}

// Значение производной полинома p'(x) = c1 + 2 c2 x + 3 c3 x^2 + ...
double polyder_eval(const Eigen::VectorXd& c, double x) {
    double dp = 0.0, p = 1.0;
    for (int k = 1; k < c.size(); ++k) { dp += k * c(k) * p; p *= x; }
    return dp;
}

} // anonymous namespace

void dmeteo1_dt(const std::vector<Observation>& observations, int n_stations,
                double t_mean, int ndeg,
                std::vector<Eigen::VectorXd>& t_coef,
                std::vector<Eigen::VectorXd>& p_coef,
                std::vector<Eigen::VectorXd>& hum_coef) {
    t_coef.assign(n_stations, Eigen::VectorXd::Zero(ndeg + 1));
    p_coef.assign(n_stations, Eigen::VectorXd::Zero(ndeg + 1));
    hum_coef.assign(n_stations, Eigen::VectorXd::Zero(ndeg + 1));

    for (int i = 0; i < n_stations; ++i) {
        std::vector<double> x, T, P, H;
        for (const Observation& o : observations) {
            double xt = (static_cast<double>(o.mjd) + o.utc) - t_mean;
            if (o.sta1 == i) { x.push_back(xt); T.push_back(o.t1); P.push_back(o.p1); H.push_back(o.e1); }
            if (o.sta2 == i) { x.push_back(xt); T.push_back(o.t2); P.push_back(o.p2); H.push_back(o.e2); }
        }
        if (x.empty()) continue;               // станция не участвует -> нули (rate = 0)
        t_coef[i]   = polyfit(x, T, ndeg);
        p_coef[i]   = polyfit(x, P, ndeg);
        hum_coef[i] = polyfit(x, H, ndeg);
    }
}

void dmeteo2_dt(int ista, int ndeg, int mjd, double utc, double t_mean,
                const std::vector<Eigen::VectorXd>& t_coef,
                const std::vector<Eigen::VectorXd>& p_coef,
                const std::vector<Eigen::VectorXd>& hum_coef,
                double& dTdt, double& dPdt, double& dHumdt) {
    (void)ndeg;
    const double x = (static_cast<double>(mjd) + utc) - t_mean;
    dTdt   = polyder_eval(t_coef[ista],   x);   // °C/сут
    dPdt   = polyder_eval(p_coef[ista],   x);   // мбар/сут
    dHumdt = polyder_eval(hum_coef[ista], x);   // %/сут
}

} // namespace ariadna
