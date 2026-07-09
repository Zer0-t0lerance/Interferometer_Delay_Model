// delay_poly.cpp
//
// Слой полиномов задержки для коррелятора: считает геоцентрическую вакуумную задержку
// станции на сетке времён сеанса, сшивает кубическим сплайном и раскладывает по блокам
// в полиномы (P0..P_degree) — формат эталонных .TXT (см. example/).
//
// Модель полинома (установлена сверкой с example/BADARY.TXT до 0.3 м): ВАКУУМНАЯ
// геометрия + релятивизм, задержка станции ОТНОСИТЕЛЬНО ЦЕНТРА ЗЕМЛИ, БЕЗ тропосферы и
// БЕЗ clock. Геоцентр моделируется как station1 (is_space в r=0), станция — station2.
// Сшивка: единый кубический сплайн по всей сетке -> C2-непрерывность на стыках блоков.

#include "functions.h"
#include "catalog_bridge.h"
#include "READ_CAT.h"
#include "../external/spline.h"
#include <fstream>
#include <cstdio>
#include <cmath>
#include <cstring>

namespace ariadna {

// --- Среда на эпоху (время, EOP, эфемериды, ориентация) ---
struct EpochEnv {
    int mjd = 0; double utc = 0, jd0 = 0, ct = 0, cent = 0, ut1_sec = 0, xp = 0, yp = 0;
    Eigen::VectorXd f, fd; Eigen::Vector2d gast;
    Eigen::Matrix3d Earth, Sun, Moon, R, dR, d2R;
    Eigen::Matrix<double, 3, 2> sun_geo, moon_geo;
};

static EpochEnv prep_epoch(int mjd, double utc, const std::vector<EOPData>& eop) {
    EpochEnv e; e.mjd = mjd; e.utc = utc; e.jd0 = cnst::MJD_OFFSET + mjd;
    double TAI, TT; tai_time(mjd, utc, TAI, TT);
    Observation o{}; o.mjd = mjd; o.utc = utc;
    double ut1_out; Eigen::VectorXd eop_int(5), deop_int(5), arg(8); Eigen::MatrixXd dd(3, 2), dl(3, 2);
    interp_eop(0, o, TT, ut1_out, eop_int, deop_int, arg, dd, dl, eop);
    e.xp = eop_int(1); e.yp = eop_int(2); double ut1_utc = eop_int(0);
    double ut1_frac = utc + ut1_utc / cnst::SECDAY;
    e.ct = t_eph40(e.jd0, TT, ut1_frac, 0, 0, 0);
    double jd_tt = e.jd0 + TT, jd_ut1 = e.jd0 + ut1_frac;
    e.cent = (jd_tt - cnst::JD2000) / cnst::JUL_CENT;
    e.ut1_sec = utc * cnst::SECDAY + ut1_utc;
    get_celestial_bodies(e.jd0, e.ct, e.Earth, e.Sun, e.Moon);
    Eigen::Matrix3d Eg; Eigen::MatrixXd sunM, moonM; jpl_eph(e.jd0, e.ct, Eg, sunM, moonM);
    e.sun_geo = sunM; e.moon_geo = moonM;
    get_r2000_matrices(jd_tt, jd_ut1, e.xp, e.yp, e.R, e.dR, e.d2R);
    e.f.resize(5); e.fd.resize(5); double c2; fund_arg(e.jd0, e.ct, c2, e.f, e.fd);
    e.gast = gast_iau2006(jd_tt, jd_ut1);
    return e;
}

// SitePrep наземной станции на эпоху: координаты + дрейф от эпохи каталога.
static SitePrep siteprep_ground(const CfxStation& st, double epoch_mjd_utc) {
    SitePrep s;
    double drift_epoch = (st.epoch_mjd > 1.0) ? st.epoch_mjd : cnst::MJD_J2000;
    double dyear = (epoch_mjd_utc - drift_epoch) / cnst::DAYS_PER_YEAR;
    s.xsta_itrf = st.xyz + st.vel * dyear;
    double r = s.xsta_itrf.norm();
    s.lat_gcen = std::asin(s.xsta_itrf.z() / r);
    s.lon_gcen = std::atan2(s.xsta_itrf.y(), s.xsta_itrf.x()); if (s.lon_gcen < 0) s.lon_gcen += cnst::TWOPI;
    double req = std::sqrt(s.xsta_itrf.x() * s.xsta_itrf.x() + s.xsta_itrf.y() * s.xsta_itrf.y());
    GEOD(req, s.xsta_itrf.z(), s.lat_geod, s.h_geod);
    double cla = std::cos(s.lat_geod), sla = std::sin(s.lat_geod), clo = std::cos(s.lon_gcen), slo = std::sin(s.lon_gcen);
    s.vw_i.col(0) << cla * clo, cla * slo, sla;
    s.vw_i.col(1) << -slo, clo, 0.0;
    s.vw_i.col(2) << -sla * clo, -sla * slo, cla;
    s.tide_data.amplitudes.setZero(); s.tide_data.phases.setZero(); s.atm_load.coef.setZero();
    s.axsty = st.mount; s.offs = st.axoff; s.pres = 1013.25; s.tC = 0.0;
    return s;
}

// Геоцентрическая задержка станции на момент. with_tropo=false — вакуумная геометрия
// (воспроизведение упрощённого эталона); true — с тропосферой (стандартная атмосфера,
// т.к. метео в задании нет). Солидные земные приливы и прилив полюса считаются всегда.
static double geocentric_delay(const SitePrep& station, const Eigen::Vector3d& K,
                               const EpochEnv& e, bool with_tropo) {
    SitePrep geo; geo.is_space = true;
    Observation obs{}; obs.sta1 = 0; obs.sta2 = 1;
    obs.p1 = 1013.25; obs.t1 = 0.0; obs.e1 = 50.0; obs.p2 = 1013.25; obs.t2 = 0.0; obs.e2 = 50.0;
    double tau, dtau;
    compute_delay_obs(geo, station, K, obs, e.mjd, e.utc, e.jd0, e.ct, e.cent, e.ut1_sec,
                      e.f, e.fd, e.gast, e.Earth, e.Sun, e.Moon, e.sun_geo, e.moon_geo,
                      e.xp, e.yp, 0, 0, e.R, e.dR, e.d2R, tau, dtau, nullptr, with_tropo);
    return tau;
}

// МНК-полином степени deg через точки (x,y). Коэффициенты P_k при x^k по возрастанию.
// x НОРМИРУЕТСЯ (u = x/xs) перед подгонкой — иначе Вандермонд с x^5 до ~1e9 чудовищно
// обусловлен и P0 теряет точность (заметно для большой задержки, напр. RASTRON).
// Затем P_k = c_k / xs^k возвращает коэффициенты в исходном базисе x.
static std::vector<double> polyfit_local(const std::vector<double>& x, const std::vector<double>& y, int deg) {
    int n = static_cast<int>(x.size()), d = deg; if (d > n - 1) d = n - 1;
    double xs = 0; for (double v : x) xs = std::max(xs, std::fabs(v));
    if (xs <= 0) xs = 1.0;
    Eigen::MatrixXd A(n, d + 1); Eigen::VectorXd b(n);
    for (int i = 0; i < n; ++i) { double u = x[i] / xs, p = 1; for (int k = 0; k <= d; ++k) { A(i, k) = p; p *= u; } b(i) = y[i]; }
    Eigen::VectorXd c = A.householderQr().solve(b);
    std::vector<double> out(deg + 1, 0.0);
    double s = 1.0; for (int k = 0; k <= d; ++k) { out[k] = c(k) / s; s *= xs; }
    return out;
}

// Полиномы задержки станции по сеансу [t_start, t_end] (MJD): сетка -> сплайн -> блоки.
StationPoly compute_station_poly(const CfxStation& st, const CfxSource& src,
                                 int mjd0, double utc0, double dur_sec,
                                 const std::vector<EOPData>& eop,
                                 double block_sec, int degree, double sample_sec,
                                 const std::vector<SpaceStation>& orbit, bool with_tropo) {
    StationPoly poly; poly.telescope = st.name; poly.source = src.name; poly.order = degree + 1;

    // Направление на источник (в J2000; для фикс. источника постоянно).
    std::vector<Source> ss(1); ss[0].ra = src.ra; ss[0].dec = src.dec;
    std::vector<Eigen::Vector3d> kv; source_vec(ss, mjd0 + utc0, kv);
    const Eigen::Vector3d K = kv[0];

    // Момент t (сек от начала) -> (mjd, utc).
    auto to_mjd_utc = [&](double t_sec, int& mjd, double& utc) {
        double tot = utc0 + t_sec / cnst::SECDAY;
        int add = static_cast<int>(std::floor(tot));
        mjd = mjd0 + add; utc = tot - add;
    };

    // SitePrep станции на момент: наземная (координаты+дрейф) или космическая (орбита).
    auto make_siteprep = [&](int m, double u) -> SitePrep {
        if (st.is_space) {
            SitePrep sp; sp.is_space = true;
            orbit_interp(orbit, m + u, sp.x_orbit, sp.v_orbit, sp.a_orbit);
            return sp;
        }
        return siteprep_ground(st, m + u);
    };

    // 1) Сетка задержек по сеансу (+ запас на края для сплайна).
    std::vector<double> xs, ys;
    for (double t = -sample_sec; t <= dur_sec + sample_sec + 1e-6; t += sample_sec) {
        int m; double u; to_mjd_utc(t, m, u);
        EpochEnv e = prep_epoch(m, u, eop);
        SitePrep sp = make_siteprep(m, u);
        xs.push_back(t); ys.push_back(geocentric_delay(sp, K, e, with_tropo));
    }
    tk::spline sp; sp.set_points(xs, ys);

    // 2) Блоки: на каждый — МНК-полином степени degree по сплайну (аргумент = сек от начала блока).
    int nblk = static_cast<int>(std::ceil(dur_sec / block_sec - 1e-9));
    const int fine = 12; // точек внутри блока для МНК
    for (int b = 0; b < nblk; ++b) {
        double t0 = b * block_sec, t1 = std::min((b + 1) * block_sec, dur_sec);
        std::vector<double> bx, by;
        for (int i = 0; i <= fine; ++i) {
            double tt = t0 + (t1 - t0) * i / fine;
            bx.push_back(tt - t0); by.push_back(sp(tt));
        }
        DelayPolyBlock blk;
        int m; double u; to_mjd_utc(t0, m, u); blk.mjd = m; blk.utc_start = u;
        int m2; double u2; to_mjd_utc(t1, m2, u2); blk.utc_stop = u2;
        blk.coef = polyfit_local(bx, by, degree);
        poly.blocks.push_back(blk);
    }
    return poly;
}

// Запись полиномов станции в .TXT (формат эталона example/).
void write_station_poly(const std::string& path, const StationPoly& poly) {
    std::ofstream f(path);
    if (!f) { std::fprintf(stderr, "write_station_poly: не открыть %s\n", path.c_str()); return; }
    f << "telescope = " << poly.telescope << "\n";
    f << "order = " << poly.order << "\n\n";
    auto wr_time = [&](std::FILE*, int mjd, double utc) {};
    (void)wr_time;
    for (const auto& b : poly.blocks) {
        // Дата/время из mjd+utc.
        int mjd = b.mjd; double utc = b.utc_start;
        // MJD -> Y/M/D (обратный Fliegel).
        long jdn = mjd + 2400001; long a = jdn + 32044, bb = (4 * a + 3) / 146097, c = a - 146097 * bb / 4;
        long dd = (4 * c + 3) / 1461, ee = c - 1461 * dd / 4, mm = (5 * ee + 2) / 153;
        int day = static_cast<int>(ee - (153 * mm + 2) / 5 + 1);
        int mon = static_cast<int>(mm + 3 - 12 * (mm / 10));
        int yr = static_cast<int>(100 * bb + dd - 4800 + mm / 10);
        auto hms = [](double u, int& h, int& mi, int& s) { double t = u * 86400.0; h = (int)(t / 3600); mi = (int)((t - h * 3600) / 60); s = (int)std::round(t - h * 3600 - mi * 60); };
        int h1, mi1, s1, h2, mi2, s2; hms(b.utc_start, h1, mi1, s1); hms(b.utc_stop, h2, mi2, s2);
        char buf[128];
        f << "source = " << poly.source << "\n";
        std::snprintf(buf, sizeof(buf), "start = %02d/%02d/%04d %02dh%02dm%02ds\n", day, mon, yr, h1, mi1, s1); f << buf;
        std::snprintf(buf, sizeof(buf), "stop = %02d/%02d/%04d %02dh%02dm%02ds\n", day, mon, yr, h2, mi2, s2); f << buf;
        for (size_t k = 0; k < b.coef.size(); ++k) {
            std::snprintf(buf, sizeof(buf), "P%zu = %.14e\n", k, b.coef[k]); f << buf;
        }
    }
}

void process_task(const std::string& cfx_path, const std::string& orbit_path,
                  const std::string& out_dir, const std::string& eop_path,
                  double block_sec, int degree, double sample_sec, bool with_tropo) {
    CfxTask task;
    if (!parse_cfx(cfx_path, task)) { std::fprintf(stderr, "process_task: не разобрать %s\n", cfx_path.c_str()); return; }

    // Орбита космической станции (если есть).
    std::vector<SpaceStation> orbit;
    bool has_space = false; for (const auto& s : task.stations) if (s.is_space) has_space = true;
    if (has_space && !orbit_path.empty()) read_scf_orbit(orbit_path, orbit);

    // Границы сеанса и источник (первого скана; предполагаем один источник на сеанс).
    if (task.scans.empty()) { std::fprintf(stderr, "process_task: нет сканов\n"); return; }
    int mjd0 = task.scans.front().mjd; double utc0 = task.scans.front().utc;
    double t_end = 0; for (const auto& sc : task.scans) {
        double e = (sc.mjd - mjd0) * 86400.0 + sc.utc * 86400.0 + sc.dur_sec; if (e > t_end) t_end = e;
    }
    double dur = t_end - utc0 * 86400.0;
    std::string src_name = task.scans.front().source;
    CfxSource src; for (const auto& s : task.sources) if (s.name == src_name) src = s;

    // EOP: 7 узлов вокруг сеанса из каталога EOPC04.
    std::vector<eop_record> raw_eop; char epath[256]; std::strncpy(epath, eop_path.c_str(), 255); epath[255] = 0;
    std::vector<EOPData> eop;
    if (ReadEOP(raw_eop, epath) == 0 || !raw_eop.empty()) select_eop_nodes(raw_eop, mjd0, cnst::EOP_NDATA, eop);
    if (eop.size() < (size_t)cnst::EOP_NDATA) { std::fprintf(stderr, "process_task: мало узлов EOP\n"); return; }

    // Каждой станции — свой файл полиномов.
    for (const auto& st : task.stations) {
        StationPoly poly = compute_station_poly(st, src, mjd0, utc0, dur, eop, block_sec, degree, sample_sec, orbit, with_tropo);
        std::string out = out_dir; if (!out.empty() && out.back() != '/' && out.back() != '\\') out += "/";
        out += st.poly_file;
        write_station_poly(out, poly);
        std::printf("  %-9s -> %s (%zu блоков)\n", st.name.c_str(), out.c_str(), poly.blocks.size());
    }
}

} // namespace ariadna
