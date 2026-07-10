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
#include <set>
#include <cctype>

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

// SitePrep наземной станции на эпоху: координаты + дрейф + физика из каталогов
// (океан/атм нагрузка, термопараметры). phys — Station с заполненными tide_data/
// atm_load/def_par (из каталогов); солидные приливы и прилив полюса считаются всегда.
static SitePrep siteprep_ground(const CfxStation& st, double epoch_mjd_utc, const Station& phys) {
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
    // Физика из каталогов: океаническая/атмосферная нагрузка + термопараметры.
    s.tide_data = phys.tide_data;   // океаническая нагрузка (SITE_TIDE_OC)
    s.atm_load = phys.atm_load;     // атмосферная нагрузка (SITE_ATM40)
    s.def_par = phys.def_par;       // термодеформация (THERM_DEF40)
    s.axsty = st.mount; s.offs = st.axoff;
    // Метео в задании нет: давление = опорное p_0 (регрессия по давлению = 0, но
    // гармоническая атм.нагрузка работает); температура = опорная t_0 (dT=0 -> термо = 0).
    s.pres = (phys.atm_load.p_0 > 0.0) ? phys.atm_load.p_0 : 1013.25;
    s.tC = phys.def_par.t_0;
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
                                 const std::vector<SpaceStation>& orbit, bool with_tropo,
                                 const Station& phys) {
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

    // Сплайны орбиты (пол/скор/уск) строим ОДИН РАЗ до цикла — иначе 39601 точка × сотни
    // вызовов orbit_interp (сильно медленно). Ось — секунды от старта орбиты (точность).
    tk::spline osx, osy, osz; double orb_t0 = 0.0; bool orb_ok = false;
    if (st.is_space && orbit.size() >= 2) {
        int n = static_cast<int>(orbit.size());
        orb_t0 = static_cast<double>(orbit[0].mjd) + orbit[0].utc;
        std::vector<double> t(n), px(n), py(n), pz(n);
        for (int i = 0; i < n; ++i) {
            t[i] = ((static_cast<double>(orbit[i].mjd) - orbit[0].mjd) + (orbit[i].utc - orbit[0].utc)) * cnst::SECDAY;
            px[i] = orbit[i].xyz.x() * 1000.0; py[i] = orbit[i].xyz.y() * 1000.0; pz[i] = orbit[i].xyz.z() * 1000.0;
        }
        osx.set_points(t, px); osy.set_points(t, py); osz.set_points(t, pz); orb_ok = true;
    }

    // SitePrep станции на момент: наземная (координаты+дрейф) или космическая (орбита).
    auto make_siteprep = [&](int m, double u) -> SitePrep {
        if (st.is_space) {
            SitePrep sp; sp.is_space = true;
            if (orb_ok) {
                double q = (m + u - orb_t0) * cnst::SECDAY; // сек от старта орбиты
                sp.x_orbit << osx(q), osy(q), osz(q);
                sp.v_orbit << osx.deriv(1, q), osy.deriv(1, q), osz.deriv(1, q);
                sp.a_orbit << osx.deriv(2, q), osy.deriv(2, q), osz.deriv(2, q);
            }
            return sp;
        }
        return siteprep_ground(st, m + u, phys);
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

// Координаты пункта приёма (ITRF XYZ, м) по имени из каталога ITRF2005 (geodetic->cartesian).
// Формат строки: DOMES NAME FULLNAME ... описание ... LON(d m s) LAT(d m s) HEIGHT.
static bool reception_itrf(const std::string& cat, const std::string& name, Eigen::Vector3d& xyz) {
    std::ifstream f(cat); if (!f) return false;
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream ss(line); std::string domes, nm;
        if (!(ss >> domes >> nm)) continue;
        if (nm != name) continue;
        // Последние 7 числовых токенов строки: lon d m s, lat d m s, height.
        std::vector<double> nums; std::istringstream s2(line); std::string tok;
        while (s2 >> tok) { try { size_t p; double v = std::stod(tok, &p); if (p == tok.size()) nums.push_back(v); } catch (...) {} }
        if (nums.size() < 7) return false;
        size_t k = nums.size() - 7;
        double lon = (nums[k] + nums[k + 1] / 60.0 + nums[k + 2] / 3600.0) * cnst::CDEGRAD;
        double lat = (nums[k + 3] + nums[k + 4] / 60.0 + nums[k + 5] / 3600.0) * cnst::CDEGRAD;
        double h = nums[k + 6];
        double a = cnst::AE, fl = 1.0 / cnst::F, e2 = fl * (2.0 - fl);
        double N = a / std::sqrt(1.0 - e2 * std::sin(lat) * std::sin(lat));
        xyz << (N + h) * std::cos(lat) * std::cos(lon), (N + h) * std::cos(lat) * std::sin(lon),
               (N * (1.0 - e2) + h) * std::sin(lat);
        return true;
    }
    return false;
}

// Григорианская дата -> MJD на 0h.
static int dp_ymd_to_mjd(int y, int m, int d) {
    long a = (14 - m) / 12, yy = y + 4800 - a, mm = m + 12 * a - 3;
    long jdn = d + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045;
    return static_cast<int>(jdn - 2400001);
}

// Время старта сегмента из имени файла данных (кодировка YYYYDDDHHMMSS, напр.
// PUSH_2014113130000). Возвращает MJD+доля суток; false если не распарсить.
static bool parse_file_mjd_utc(const std::string& fileval, double& mjd_utc) {
    for (size_t i = 0; i < fileval.size();) {
        if (std::isdigit(static_cast<unsigned char>(fileval[i]))) {
            size_t j = i; while (j < fileval.size() && std::isdigit(static_cast<unsigned char>(fileval[j]))) ++j;
            if (j - i >= 13) {
                std::string d = fileval.substr(i, 13);
                int Y = std::stoi(d.substr(0, 4)), DDD = std::stoi(d.substr(4, 3));
                int HH = std::stoi(d.substr(7, 2)), MM = std::stoi(d.substr(9, 2)), SS = std::stoi(d.substr(11, 2));
                mjd_utc = dp_ymd_to_mjd(Y, 1, 1) + (DDD - 1) + (HH * 3600.0 + MM * 60.0 + SS) / 86400.0;
                return true;
            }
            i = j;
        } else ++i;
    }
    return false;
}

// Даунлинк-задержка космос->пункт приёма на момент mjd_utc (пункт приёма движется с Землёй).
static double timeofs_at(double mjd_utc, const std::vector<SpaceStation>& orbit,
                         const std::vector<EOPData>& eop, const Eigen::Vector3d& recv_itrf) {
    int m = static_cast<int>(std::floor(mjd_utc)); double u = mjd_utc - m;
    EpochEnv e = prep_epoch(m, u, eop);
    Eigen::Vector3d xr, vr, aa; orbit_interp(orbit, mjd_utc, xr, vr, aa);
    Eigen::Vector3d recv_j2000 = e.R * recv_itrf; // ITRF -> J2000 на время t (учёт вращения Земли)
    return -((xr - recv_j2000).norm()) / cnst::C;
}

// Перезапись задания cfx с пересчитанными TIMEOFS (даунлинк космос->пункт приёма).
// TIMEOFS ГЕНЕРИРУЮТСЯ из времени файлов данных (FILExx, кодировка YYYYDDDHHMMSS) — работает
// и «с нуля» (когда строк TIMEOFS в задании нет). Старые строки TIMEOFS заменяются.
void write_timeofs_cfx(const std::string& cfx_in, const std::string& cfx_out,
                       const std::vector<SpaceStation>& orbit, const std::vector<EOPData>& eop,
                       const Eigen::Vector3d& recv_itrf, const std::string& space_name) {
    std::ifstream fin(cfx_in); std::ofstream fout(cfx_out);
    if (!fin || !fout) { std::fprintf(stderr, "write_timeofs_cfx: ошибка файла\n"); return; }
    std::string line; bool in_tlsc = false, is_space = false; int count = 0;
    std::set<std::string> written; // индексы TIMEOFS, уже сгенерированные (после FILE)
    while (std::getline(fin, line)) {
        std::string t = line; size_t a = t.find_first_not_of(" \t\r\n"); std::string tl = (a == std::string::npos) ? "" : t.substr(a);
        if (tl.rfind("[$TLSC]", 0) == 0) { in_tlsc = true; is_space = false; }
        else if (tl.rfind("[$end]", 0) == 0 || tl.rfind("[$END]", 0) == 0) { in_tlsc = false; is_space = false; }
        else if (in_tlsc && tl.rfind("name", 0) == 0 && tl.find(space_name) != std::string::npos) is_space = true;
        else if (in_tlsc && tl.rfind("ORB_FILE", 0) == 0) is_space = true;

        // FILExx в космоблоке: копируем строку и СРАЗУ пишем сгенерированный TIMEOFSxx.
        if (in_tlsc && is_space && tl.rfind("FILE", 0) == 0 && tl.size() > 4 && std::isdigit((unsigned char)tl[4])) {
            fout << line << "\n";
            std::string idx; for (size_t k = 4; k < tl.size() && std::isdigit((unsigned char)tl[k]); ++k) idx += tl[k];
            std::string val = tl.substr(tl.find('=') + 1);
            double mjd_utc;
            if (parse_file_mjd_utc(val, mjd_utc)) {
                double tof = timeofs_at(mjd_utc, orbit, eop, recv_itrf);
                char buf[160]; std::snprintf(buf, sizeof(buf), "   TIMEOFS%s = %.15e, %.15e", idx.c_str(), tof, mjd_utc);
                fout << buf << "\n"; written.insert(idx); ++count;
            }
            continue;
        }
        // Старые TIMEOFS: если уже сгенерировали для этого индекса — пропускаем (заменили).
        if (in_tlsc && is_space && tl.rfind("TIMEOFS", 0) == 0) {
            std::string idx; for (size_t k = 7; k < tl.size() && std::isdigit((unsigned char)tl[k]); ++k) idx += tl[k];
            if (written.count(idx)) continue;
            fout << line << "\n"; // не смогли сгенерировать -> оставляем как есть
            continue;
        }
        fout << line << "\n";
    }
    std::printf("  TIMEOFS сгенерировано: %d -> %s (пункт приёма |R_ITRF|=%.1f км, повёрнут в J2000 на каждый момент)\n",
                count, cfx_out.c_str(), recv_itrf.norm() / 1e3);
}

void process_task(const std::string& cfx_path, const std::string& orbit_path,
                  const std::string& out_dir, const std::string& eop_path,
                  double block_sec, int degree, double sample_sec, bool with_tropo,
                  const std::string& recv_name) {
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

    // Физика станций из каталогов (та же папка, что EOP): океан/атм нагрузка + термопараметры.
    std::string catdir = eop_path;
    size_t sl = catdir.find_last_of("/\\"); catdir = (sl == std::string::npos) ? "" : catdir.substr(0, sl + 1);
    std::vector<Station> phys(task.stations.size());
    for (size_t i = 0; i < task.stations.size(); ++i) phys[i].name = task.stations[i].name;
    {
        char p_oc[256], p_at[256], p_an[256];
        std::snprintf(p_oc, 256, "%sVLBI_ocload_40.cat", catdir.c_str());
        std::snprintf(p_at, 256, "%sVLBI_atmload4_12.cat", catdir.c_str());
        std::snprintf(p_an, 256, "%santenna-info.cat", catdir.c_str());
        std::vector<oc_record> oc; std::vector<atm_record> at; std::vector<ant_record> an;
        ReadOC(oc, p_oc);        map_ocean_tides_to_stations(oc, phys);
        ReadATM(at, p_at);       map_atm_loading_to_stations(at, phys);
        ReadANT_INFO(an, p_an);  build_def_par_from_ant_info(an, phys);
    }

    // Каждой станции — свой файл полиномов.
    for (size_t i = 0; i < task.stations.size(); ++i) {
        const CfxStation& st = task.stations[i];
        StationPoly poly = compute_station_poly(st, src, mjd0, utc0, dur, eop, block_sec, degree, sample_sec, orbit, with_tropo, phys[i]);
        std::string out = out_dir; if (!out.empty() && out.back() != '/' && out.back() != '\\') out += "/";
        out += st.poly_file;
        write_station_poly(out, poly);
        std::printf("  %-9s -> %s (%zu блоков)\n", st.name.c_str(), out.c_str(), poly.blocks.size());
    }

    // TIMEOFS (даунлинк космос->пункт приёма): пересчитать и записать новый *_p.cfx.
    if (has_space && !orbit.empty()) {
        char itrf[256]; std::snprintf(itrf, 256, "%sITRF2005_2.CAT", catdir.c_str());
        Eigen::Vector3d recv;
        std::string space_name; for (const auto& s : task.stations) if (s.is_space) space_name = s.name;
        if (reception_itrf(itrf, recv_name, recv)) {
            std::string cfx_out = cfx_path;
            size_t dot = cfx_out.rfind(".cfx");
            cfx_out = (dot == std::string::npos) ? (cfx_out + "_p.cfx") : (cfx_out.substr(0, dot) + "_p.cfx");
            write_timeofs_cfx(cfx_path, cfx_out, orbit, eop, recv, space_name);
        } else {
            std::fprintf(stderr, "process_task: пункт приёма '%s' не найден в %s (TIMEOFS пропущены)\n", recv_name.c_str(), itrf);
        }
    }
}

} // namespace ariadna
