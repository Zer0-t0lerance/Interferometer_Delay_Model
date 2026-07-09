// cfx_parser.cpp
//
// Разбор файла задания коррелятора (.cfx) в структуру CfxTask.
// Формат — INI-подобные блоки:
//   [$TLSC] ... [$end]      — станция (наземная: TLSC_PAR = X,Y,Z,VX,VY,VZ,AXOFF,MJD,MOUNT;
//                             космическая RASTRON: ORB_FILE + TIMEOFS вместо координат);
//   [$SOURCE] ... [$END]    — источник (RA, DEC в ГРАДУСАХ);
//   [$skan] ... [$end]      — скан (start = <дата> , <длит>s; source; telescopes = iam,...);
//   [$OUTPAR] ... [$END]    — выходные параметры (%W — папка полиномов).
// Даты вида 23d04m2014y13h00m00s.

#include "functions.h"
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cctype>
#include <algorithm>

namespace ariadna {

static int ymd_to_mjd_cfx(int y, int m, int d) {
    long a = (14 - m) / 12, yy = y + 4800 - a, mm = m + 12 * a - 3;
    long jdn = d + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045;
    return static_cast<int>(jdn - 2400001);
}

// "23d04m2014y13h00m00s" -> mjd, utc(доля суток).
static bool parse_cfx_datetime(const std::string& s, int& mjd, double& utc) {
    int D, M, Y, h, mi, sec;
    if (std::sscanf(s.c_str(), "%dd%dm%dy%dh%dm%ds", &D, &M, &Y, &h, &mi, &sec) != 6) return false;
    mjd = ymd_to_mjd_cfx(Y, M, D);
    utc = (h * 3600.0 + mi * 60.0 + sec) / cnst::SECDAY;
    return true;
}

static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    return (a == std::string::npos) ? "" : s.substr(a, b - a + 1);
}

// Значение после первого '=' (ключ может содержать пробелы: "CLOCK DELAY = ...").
static std::string value_after_eq(const std::string& line) {
    size_t p = line.find('=');
    return (p == std::string::npos) ? "" : trim(line.substr(p + 1));
}

bool parse_cfx(const std::string& path, CfxTask& task) {
    task = CfxTask{};
    std::ifstream f(path);
    if (!f) { std::fprintf(stderr, "parse_cfx: не открыть %s\n", path.c_str()); return false; }

    std::string line;
    enum { NONE, TLSC, SOURCE, SKAN, OUTPAR, OTHER } blk = NONE;
    CfxStation st; CfxSource src; CfxScan sc;

    auto lc = [](std::string s){ std::transform(s.begin(), s.end(), s.begin(), ::tolower); return s; };

    while (std::getline(f, line)) {
        std::string t = trim(line);
        if (t.empty() || t[0] == '#') continue;

        // Начало/конец блоков.
        std::string tl = lc(t);
        if (tl == "[$tlsc]")      { blk = TLSC;   st = CfxStation{}; continue; }
        if (tl == "[$source]")    { blk = SOURCE; src = CfxSource{}; continue; }
        if (tl == "[$skan]")      { blk = SKAN;   sc = CfxScan{};    continue; }
        if (tl == "[$outpar]")    { blk = OUTPAR; continue; }
        if (tl == "[$clock]")     { blk = OTHER;  continue; }
        if (tl == "[$end]") {
            if (blk == TLSC)   task.stations.push_back(st);
            if (blk == SOURCE) task.sources.push_back(src);
            if (blk == SKAN)   task.scans.push_back(sc);
            blk = NONE; continue;
        }

        std::string key = trim(t.substr(0, t.find('=')));
        std::string kl = lc(key);
        std::string val = value_after_eq(t);

        if (blk == TLSC) {
            if (kl == "name") st.name = val;
            else if (kl == "iam_name") st.iam = val;
            else if (kl == "poly_file") { size_t c = val.find(':'); st.poly_file = (c==std::string::npos)?val:trim(val.substr(c+1)); }
            else if (kl == "orb_file") { st.orb_file = val; st.is_space = true; }
            else if (kl == "tlsc_par") {
                std::vector<double> v; std::string tok; std::stringstream ss(val);
                std::string mount;
                std::vector<std::string> fields;
                while (std::getline(ss, tok, ',')) fields.push_back(trim(tok));
                // X,Y,Z,VX,VY,VZ,AXOFF,MJD,MOUNT
                if (fields.size() >= 9) {
                    st.xyz << std::stod(fields[0]), std::stod(fields[1]), std::stod(fields[2]);
                    st.vel << std::stod(fields[3]), std::stod(fields[4]), std::stod(fields[5]);
                    st.axoff = std::stod(fields[6]);
                    st.epoch_mjd = std::stod(fields[7]);
                    st.mount = fields[8];
                }
            }
            else if (kl == "clock delay") { try { st.clock_delay = std::stod(val); } catch (...) {} }
            else if (kl == "clock rate")  { try { st.clock_rate  = std::stod(val); } catch (...) {} }
        } else if (blk == SOURCE) {
            if (kl == "name") src.name = val;
            else if (kl == "ra")  src.ra  = std::stod(val) * cnst::CDEGRAD; // градусы -> рад
            else if (kl == "dec") src.dec = std::stod(val) * cnst::CDEGRAD;
        } else if (blk == SKAN) {
            if (kl == "start") {
                // "23d04m2014y13h00m00s , 870s"
                size_t comma = val.find(',');
                std::string dt = trim(val.substr(0, comma));
                std::string dur = (comma==std::string::npos) ? "" : trim(val.substr(comma+1));
                parse_cfx_datetime(dt, sc.mjd, sc.utc);
                if (!dur.empty()) { size_t sp = dur.find('s'); sc.dur_sec = std::stod(dur.substr(0, sp)); }
            } else if (kl == "source") sc.source = val;
            else if (kl == "telescopes") {
                std::stringstream ss(val); std::string tok;
                while (std::getline(ss, tok, ',')) { std::string x = trim(tok); if (!x.empty()) sc.tel_iam.push_back(x); }
            }
        } else if (blk == OUTPAR) {
            if (key == "%W") task.out_path = val;
        }
    }
    return !task.stations.empty() && !task.scans.empty();
}

} // namespace ariadna
