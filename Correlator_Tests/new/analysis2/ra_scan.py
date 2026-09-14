#!/usr/bin/env python
# ra_scan.py - корректный учёт откликов на РадиоАстроне.
#
# ПОЧЕМУ ОТДЕЛЬНЫЙ СЧЁТ. В выходном файле поиска отклика опорная антенна назначается
# СВОЯ В КАЖДОМ БЛОКЕ, и в разных блоках она разная. Запись "CLOCK = RA" в блоке -
# это поиск отклика на РадиоАстроне относительно опорной антенны ЭТОГО блока. Поэтому
# пара (опорная, RA) не является базой, которая прослеживается через весь сеанс:
# она возникает и исчезает вместе с назначением опорной. Более того, в 6 блоках
# опорная антенна в старом и новом прогоне РАЗНАЯ, и такие блоки при сравнении
# "по базам" вообще выпадают.
#
# Правильная единица счёта для космического плеча: СКАН. Вопрос ставится так -
# "был ли в этом скане обнаружен отклик на РадиоАстроне", без привязки к тому,
# какая антенна оказалась опорной. Ответ от выбора опорной не зависит.
#
# Выход: ra_scans.csv, REPORT_ra.txt

import re, os, glob, csv
from collections import defaultdict
import uvwlib

H = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(H)

re_hdr = re.compile(r'#\s*(\w+)-POL\s+Time interval\s*=\s*\[(\d+:\d+:\d+)\s*-\s*(\d+:\d+:\d+)\[\s*'
                    r'reference antenna:\s*(\S+)')
re_clk = re.compile(r'^(#?)\s*CLOCK\s*=\s*(\S+?),\s*(\S+?),\s*([-\d.eE+]+),\s*([-\d.eE+]+),\s*([-\d.eE+]+)')
re_snr = re.compile(r'#\s*snr\s*=\s*([-\d.eE+]+)\s*,\s*IF\s*=\s*(-?\d+)')


def parse(path):
    out = []; pol = ref = t0 = None
    for ln in open(path, encoding="utf-8", errors="replace").read().replace("\r\n", "\n").split("\n"):
        h = re_hdr.search(ln)
        if h:
            pol, t0, ref = h.group(1), h.group(2), h.group(4); continue
        c = re_clk.match(ln)
        if c:
            out.append(dict(pol=pol, ref=ref, t=t0, tfull=c.group(3), sta=c.group(2),
                            det=(c.group(1) != "#"), delay=float(c.group(4)),
                            rate=float(c.group(5)), snr=None, iff=None)); continue
        s = re_snr.search(ln)
        if s and out and out[-1]['snr'] is None:
            out[-1]['snr'] = float(s.group(1)); out[-1]['iff'] = int(s.group(2))
    return out


UV = uvwlib.build_index()
rows = []
for d in sorted(glob.glob(os.path.join(ROOT, "*"))):
    if not os.path.isdir(d) or os.path.basename(d).startswith("analysis"): continue
    exp = os.path.basename(d)
    for old in sorted(glob.glob(os.path.join(d, "*_ffr.txt"))):
        if old.endswith("_p_ffr.txt"): continue
        new = old[:-len("_ffr.txt")] + "_p_ffr.txt"
        if not os.path.exists(new): continue
        band = (re.search(r'_([CLKSXQ])_', os.path.basename(old)) or [None, '?'])[1]
        O, N = parse(old), parse(new)
        ra_o = {(r['pol'], r['t']): r for r in O if r['sta'] == 'RA' and r['iff'] != -100}
        ra_n = {(r['pol'], r['t']): r for r in N if r['sta'] == 'RA' and r['iff'] != -100}
        for k in sorted(set(ra_o) & set(ra_n)):
            o, n = ra_o[k], ra_n[k]
            t = uvwlib.parse_tfull(o['tfull'])
            bo = uvwlib.bproj(UV, exp, band, 'RA', o['ref'], t)
            bn = uvwlib.bproj(UV, exp, band, 'RA', n['ref'], t)
            rows.append(dict(exp=exp, band=band, pol=k[0], t=k[1], tfull=o['tfull'],
                             bproj_old=("%.4f" % bo) if bo else "",
                             bproj_new=("%.4f" % bn) if bn else "",
                             ref_old=o['ref'], ref_new=n['ref'],
                             ref_same=int(o['ref'] == n['ref']),
                             det_old=int(o['det']), det_new=int(n['det']),
                             snr_old=o['snr'], snr_new=n['snr'],
                             delay_old=o['delay'], delay_new=n['delay'],
                             rate_old=o['rate'], rate_new=n['rate']))

with open(os.path.join(H, "ra_scans.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

BUG = lambda r: r['exp'] == 'raks08ad' and r['band'] == 'C'   # блок с багом именования полиномов
D = [r for r in rows if not (r['exp'] == 'raks08ad' and r['band'] == 'C' and
                             (r['ref_old'] == 'Ef' or r['ref_new'] == 'Ef'))]

L = []; A = L.append
A("=" * 100)
A("ОТКЛИК НА РАДИОАСТРОНЕ, СЧЁТ ПО СКАНАМ")
A("Единица - скан сеанса в одной поляризации. Вопрос: обнаружен ли в этом скане отклик")
A("на РадиоАстроне (относительно той опорной антенны, которая назначена в этом блоке).")
A("=" * 100)
A("")
A("всего сканов с записью по РадиоАстрону: %d   (после исключения блока с багом: %d)"
  % (len(rows), len(D)))
nd = [r for r in D if not r['ref_same']]
A("сканов, где опорная антенна в старом и новом прогоне РАЗНАЯ: %d" % len(nd))
for r in nd:
    A("    %-9s %s %s %-8s  опорная %s -> %s   отклик %d -> %d   snr %.2f -> %.2f"
      % (r['exp'], r['band'], r['pol'], r['t'], r['ref_old'], r['ref_new'],
         r['det_old'], r['det_new'], r['snr_old'], r['snr_new']))
o = sum(r['det_old'] for r in D); n = sum(r['det_new'] for r in D)
lost = [r for r in D if r['det_old'] and not r['det_new']]
gain = [r for r in D if not r['det_old'] and r['det_new']]
A("")
A("ОТКЛИК НА РАДИОАСТРОНЕ: было %d -> стало %d   (%+d) из %d сканов" % (o, n, n - o, len(D)))
A("   в обеих моделях: %d | ПРОПАЛ: %d | ПОЯВИЛСЯ: %d | нет ни в одной: %d"
  % (sum(1 for r in D if r['det_old'] and r['det_new']), len(lost), len(gain),
     sum(1 for r in D if not r['det_old'] and not r['det_new'])))

RK = lambda r: (r['exp'], r['band'])
RUN = defaultdict(list)
for r in D: RUN[RK(r)].append(r)
tot = lambda a: (sum(x['det_old'] for x in a), sum(x['det_new'] for x in a))
rl = [k for k, v in RUN.items() if tot(v)[0] > 0 and tot(v)[1] == 0]
rg = [k for k, v in RUN.items() if tot(v)[0] == 0 and tot(v)[1] > 0]
r0 = [k for k, v in RUN.items() if tot(v)[0] == 0 and tot(v)[1] == 0]
A("")
A("ПО СЕАНСАМ (всего сеансов с записями по РадиоАстрону: %d)" % len(RUN))
A("   полностью потеряли отклик: %d  %s" % (len(rl), ["%s %s" % k for k in rl]))
A("   полностью приобрели отклик: %d  %s" % (len(rg), ["%s %s" % k for k in rg]))
A("   отклика нет ни в одной модели: %d  %s" % (len(r0), ["%s %s" % k for k in r0]))
A("")
A("СМЕНА СТАТУСА ПОСКАННО")
for lab, arr in (("ПРОПАЛ", lost), ("ПОЯВИЛСЯ", gain)):
    for r in arr:
        ro, rn = tot(RUN[RK(r)])
        A("   %-8s %-9s %s %s %-8s опорная %-3s/%-3s  snr %6.2f -> %6.2f   сеанс %d -> %d из %d"
          % (lab, r['exp'], r['band'], r['pol'], r['t'], r['ref_old'], r['ref_new'],
             r['snr_old'], r['snr_new'], ro, rn, len(RUN[RK(r)])))
A("")
A("СЕАНСЫ ЦЕЛИКОМ (скан за сканом; . = нет отклика, X = есть)")
for k in sorted(RUN):
    v = sorted(RUN[k], key=lambda r: (r['pol'], r['t']))
    A("   %-9s %s  было %s   стало %s   (%d -> %d)"
      % (k[0], k[1], "".join("X" if r['det_old'] else "." for r in v),
         "".join("X" if r['det_new'] else "." for r in v),
         tot(v)[0], tot(v)[1]))
open(os.path.join(H, "REPORT_ra.txt"), "w", encoding="utf-8").write("\n".join(L))
print("ra_scans.csv, REPORT_ra.txt: сканов %d, отклик %d -> %d" % (len(D), o, n))
