#!/usr/bin/env python
# figs.py - графики для отчёта. Прямое сравнение старого прогона (без _p) и нового (_p).
# Везде: было = сплошная синяя, стало = красный пунктир.
import csv, os, statistics as st
from collections import defaultdict, Counter
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

H = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(H, "figures"); os.makedirs(OUT, exist_ok=True)
C_OLD, C_NEW = "tab:blue", "tab:red"
THR = 11.0

R = []
for r in csv.DictReader(open(os.path.join(H, "compare.csv"), encoding="utf-8")):
    r['det_old'] = int(r['det_old']); r['det_new'] = int(r['det_new'])
    for k in ('snr_old','snr_new','delay_old','delay_new','rate_old','rate_new'): r[k] = float(r[k])
    r['bproj'] = float(r['bproj']) if r['bproj'] else None
    R.append(r)
BUG = lambda r: r['exp'] == 'raks08ad' and r['band'] == 'C' and r['ref'] == 'Ef'
bug = [r for r in R if BUG(r)]
D   = [r for r in R if not BUG(r)]
kept = [r for r in D if r['det_old'] and r['det_new']]
lost = [r for r in D if r['det_old'] and not r['det_new']]
gain = [r for r in D if not r['det_old'] and r['det_new']]
SPACE = lambda r: r['sta'] == 'RA' or r['ref'] == 'RA'
BK = lambda r: (r['exp'], r['band'], r['pol'], r['ref'], r['sta'])
RK = lambda r: (r['exp'], r['band'])
RUN, BASE = defaultdict(list), defaultdict(list)
for r in D: RUN[RK(r)].append(r); BASE[BK(r)].append(r)
tot = lambda a: (sum(x['det_old'] for x in a), sum(x['det_new'] for x in a))
def verdict(r):
    ro, rn = tot(RUN[RK(r)]); bo, bn = tot(BASE[BK(r)])
    if rn == 0 and ro > 0: return "весь прогон потерял"
    if bo == 0 and bn > 0: return "вся база приобрела"
    if bn == 0 and bo > 0: return "вся база потеряла"
    return "только этот интервал"

BINS = [(0,.25,"0-0.25"),(.25,.5,"0.25-0.5"),(.5,1,"0.5-1"),(1,3,"1-3"),(3,6,"3-6"),(6,99,"6-10")]

# ============ 1. СВОДКА ============
bo_, bn_ = sum(r['det_old'] for r in D), sum(r['det_new'] for r in D)
def cnt(f, g): return sum(1 for r in kept if g(abs(f(r,'new'))) ) if False else 0
def direct(fo, fn, absv=True):
    up = dn = eq = 0
    for r in kept:
        a, b = (abs(fo(r)), abs(fn(r))) if absv else (fo(r), fn(r))
        up += b > a; dn += b < a; eq += b == a
    return up, dn, eq
s_up, s_dn, s_eq = direct(lambda r: r['snr_old'], lambda r: r['snr_new'], False)
d_up, d_dn, d_eq = direct(lambda r: r['delay_old'], lambda r: r['delay_new'])
r_up, r_dn, r_eq = direct(lambda r: r['rate_old'], lambda r: r['rate_new'])
runs_lost = sum(1 for k, v in RUN.items() if tot(v)[0] > 0 and tot(v)[1] == 0)
runs_gain = sum(1 for k, v in RUN.items() if tot(v)[0] == 0 and tot(v)[1] > 0)
bas_lost  = sum(1 for k, v in BASE.items() if tot(v)[0] > 0 and tot(v)[1] == 0)
bas_gain  = sum(1 for k, v in BASE.items() if tot(v)[0] == 0 and tot(v)[1] > 0)
rows = [
 ["ВЫБОРКА", "", "", ""],
 ["записей (база x момент)", "%d" % len(D), "", "исключено %d (raks08ad C, затёртый полином)" % len(bug)],
 ["  наземные / с РадиоАстроном", "%d / %d" % (sum(1 for r in D if not SPACE(r)), sum(1 for r in D if SPACE(r))), "", "проекция базы 0.02 - 9.76 диаметров Земли"],
 ["КОРРЕЛЯЦИЯ ЕСТЬ / НЕТ", "было", "стало", ""],
 ["записей с корреляцией", "%d" % bo_, "%d" % bn_, "изменение %+d" % (bn_ - bo_)],
 ["  есть в обеих", "%d" % len(kept), "", ""],
 ["  ПРОПАЛА", "%d" % len(lost), "", "из них: вся база 2, один интервал 4"],
 ["  ПОЯВИЛАСЬ", "", "%d" % len(gain), "из них: вся база 7, один интервал 3"],
 ["ОХВАТ ИЗМЕНЕНИЙ", "", "", ""],
 ["прогонов потеряло корреляцию ЦЕЛИКОМ", "", "%d" % runs_lost, "ни один прогон не остался без корреляции"],
 ["прогонов ПРИОБРЕЛО корреляцию целиком", "", "%d" % runs_gain, "raes03fr K и raks01ad K: было 0 корреляций, стало 2"],
 ["баз потеряло корреляцию целиком", "", "%d" % bas_lost, "обе базы имели 1-2 интервала"],
 ["баз ПРИОБРЕЛО корреляцию целиком", "", "%d" % bas_gain, "5 из 6 - базы с РадиоАстроном"],
 ["ВЕЛИЧИНЫ (прямой счёт по %d записям)" % len(kept), "стало больше", "стало меньше", "без изменений"],
 ["SNR", "%d" % s_up, "%d" % s_dn, "%d" % s_eq],
 ["|остаточная задержка|", "%d" % d_up, "%d" % d_dn, "%d" % d_eq],
 ["|остаточный fringe rate|", "%d" % r_up, "%d" % r_dn, "%d  <- единственный явный сдвиг" % r_eq],
]
fig, ax = plt.subplots(figsize=(15.5, 8.2)); ax.axis('off')
t = ax.table(cellText=rows, colLabels=["показатель", "было", "стало", "комментарий"],
             cellLoc='left', colWidths=[.32, .12, .12, .44], loc='center')
t.auto_set_font_size(False); t.set_fontsize(10); t.scale(1, 1.6)
for (i, j), c in t.get_celld().items():
    c.set_edgecolor('lightgray')
    if i == 0: c.set_facecolor('#dfe6ef'); c.set_text_props(weight='bold')
    elif rows[i-1][1] in ("", "было", "стало больше"): c.set_facecolor('#eef2f7'); c.set_text_props(weight='bold')
ax.set_title("СВОДКА. Старая модель задержки (БЫЛО) -> новая (СТАЛО). 15 сеансов, 30 прогонов коррелятора",
             fontsize=14, weight='bold', pad=24)
fig.text(.5, .02, "Корреляция = строка CLOCK в выходном файле не закомментирована '#'. "
                  "SNR в программе считается некорректно (сигнал/шум + шум, пол ~10), отклик признаётся примерно с 11.",
         ha='center', fontsize=9, style='italic')
fig.tight_layout(rect=[0, .035, 1, .95]); fig.savefig(os.path.join(OUT, "1_svodka.png"), dpi=115); plt.close(fig)

# ============ 2. ЗАВИСИМОСТЬ ОТ ПРОЕКЦИИ БАЗЫ ============
fig, ax = plt.subplots(1, 2, figsize=(15, 6))
labs, fo, fn, ns = [], [], [], []
for lo, hi, lab in BINS:
    s = [r for r in D if r['bproj'] is not None and lo <= r['bproj'] < hi]
    if not s: continue
    labs.append(lab); ns.append(len(s))
    fo.append(100 * sum(r['det_old'] for r in s) / len(s)); fn.append(100 * sum(r['det_new'] for r in s) / len(s))
x = np.arange(len(labs)); w = .38
a = ax[0]
a.bar(x - w/2, fo, w, color=C_OLD, alpha=.85, label="было")
a.bar(x + w/2, fn, w, color=C_NEW, alpha=.85, label="стало")
for i, (o, n, m) in enumerate(zip(fo, fn, ns)):
    a.text(i - w/2, o, "%.0f%%" % o, ha='center', va='bottom', fontsize=9)
    a.text(i + w/2, n, "%.0f%%" % n, ha='center', va='bottom', fontsize=9)
    a.text(i, -6, "n=%d" % m, ha='center', fontsize=8, color='dimgray')
a.set_xticks(x); a.set_xticklabels(labs); a.set_ylim(0, 100)
a.set_xlabel("проекция базы, диаметры Земли"); a.set_ylabel("доля записей с корреляцией, %")
a.set_title("(а) Доля найденных откликов по длине базы", fontsize=12)
a.legend(fontsize=10); a.grid(axis='y', alpha=.3)
a.axvline(2.5, color='gray', ls=':', lw=1.2)
a.text(2.55, 95, "правее — базы с РадиоАстроном", fontsize=9, color='dimgray')
a = ax[1]
a.scatter([r['bproj'] for r in kept], [r['snr_new']/r['snr_old'] for r in kept], s=16,
          c='lightsteelblue', label="корреляция в обеих", edgecolors='none')
a.scatter([r['bproj'] for r in lost], [r['snr_new']/r['snr_old'] for r in lost], s=90, marker='v',
          c='tab:red', edgecolors='k', label="ПРОПАЛА", zorder=4)
a.scatter([r['bproj'] for r in gain], [r['snr_new']/r['snr_old'] for r in gain], s=90, marker='^',
          c='tab:cyan', edgecolors='k', label="ПОЯВИЛАСЬ", zorder=4)
a.axhline(1, color='k', lw=1.2); a.set_xscale('log'); a.set_yscale('log')
a.set_xlabel("проекция базы, диаметры Земли"); a.set_ylabel("SNR стало / было")
a.set_title("(б) Изменение SNR и смены статуса по длине базы", fontsize=12)
a.legend(fontsize=9); a.grid(alpha=.3, which='both')
fig.suptitle("ЗАВИСИМОСТЬ ОТ ПРОЕКЦИИ БАЗЫ", fontsize=15, weight='bold')
fig.tight_layout(rect=[0, 0, 1, .93]); fig.savefig(os.path.join(OUT, "2_proekciya_bazy.png"), dpi=115); plt.close(fig)

# ============ 3. ОХВАТ ИЗМЕНЕНИЙ ============
fig, ax = plt.subplots(1, 2, figsize=(15, 5.6))
cl = Counter(verdict(r) for r in lost); cg = Counter(verdict(r) for r in gain)
order = ["только этот интервал", "вся база потеряла", "вся база приобрела", "весь прогон потерял"]
a = ax[0]
la = [k for k in order if cl.get(k) or cg.get(k)]
xx = np.arange(len(la))
a.bar(xx - .2, [cl.get(k, 0) for k in la], .4, color='tab:red', alpha=.85, label="ПРОПАЛА")
a.bar(xx + .2, [cg.get(k, 0) for k in la], .4, color='tab:cyan', alpha=.85, label="ПОЯВИЛАСЬ")
for i, k in enumerate(la):
    if cl.get(k): a.text(i - .2, cl[k], str(cl[k]), ha='center', va='bottom', fontsize=11, weight='bold')
    if cg.get(k): a.text(i + .2, cg[k], str(cg[k]), ha='center', va='bottom', fontsize=11, weight='bold')
a.set_xticks(xx); a.set_xticklabels([k.replace(" ", "\n", 1) for k in la], fontsize=9)
a.set_ylabel("число записей"); a.set_title("(а) Какого масштаба изменение", fontsize=12)
a.legend(fontsize=10); a.grid(axis='y', alpha=.3)
a = ax[1]
cats = ["прогоны\nпотеряли", "прогоны\nПРИОБРЕЛИ", "базы\nпотеряли", "базы\nПРИОБРЕЛИ"]
vals = [runs_lost, runs_gain, bas_lost, bas_gain]
cols = ['tab:red', 'tab:cyan', 'tab:red', 'tab:cyan']
b_ = a.bar(cats, vals, color=cols, alpha=.85)
for r_, v_ in zip(b_, vals): a.text(r_.get_x()+r_.get_width()/2, v_, str(v_), ha='center', va='bottom', fontsize=13, weight='bold')
a.set_ylabel("количество"); a.set_title("(б) Целиком потеряно / приобретено", fontsize=12); a.grid(axis='y', alpha=.3)
a.text(.5, .78, "ни один прогон не остался без корреляции;\nдва прогона (raes03fr K, raks01ad K) её ПРИОБРЕЛИ",
       transform=a.transAxes, ha='center', fontsize=10, bbox=dict(fc='lightyellow', ec='gray'))
fig.suptitle("ОХВАТ ИЗМЕНЕНИЙ: весь прогон / вся база / один интервал", fontsize=15, weight='bold')
fig.tight_layout(rect=[0, 0, 1, .92]); fig.savefig(os.path.join(OUT, "3_ohvat.png"), dpi=115); plt.close(fig)

# ============ 4/5. ДЕТАЛИ КАЖДОГО СЛУЧАЯ ============
def cases_fig(cases, title, fname, mark):
    n = len(cases); nc = 2; nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(14, 3.1*nr), squeeze=False)
    for i, r in enumerate(sorted(cases, key=lambda r: (r['exp'], r['band'], r['t']))):
        a = axs[i//nc][i%nc]
        ser = sorted(BASE[BK(r)], key=lambda z: z['t'])
        xs = [z['t'][:5] for z in ser]
        a.plot(range(len(ser)), [z['snr_old'] for z in ser], '-o', color=C_OLD, lw=1.6, ms=6, label="было")
        a.plot(range(len(ser)), [z['snr_new'] for z in ser], '--s', color=C_NEW, lw=1.6, ms=6, label="стало")
        for j, z in enumerate(ser):
            if not z['det_old']: a.plot(j, z['snr_old'], 'o', mfc='white', mec=C_OLD, ms=9, mew=1.8, zorder=5)
            if not z['det_new']: a.plot(j, z['snr_new'], 's', mfc='white', mec=C_NEW, ms=9, mew=1.8, zorder=5)
        a.axhline(THR, color='k', ls=':', lw=1.2)
        a.axvline([z['t'] for z in ser].index(r['t']), color=mark, lw=2.4, alpha=.3)
        a.set_xticks(range(len(ser))); a.set_xticklabels(xs, fontsize=8)
        ro, rn = tot(RUN[RK(r)]); bo, bn = tot(BASE[BK(r)])
        a.set_title("%s %s   база %s-%s   B=%.2f диам.Земли\nпрогон: корреляций %d->%d   база: %d->%d   [%s]" % (
            r['exp'], r['band'], r['sta'], r['ref'], r['bproj'], ro, rn, bo, bn, verdict(r)), fontsize=9)
        a.set_ylabel("SNR", fontsize=9); a.grid(alpha=.3)
        if i == 0: a.legend(fontsize=8)
    for i in range(n, nr*nc): axs[i//nc][i%nc].axis('off')
    fig.suptitle(title, fontsize=14, weight='bold')
    fig.tight_layout(rect=[0, 0, 1, 1 - .32/nr]); fig.savefig(os.path.join(OUT, fname), dpi=110); plt.close(fig)

cases_fig(lost, "ГДЕ КОРРЕЛЯЦИЯ ПРОПАЛА (%d записей). Полый маркер = корреляции нет, пунктир по горизонтали = порог 11.\n"
                "Показана вся серия базы, чтобы видеть: пропало на одном интервале или на всей базе." % len(lost),
          "4_propala.png", "tab:red")
cases_fig(gain, "ГДЕ КОРРЕЛЯЦИЯ ПОЯВИЛАСЬ (%d записей). Полый маркер = корреляции нет.\n"
                "Показана вся серия базы: видно, приобрела вся база или только один интервал." % len(gain),
          "5_poyavilas.png", "tab:cyan")

# ============ 6. ВЕЛИЧИНЫ: прямое сравнение ============
def pair_panel(a, xs, ys, xlab, ylab, title, note, log=True):
    g = [(x, y) for x, y, r in zip(xs, ys, kept) if not SPACE(r) and (x > 0 and y > 0 or not log)]
    s = [(x, y) for x, y, r in zip(xs, ys, kept) if SPACE(r) and (x > 0 and y > 0 or not log)]
    a.scatter([p[0] for p in g], [p[1] for p in g], s=20, alpha=.7, c=C_OLD, label="наземные", edgecolors='none')
    a.scatter([p[0] for p in s], [p[1] for p in s], s=20, alpha=.7, c=C_NEW, label="с РадиоАстроном", edgecolors='none')
    allv = [v for p in g + s for v in p]
    if allv:
        lim = [min(allv)*.6, max(allv)*1.7]; a.plot(lim, lim, 'k-', lw=1.2); a.set_xlim(lim); a.set_ylim(lim)
    if log: a.set_xscale('log'); a.set_yscale('log')
    a.set_xlabel(xlab); a.set_ylabel(ylab); a.set_title(title, fontsize=11); a.grid(alpha=.3, which='both')
    a.legend(fontsize=9, loc='upper left')
    a.text(.97, .04, note, transform=a.transAxes, ha='right', fontsize=9,
           bbox=dict(fc='lightyellow', ec='gray', alpha=.92))

fig, ax = plt.subplots(1, 3, figsize=(17.5, 5.8))
pair_panel(ax[0], [r['snr_old'] for r in kept], [r['snr_new'] for r in kept],
           "SNR было", "SNR стало", "(а) SNR",
           "выросло %d, упало %d\nиз %d записей" % (s_up, s_dn, len(kept)))
pair_panel(ax[1], [abs(r['rate_old']) for r in kept], [abs(r['rate_new']) for r in kept],
           "|fringe rate| было", "|fringe rate| стало", "(б) Остаточный fringe rate",
           "МЕНЬШЕ стало у %d,\nбольше у %d  <- главный выигрыш" % (r_dn, r_up))
pair_panel(ax[2], [abs(r['delay_old'])*1e6 for r in kept], [abs(r['delay_new'])*1e6 for r in kept],
           "|задержка| было, мкс", "|задержка| стало, мкс", "(в) Остаточная задержка",
           "меньше %d, больше %d,\nбез изменений %d" % (d_dn, d_up, d_eq))
fig.suptitle("ВЕЛИЧИНЫ: прямое сравнение по %d записям, где корреляция есть в обеих моделях" % len(kept),
             fontsize=15, weight='bold')
fig.tight_layout(rect=[0, 0, 1, .93]); fig.savefig(os.path.join(OUT, "6_velichiny.png"), dpi=115); plt.close(fig)

print("записей %d | корреляция %d -> %d | пропала %d | появилась %d" % (len(D), bo_, bn_, len(lost), len(gain)))
print("SNR: +%d/-%d | rate: меньше %d/больше %d | задержка: меньше %d/больше %d" % (s_up, s_dn, r_dn, r_up, d_dn, d_up))
for f in sorted(os.listdir(OUT)): print("  figures/" + f)
