#!/usr/bin/env python
# figs2.py - графики с РАЗДЕЛЕНИЕМ на наземно-космические (с РадиоАстроном) и наземные базы.
# Везде: было = сплошная синяя, стало = красный пунктир/квадрат.
import csv, os
from collections import defaultdict
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
BUG   = lambda r: r['exp'] == 'raks08ad' and r['band'] == 'C' and r['ref'] == 'Ef'
SPACE = lambda r: r['sta'] == 'RA' or r['ref'] == 'RA'
D  = [r for r in R if not BUG(r)]
SP = [r for r in D if SPACE(r)]
GR = [r for r in D if not SPACE(r)]
BK = lambda r: (r['exp'], r['band'], r['pol'], r['ref'], r['sta'])
RK = lambda r: (r['exp'], r['band'])
tot = lambda a: (sum(x['det_old'] for x in a), sum(x['det_new'] for x in a))
BASE, RUN = defaultdict(list), defaultdict(list)
for r in D: BASE[BK(r)].append(r); RUN[RK(r)].append(r)

def split3(S):
    return ([r for r in S if r['det_old'] and r['det_new']],
            [r for r in S if r['det_old'] and not r['det_new']],
            [r for r in S if not r['det_old'] and r['det_new']])

def scope(r):
    bo, bn = tot(BASE[BK(r)])
    if bn == 0 and bo > 0: return "вся база"
    if bo == 0 and bn > 0: return "вся база"
    return "один скан"

SETS = [("Наземно-космические базы (РадиоАстрон)", SP), ("Наземные базы", GR)]

# ------------------------------------------------------------------ 1. Сводка
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.4))
for ax, (name, S) in zip(axs, SETS):
    kept, lost, gain = split3(S)
    o, n = tot(S)
    x = np.arange(2)
    ax.bar(x - 0.18, [o, len(S) - o], 0.36, color=C_OLD, label="было (старая модель)")
    ax.bar(x + 0.18, [n, len(S) - n], 0.36, color=C_NEW, hatch="//", edgecolor="k",
           linewidth=.6, label="стало (новая модель)")
    for xi, a, b in zip(x, [o, len(S) - o], [n, len(S) - n]):
        ax.text(xi - 0.18, a + 1, str(a), ha="center", fontsize=10)
        ax.text(xi + 0.18, b + 1, str(b), ha="center", fontsize=10)
    ax.set_xticks(x); ax.set_xticklabels(["есть корреляция", "нет корреляции"])
    ax.set_ylabel("число записей (база x скан)")
    ax.set_title("%s\nвсего %d записей, %+d корреляций\nПРОПАЛА %d, ПОЯВИЛАСЬ %d"
                 % (name, len(S), n - o, len(lost), len(gain)), fontsize=10)
    ax.legend(fontsize=8); ax.grid(axis="y", alpha=.3)
fig.suptitle("Корреляционные отклики: старая модель задержки против новой", fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "S1_svodka.png"), dpi=140); plt.close(fig)

# ------------------------------------------- 2. Зависимость от проекции базы
fig, axs = plt.subplots(1, 2, figsize=(13, 5.4))
BINS = {"Наземно-космические базы (РадиоАстрон)": [1.5, 3, 5, 7, 10],
        "Наземные базы": [0, .2, .4, .6, .9]}
for ax, (name, S) in zip(axs, SETS):
    edges = BINS[name]
    lab, fo, fn = [], [], []
    for i in range(len(edges) - 1):
        b = [r for r in S if r['bproj'] is not None and edges[i] <= r['bproj'] < edges[i + 1]]
        if not b: continue
        o, n = tot(b)
        lab.append("%.2g-%.2g\n(%d)" % (edges[i], edges[i + 1], len(b)))
        fo.append(100. * o / len(b)); fn.append(100. * n / len(b))
    x = np.arange(len(lab))
    ax.plot(x, fo, "-o", color=C_OLD, lw=2, label="было")
    ax.plot(x, fn, "--s", color=C_NEW, lw=2, label="стало")
    for xi, a, b in zip(x, fo, fn):
        ax.annotate("%+.0f%%" % (b - a), (xi, max(a, b) + 3), ha="center", fontsize=9,
                    color=("green" if b >= a else "darkred"))
    ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=9)
    ax.set_xlabel("проекция базы, диаметры Земли (в скобках - число записей)")
    ax.set_ylabel("доля записей с корреляцией, %")
    ax.set_ylim(0, 108); ax.set_title(name, fontsize=11)
    ax.legend(); ax.grid(alpha=.3)
fig.suptitle("Доля обнаруженных откликов в зависимости от длины проекции базы", fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "S2_proekciya.png"), dpi=140); plt.close(fig)

# ------------------------------------------------ 3. SNR было-стало, все точки
fig, axs = plt.subplots(1, 2, figsize=(13, 6))
for ax, (name, S) in zip(axs, SETS):
    kept, lost, gain = split3(S)
    none = [r for r in S if not r['det_old'] and not r['det_new']]
    ax.scatter([r['snr_old'] for r in kept], [r['snr_new'] for r in kept], s=18,
               c="tab:green", alpha=.6, label="корреляция в обеих (%d)" % len(kept))
    ax.scatter([r['snr_old'] for r in none], [r['snr_new'] for r in none], s=14,
               c="0.65", alpha=.6, label="нет ни в одной (%d)" % len(none))
    ax.scatter([r['snr_old'] for r in lost], [r['snr_new'] for r in lost], s=95,
               marker="v", c="darkred", edgecolors="k", zorder=5, label="ПРОПАЛА (%d)" % len(lost))
    ax.scatter([r['snr_old'] for r in gain], [r['snr_new'] for r in gain], s=95,
               marker="^", c="tab:orange", edgecolors="k", zorder=5, label="ПОЯВИЛАСЬ (%d)" % len(gain))
    m = max([r['snr_old'] for r in S] + [r['snr_new'] for r in S])
    ax.plot([8, m * 1.1], [8, m * 1.1], "k:", lw=1, label="без изменений")
    # подписать три самых сильных относительных изменения среди устоявших записей
    strong = sorted(kept, key=lambda r: abs(np.log(r['snr_new'] / r['snr_old'])), reverse=True)[:3]
    for r in strong:
        ax.annotate("%s %s %s-%s  %.0f%%" % (r['exp'], r['band'], r['sta'], r['ref'],
                                             100. * r['snr_new'] / r['snr_old']),
                    (r['snr_old'], r['snr_new']), textcoords="offset points", xytext=(8, -12),
                    fontsize=7.5, color="0.25",
                    arrowprops=dict(arrowstyle="-", lw=.6, color="0.5"))
    ax.axhline(THR, color=C_NEW, ls="--", lw=1); ax.axvline(THR, color=C_OLD, lw=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("SNR (условный), старая модель"); ax.set_ylabel("SNR (условный), новая модель")
    ax.set_title(name, fontsize=11); ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=.3, which="both")
fig.suptitle("Отклик поштучно. Линии - порог обнаружения 11 "
             "(SNR в программе завышен аддитивным шумом)", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "S3_snr_tochki.png"), dpi=140); plt.close(fig)

# ------------------------------------------------------- 4. Случаи смены статуса
ch = [r for r in D if r['det_old'] != r['det_new']]
ch.sort(key=lambda r: (0 if SPACE(r) else 1, r['bproj'] or 0))
fig, ax = plt.subplots(figsize=(12.5, 6.5))
y = list(range(len(ch)))[::-1]
for yi, r in zip(y, ch):
    ax.plot([r['snr_old'], r['snr_new']], [yi, yi], color="0.5", lw=1, zorder=1)
    ax.scatter(r['snr_old'], yi, s=55, color=C_OLD, zorder=3)
    ax.scatter(r['snr_new'], yi, s=55, color=C_NEW, marker="s", zorder=3)
ax.axvline(THR, color="k", ls="--", lw=1.5)
ax.set_yticks(y)
labs = ["%-9s %s %-2s-%-2s %s  B=%4.2f  %-9s [%s]" %
        (r['exp'], r['band'], r['sta'], r['ref'], r['t'][:5], r['bproj'] or 0,
         "ПРОПАЛА" if r['det_old'] else "ПОЯВИЛАСЬ", scope(r)) for r in ch]
ax.set_yticklabels(labs, fontsize=8, family="monospace")
for tl, r in zip(ax.get_yticklabels(), ch):
    tl.set_color("darkred" if r['det_old'] else "darkgreen")
n_sp = sum(1 for r in ch if SPACE(r))
if 0 < n_sp < len(ch):
    ax.axhline(y[n_sp] + .5, color="k", lw=1.2)
ax.set_xlim(9.1, 12.9)
ax.text(12.85, y[0] + .4, "наземно-космические", ha="right", fontsize=10, weight="bold")
if n_sp < len(ch):
    ax.text(12.85, y[n_sp] + .0, "наземные", ha="right", fontsize=10, weight="bold")
ax.text(THR + .05, 0.2, "порог 11", ha="left", fontsize=9)
ax.set_xlabel("SNR (условный)")
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=C_OLD, label="было"),
                   Line2D([], [], marker="s", ls="", color=C_NEW, label="стало")],
          loc="lower right", fontsize=9)
ax.set_title("Все случаи смены статуса корреляции\n"
             "B - проекция базы в диаметрах Земли; сверху - наземно-космические (РадиоАстрон)",
             fontsize=12)
ax.grid(axis="x", alpha=.3)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "S4_smena_statusa.png"), dpi=140); plt.close(fig)

# --------------------------------------------------------- 5. Остаточный fringe rate
fig, axs = plt.subplots(1, 2, figsize=(13, 5.6))
for ax, (name, S) in zip(axs, SETS):
    kept, _, _ = split3(S)
    ro = np.abs(np.array([r['rate_old'] for r in kept]))
    rn = np.abs(np.array([r['rate_new'] for r in kept]))
    floor = 1e-15
    ro = np.maximum(ro, floor); rn = np.maximum(rn, floor)
    ax.scatter(ro, rn, s=20, c="tab:purple" if S is SP else "tab:cyan", alpha=.7)
    lo, hi = floor * .5, max(ro.max(), rn.max()) * 3
    ax.plot([lo, hi], [lo, hi], "k:", lw=1)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    better = int((rn < ro).sum()); worse = int((rn > ro).sum())
    mo = float(np.median(np.abs([r['rate_old'] for r in kept])))
    mn = float(np.median(np.abs([r['rate_new'] for r in kept])))
    ax.axvline(mo, color=C_OLD, lw=1.2); ax.axhline(mn, color=C_NEW, ls="--", lw=1.2)
    nz = sum(1 for r in kept if r['rate_new'] == 0.0)
    if nz:
        ax.text(hi * .9, floor, "нижний ряд - ровно 0 (%d зап.)" % nz,
                ha="right", va="bottom", fontsize=8, color="0.35")
    ax.set_xlabel("|остаточный fringe rate|, старая модель")
    ax.set_ylabel("|остаточный fringe rate|, новая модель")
    ax.set_title("%s\nниже диагонали (стало лучше): %d,  выше: %d  из %d\n"
                 "медиана (линии): %.2e -> %.2e"
                 % (name, better, worse, len(kept), mo, mn), fontsize=10)
    ax.grid(alpha=.3, which="both")
fig.suptitle("Остаточная частота интерференции по записям с корреляцией в обеих моделях", fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "S5_fringe_rate.png"), dpi=140); plt.close(fig)

# ------------------------------------------------------ 6. Прямой счёт величин
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.2))
for ax, (name, S) in zip(axs, SETS):
    kept, _, _ = split3(S)
    # SNR намеренно не показан: у порога он определяется шумом, а не моделью
    labels = ["|fringe rate|", "|задержка|"]
    fns = [lambda r: (abs(r['rate_new']), abs(r['rate_old'])),
           lambda r: (abs(r['delay_new']), abs(r['delay_old']))]
    up = [sum(1 for r in kept if f(r)[0] > f(r)[1]) for f in fns]
    dn = [sum(1 for r in kept if f(r)[0] < f(r)[1]) for f in fns]
    x = np.arange(len(labels))
    ax.bar(x - .18, up, .36, color="tab:orange", label="стало больше")
    ax.bar(x + .18, dn, .36, color="tab:green", label="стало меньше")
    for xi, a, b in zip(x, up, dn):
        ax.text(xi - .18, a + .8, str(a), ha="center", fontsize=9)
        ax.text(xi + .18, b + .8, str(b), ha="center", fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("число записей"); ax.set_title("%s (%d записей)" % (name, len(kept)), fontsize=11)
    ax.set_ylim(0, max(up + dn) * 1.3)
    ax.legend(fontsize=8, loc="upper center", ncol=2); ax.grid(axis="y", alpha=.3)
fig.suptitle("Прямой счёт: у скольких записей остаточная величина выросла, у скольких упала\n"
             "меньше = лучше", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "S6_velichiny.png"), dpi=140); plt.close(fig)

print("готово:", ", ".join(sorted(f for f in os.listdir(OUT) if f.startswith("S"))))

# ---------------------------------- 7. По сеансам (эксперимент+диапазон), только космос
runs = sorted({RK(r) for r in SP})
vals = [(k, tot([r for r in SP if RK(r) == k])) for k in runs]
vals.sort(key=lambda kv: (kv[1][1] - kv[1][0], kv[1][0]))
fig, ax = plt.subplots(figsize=(11, 9))
y = np.arange(len(vals))
ax.barh(y + .2, [v[1][0] for v in vals], .38, color=C_OLD, label="было")
ax.barh(y - .2, [v[1][1] for v in vals], .38, color=C_NEW, hatch="//", edgecolor="k",
        linewidth=.5, label="стало")
ax.set_yticks(y)
ax.set_yticklabels(["%s %s" % k for k, _ in vals], fontsize=8, family="monospace")
for yi, (k, (o, n)) in zip(y, vals):
    if n != o:
        ax.text(max(o, n) + .15, yi, "%+d" % (n - o), va="center", fontsize=8,
                color=("green" if n > o else "darkred"), weight="bold")
ax.set_xlabel("число сканов с корреляцией на наземно-космических базах")
ax.set_title("Наземно-космические базы: корреляция по сеансам (эксперимент + диапазон)", fontsize=12)
ax.legend(loc="lower right"); ax.grid(axis="x", alpha=.3)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "S7_po_seansam.png"), dpi=140); plt.close(fig)
print("+ S7_po_seansam.png")
