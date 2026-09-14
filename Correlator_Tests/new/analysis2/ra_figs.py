#!/usr/bin/env python
# ra_figs.py - графики по космическому плечу, счёт ПО СКАНАМ (см. ra_scan.py).
import csv, os
from collections import defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

H = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(H, "figures"); os.makedirs(OUT, exist_ok=True)
C_OLD, C_NEW = "tab:blue", "tab:red"
THR = 11.0

R = []
for r in csv.DictReader(open(os.path.join(H, "ra_scans.csv"), encoding="utf-8")):
    for k in ('snr_old','snr_new','rate_old','rate_new','delay_old','delay_new'): r[k] = float(r[k])
    for k in ('det_old','det_new','ref_same'): r[k] = int(r[k])
    for k in ('bproj_old','bproj_new'): r[k] = float(r[k]) if r[k] else None
    R.append(r)
D = [r for r in R if not (r['exp'] == 'raks08ad' and r['band'] == 'C'
                          and (r['ref_old'] == 'Ef' or r['ref_new'] == 'Ef'))]
RK = lambda r: (r['exp'], r['band'])
RUN = defaultdict(list)
for r in D: RUN[RK(r)].append(r)
tot = lambda a: (sum(x['det_old'] for x in a), sum(x['det_new'] for x in a))
lost = [r for r in D if r['det_old'] and not r['det_new']]
gain = [r for r in D if not r['det_old'] and r['det_new']]

# ------------------------------------------------- R1: карта сеансов скан за сканом
keys = sorted(RUN, key=lambda k: (tot(RUN[k])[1] - tot(RUN[k])[0], -tot(RUN[k])[0], k))
nmax = max(len(v) for v in RUN.values())
fig, ax = plt.subplots(figsize=(11.5, 9))
for i, k in enumerate(keys):
    v = sorted(RUN[k], key=lambda r: (r['pol'], r['t']))
    o, n = tot(v)
    for j, r in enumerate(v):
        for dy, det, c in ((-0.20, r['det_old'], C_OLD), (0.20, r['det_new'], C_NEW)):
            ax.add_patch(plt.Rectangle((j + .06, i + dy - .17), .88, .34,
                                       facecolor=(c if det else "white"),
                                       edgecolor=c, linewidth=1.1,
                                       hatch=("//" if (det and c == C_NEW) else None)))
    if n != o:
        ax.text(nmax + .45, i, "%+d" % (n - o), va="center", fontsize=10, weight="bold",
                color=("green" if n > o else "darkred"))
ax.set_xlim(0, nmax + 1.5); ax.set_ylim(-.6, len(keys) - .4)
ax.set_yticks(range(len(keys)))
ax.set_yticklabels(["%s %s" % k for k in keys], fontsize=8.5, family="monospace")
ax.set_xticks([j + .5 for j in range(nmax)])
ax.set_xticklabels(["скан %d" % (j + 1) for j in range(nmax)], fontsize=9)
ax.set_title("Отклик на РадиоАстроне по сканам каждого сеанса", fontsize=13, pad=12)
# пояснение и легенда - под осью, чтобы не наезжать на заголовок
ax.legend(handles=[Patch(facecolor=C_OLD, edgecolor=C_OLD, label="было"),
                   Patch(facecolor=C_NEW, edgecolor=C_NEW, hatch="//", label="стало"),
                   Patch(facecolor="white", edgecolor="0.4", label="отклика нет")],
          loc="upper center", bbox_to_anchor=(.5, -.045), ncol=3, fontsize=10, frameon=False,
          handlelength=1.6, columnspacing=2.2)
fig.text(.5, .012, "в каждой клетке верхняя половина - было, нижняя - стало; "
                   "заливка = отклик обнаружен",
         ha="center", fontsize=9.5, color="0.3")
ax.invert_yaxis()
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(rect=(0, .035, 1, 1))
fig.savefig(os.path.join(OUT, "R1_seansy.png"), dpi=140); plt.close(fig)

# ------------------------------------------------------------- R2: сводка + смены
fig, axs = plt.subplots(1, 2, figsize=(13, 5.6))
ax = axs[0]
o, n = tot(D)
x = np.arange(2)
ax.bar(x - .18, [o, len(D) - o], .36, color=C_OLD, label="было")
ax.bar(x + .18, [n, len(D) - n], .36, color=C_NEW, hatch="//", edgecolor="k",
       linewidth=.6, label="стало")
for xi, a, b in zip(x, [o, len(D) - o], [n, len(D) - n]):
    ax.text(xi - .18, a + .7, str(a), ha="center", fontsize=10)
    ax.text(xi + .18, b + .7, str(b), ha="center", fontsize=10)
ax.set_xticks(x); ax.set_xticklabels(["отклик обнаружен", "отклика нет"])
ax.set_ylabel("число сканов")
ax.set_title("Отклик на РадиоАстроне: %d сканов\nбыло %d, стало %d (%+d)\nПРОПАЛ %d, ПОЯВИЛСЯ %d"
             % (len(D), o, n, n - o, len(lost), len(gain)), fontsize=10)
ax.legend(fontsize=9); ax.grid(axis="y", alpha=.3)

ax = axs[1]
ch = sorted(lost + gain, key=lambda r: (0 if r['det_old'] else 1, r['bproj_new'] or 0))
y = list(range(len(ch)))[::-1]
for yi, r in zip(y, ch):
    ax.plot([r['snr_old'], r['snr_new']], [yi, yi], color="0.5", lw=1, zorder=1)
    ax.scatter(r['snr_old'], yi, s=55, color=C_OLD, zorder=3)
    ax.scatter(r['snr_new'], yi, s=55, color=C_NEW, marker="s", zorder=3)
ax.axvline(THR, color="k", ls="--", lw=1.5)
ax.set_yticks(y)
ax.set_yticklabels(["%-9s %s %-5s оп.%-3s B=%4.1f  %d/%d" %
                    (r['exp'], r['band'], r['t'][:5], r['ref_new'], r['bproj_new'] or 0,
                     *tot(RUN[RK(r)])) for r in ch], fontsize=8, family="monospace")
for tl, r in zip(ax.get_yticklabels(), ch):
    tl.set_color("darkred" if r['det_old'] else "darkgreen")
ax.set_xlabel("SNR (условный)")
ax.set_title("Сканы, сменившие статус\nкрасным - пропал, зелёным - появился\n"
             "оп. = опорная антенна, B - проекция базы, N/M - отклики в сеансе", fontsize=10)
ax.grid(axis="x", alpha=.3)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "R2_svodka.png"), dpi=140); plt.close(fig)

# ------------------------------------------------ R3: зависимость от проекции базы
GR = []
for r in csv.DictReader(open(os.path.join(H, "compare.csv"), encoding="utf-8")):
    if r['sta'] == 'RA' or r['ref'] == 'RA': continue
    if r['exp'] == 'raks08ad' and r['band'] == 'C' and r['ref'] == 'Ef': continue
    if not r['bproj']: continue
    GR.append(dict(b=float(r['bproj']), o=int(r['det_old']), n=int(r['det_new'])))

fig, axs = plt.subplots(1, 2, figsize=(13, 5.6))
edges = [1.5, 3, 5, 7, 10]
lab, fo, fn = [], [], []
for i2 in range(len(edges) - 1):
    bo = [r for r in D if r['bproj_old'] is not None and edges[i2] <= r['bproj_old'] < edges[i2+1]]
    bn = [r for r in D if r['bproj_new'] is not None and edges[i2] <= r['bproj_new'] < edges[i2+1]]
    if not bo or not bn: continue
    lab.append("%.2g-%.2g\n(%d/%d)" % (edges[i2], edges[i2+1], len(bo), len(bn)))
    fo.append(100. * sum(r['det_old'] for r in bo) / len(bo))
    fn.append(100. * sum(r['det_new'] for r in bn) / len(bn))
ax = axs[0]
x = np.arange(len(lab))
ax.plot(x, fo, "-o", color=C_OLD, lw=2, label="было")
ax.plot(x, fn, "--s", color=C_NEW, lw=2, label="стало")
for xi, a, b in zip(x, fo, fn):
    ax.annotate("%+.0f%%" % (b - a), (xi, max(a, b) + 3), ha="center", fontsize=9,
                color=("green" if b >= a else "darkred"))
ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=9)
ax.set_xlabel("проекция базы РадиоАстрон - опорная, диаметры Земли\n"
              "(в скобках - сканов в старом / новом прогоне)")
ax.set_ylabel("доля сканов с откликом, %")
ax.set_ylim(0, 108); ax.legend(); ax.grid(alpha=.3)
ax.set_title("Космическое плечо (счёт по сканам)", fontsize=11)

ax = axs[1]
ge = [0, .2, .4, .6, .9]
lab, fo, fn = [], [], []
for i2 in range(len(ge) - 1):
    b = [r for r in GR if ge[i2] <= r['b'] < ge[i2+1]]
    if not b: continue
    lab.append("%.2g-%.2g\n(%d)" % (ge[i2], ge[i2+1], len(b)))
    fo.append(100. * sum(r['o'] for r in b) / len(b))
    fn.append(100. * sum(r['n'] for r in b) / len(b))
x = np.arange(len(lab))
ax.plot(x, fo, "-o", color=C_OLD, lw=2, label="было")
ax.plot(x, fn, "--s", color=C_NEW, lw=2, label="стало")
for xi, a, b in zip(x, fo, fn):
    ax.annotate("%+.0f%%" % (b - a), (xi, max(a, b) + 3), ha="center", fontsize=9,
                color=("green" if b >= a else "darkred"))
ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=9)
ax.set_xlabel("проекция базы, диаметры Земли\n(в скобках - число записей)")
ax.set_ylabel("доля записей с корреляцией, %")
ax.set_ylim(0, 108); ax.legend(); ax.grid(alpha=.3)
ax.set_title("Наземные базы (счёт по записям база x скан)", fontsize=11)
fig.suptitle("Доля обнаруженных откликов в зависимости от длины проекции базы", fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "R3_proekciya.png"), dpi=140); plt.close(fig)

print("R1_seansy.png, R2_svodka.png, R3_proekciya.png")
print("сканов %d, отклик %d -> %d, пропал %d, появился %d" % (len(D), o, n, len(lost), len(gain)))
for i in range(len(edges) - 1):
    bo = [r for r in D if r['bproj_old'] is not None and edges[i] <= r['bproj_old'] < edges[i+1]]
    bn = [r for r in D if r['bproj_new'] is not None and edges[i] <= r['bproj_new'] < edges[i+1]]
    if bo and bn:
        print("  B %g-%g: было %d/%d = %.1f%%, стало %d/%d = %.1f%%" % (
            edges[i], edges[i+1], sum(r['det_old'] for r in bo), len(bo),
            100.*sum(r['det_old'] for r in bo)/len(bo),
            sum(r['det_new'] for r in bn), len(bn),
            100.*sum(r['det_new'] for r in bn)/len(bn)))
