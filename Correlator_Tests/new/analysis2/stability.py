#!/usr/bin/env python
# stability.py - устойчивость остаточного fringe rate от скана к скану внутри сеанса.
#
# Идея: на одной базе внутри одного сеанса остаточный fringe rate должен быть примерно
# постоянным (он определяется уходом станционных часов). Насколько он "гуляет" от скана
# к скану - это и есть мера того, насколько хорошо модель задержки описывает изменение
# геометрии во времени. Сравниваем РАЗМАХ (max - min) по одной и той же базе в старой
# и новой модели. Сравнение парное: берутся только те базы, где корреляция есть в обеих
# моделях, и размах считается по одним и тем же сканам.
#
# Выход: REPORT_stability.txt, figures/S8_stabilnost.png, figures/S9_primery.png

import csv, os, statistics as st
from collections import defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

H = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(H, "figures"); os.makedirs(OUT, exist_ok=True)
C_OLD, C_NEW = "tab:blue", "tab:red"

R = []
for r in csv.DictReader(open(os.path.join(H, "compare.csv"), encoding="utf-8")):
    for k in ('snr_old','snr_new','rate_old','rate_new'): r[k] = float(r[k])
    R.append(r)
BUG = lambda r: r['exp'] == 'raks08ad' and r['band'] == 'C' and r['ref'] == 'Ef'
SPk = lambda k: k[3] == 'RA' or k[4] == 'RA'

# только записи с корреляцией В ОБЕИХ моделях - иначе сравнение непарное
D = [r for r in R if not BUG(r) and r['det_old'] == '1' and r['det_new'] == '1']
B = defaultdict(list)
for r in D:
    B[(r['exp'], r['band'], r['pol'], r['ref'], r['sta'])].append(r)

rows = []
for k, v in B.items():
    if len(v) < 2: continue                      # размах требует минимум двух сканов
    v.sort(key=lambda r: r['t'])
    ro = [r['rate_old'] for r in v]; rn = [r['rate_new'] for r in v]
    rows.append(dict(k=k, v=v, n=len(v), sp=SPk(k),
                     po=max(ro) - min(ro), pn=max(rn) - min(rn),
                     so=st.pstdev(ro), sn=st.pstdev(rn)))

SETS = [("Наземно-космические базы (РадиоАстрон)", [x for x in rows if x['sp']]),
        ("Наземные базы", [x for x in rows if not x['sp']])]

# ------------------------------------------------------------------- текстовая сводка
L = []; A = L.append
A("=" * 96)
A("УСТОЙЧИВОСТЬ ОСТАТОЧНОГО FRINGE RATE ОТ СКАНА К СКАНУ")
A("Размах = (max - min) остаточного fringe rate по сканам одной базы внутри одного сеанса.")
A("Сравнение парное: только базы, где корреляция есть в обеих моделях, и не менее 2 сканов.")
A("=" * 96)
for name, S in SETS:
    better = sum(1 for x in S if x['pn'] < x['po'])
    worse  = sum(1 for x in S if x['pn'] > x['po'])
    same   = len(S) - better - worse
    rat = sorted(x['pn'] / x['po'] for x in S if x['po'] > 0)
    A("")
    A("### %s ###" % name)
    A("  баз: %d   сканов в них: %d" % (len(S), sum(x['n'] for x in S)))
    A("  стало стабильнее: %d | стало хуже: %d | без изменений: %d" % (better, worse, same))
    A("  медиана размаха: %.3g -> %.3g   (в %.2f раза)"
      % (st.median([x['po'] for x in S]), st.median([x['pn'] for x in S]),
         st.median([x['po'] for x in S]) / max(st.median([x['pn'] for x in S]), 1e-30)))
    A("  медиана СКО:     %.3g -> %.3g"
      % (st.median([x['so'] for x in S]), st.median([x['sn'] for x in S])))
    if rat:
        A("  отношение размаха новая/старая: медиана %.3f" % rat[len(rat) // 2])
        A("     стал вдвое и более узким (<0.5): %d баз" % sum(1 for x in rat if x < 0.5))
        A("     стал вдвое и более широким (>2): %d баз" % sum(1 for x in rat if x > 2))
S4 = sorted([x for x in rows if x['n'] >= 4 and x['po'] > 0], key=lambda x: x['pn'] / x['po'])
A("")
A("### Базы с 4 и более сканами, по возрастанию отношения размаха ###")
for x in S4:
    A("  %-4s %-9s %s %-3s-%-3s n=%d   %.3g -> %.3g   x%.2f"
      % ("КОСМ" if x['sp'] else "земл", x['k'][0], x['k'][1], x['k'][4], x['k'][3],
         x['n'], x['po'], x['pn'], x['pn'] / x['po']))
open(os.path.join(H, "REPORT_stability.txt"), "w", encoding="utf-8").write("\n".join(L))

# --------------------------------------------------- S8: статистика по всем базам
fig, axs = plt.subplots(2, 2, figsize=(13, 10))
FLOOR = 1e-15
for col, (name, S) in enumerate(SETS):
    # верхний ряд - парная диаграмма размаха
    ax = axs[0][col]
    po = np.maximum([x['po'] for x in S], FLOOR)
    pn = np.maximum([x['pn'] for x in S], FLOOR)
    ax.scatter(po, pn, s=np.array([x['n'] for x in S]) * 14,
               c=("tab:purple" if col == 0 else "tab:cyan"), alpha=.75, edgecolors="k", linewidths=.4)
    lo, hi = FLOOR * .5, max(po.max(), pn.max()) * 3
    ax.plot([lo, hi], [lo, hi], "k:", lw=1)
    ax.fill_between([lo, hi], [lo, lo], [lo, hi], color="green", alpha=.05)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    better = sum(1 for x in S if x['pn'] < x['po']); worse = sum(1 for x in S if x['pn'] > x['po'])
    ax.text(.04, .93, "ниже диагонали =\nстало стабильнее", transform=ax.transAxes,
            fontsize=9, color="darkgreen", va="top")
    ax.set_xlabel("размах fringe rate, старая модель")
    ax.set_ylabel("размах fringe rate, новая модель")
    ax.set_title("%s\nбаз %d: стабильнее %d, хуже %d\n(размер точки - число сканов)"
                 % (name, len(S), better, worse), fontsize=10)
    ax.grid(alpha=.3, which="both")

    # нижний ряд - распределение отношения
    ax = axs[1][col]
    rat = np.array([x['pn'] / x['po'] for x in S if x['po'] > 0])
    rat = np.clip(rat, 2.0 ** -6, 2.0 ** 6)
    bins = 2.0 ** np.arange(-6, 6.5, .5)
    ax.hist(rat, bins=bins, color=("tab:purple" if col == 0 else "tab:cyan"),
            edgecolor="k", linewidth=.5)
    ax.axvline(1, color="k", lw=1.5, label="без изменений")
    med = float(np.median(rat))
    ax.axvline(med, color=C_NEW, ls="--", lw=2, label="медиана %.2f" % med)
    ax.set_xscale("log", base=2); ax.set_xlim(2.0 ** -6.3, 2.0 ** 6.3)
    ax.set_xticks([2.0 ** k for k in range(-6, 7, 2)])
    ax.set_xticklabels(["1/64", "1/16", "1/4", "1", "4", "16", "64"])
    ax.set_xlabel("во сколько раз изменился размах (новая / старая)")
    ax.set_ylabel("число баз")
    ax.set_title("слева от единицы - стало стабильнее", fontsize=10)
    ax.legend(fontsize=9); ax.grid(axis="y", alpha=.3)
fig.suptitle("Устойчивость остаточного fringe rate от скана к скану\n"
             "размах = (max - min) по сканам одной базы внутри сеанса", fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "S8_stabilnost.png"), dpi=140); plt.close(fig)

# --------------------------------------------------- S9: четыре наглядных примера
PICK = [(('raks01ho', 'C', 'LL', 'Wb', 'RA'), "наземно-космическая: размах упал в 8 раз"),
        (('raks01ab', 'C', 'LL', 'Ev', 'RA'), "наземно-космическая: без изменений"),
        (('raks01ho', 'L', 'RR', 'Ef', 'Mc'), "наземная: систематический уход убран полностью"),
        (('raes03fr', 'C', 'LL', 'Wb', 'Bd'), "наземная: проигрыш (но уровень крайне мал)")]
by_key = {x['k']: x for x in rows}
fig, axs = plt.subplots(2, 2, figsize=(13, 8.5))
for ax, (key, caption) in zip(axs.ravel(), PICK):
    x = by_key.get(key)
    if x is None:
        ax.set_visible(False); continue
    t = [r['t'].rsplit(':', 1)[0] for r in x['v']]
    ii = np.arange(len(t))
    yo = [r['rate_old'] for r in x['v']]; yn = [r['rate_new'] for r in x['v']]
    ax.plot(ii, yo, "-o", color=C_OLD, lw=2, label="было")
    ax.plot(ii, yn, "--s", color=C_NEW, lw=2, label="стало")
    if min(yo + yn) < 0 < max(yo + yn):
        ax.axhline(0, color="0.7", lw=.8)
    ax.set_xticks(ii); ax.set_xticklabels(t, fontsize=9)
    ax.set_xlabel("скан"); ax.set_ylabel("остаточный fringe rate")
    ax.set_title("%s %s, база %s-%s\n%s\nразмах %.3g -> %.3g"
                 % (key[0], key[1], key[4], key[3], caption, x['po'], x['pn']), fontsize=10)
    ax.legend(fontsize=9); ax.grid(alpha=.3)
fig.suptitle("Что такое размах: остаточный fringe rate по сканам на четырёх конкретных базах",
             fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "S9_primery.png"), dpi=140); plt.close(fig)

print("REPORT_stability.txt, S8_stabilnost.png, S9_primery.png")
