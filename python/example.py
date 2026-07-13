#!/usr/bin/env python
# example.py — использование модели задержки из Python (нативный модуль ariadna).
#
# Сборка модуля:  sh python/build_pybind.sh   (или  .\python\build_pybind.ps1)
# Запуск (из корня репозитория, где лежит ariadna.pyd):  python python/example.py

import os, sys
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)                # чтобы найти ariadna.pyd в корне репозитория
import ariadna

EPH = os.path.join(REPO, "external/dephem-master/linux_p1550p2650.440t")
EOP = os.path.join(REPO, "external/catalogs/EOPC04_14_IAU2000_62-now.cat")
CFX = os.path.join(REPO, "example/GVLBI_RAKS01A0_L_20140423T130000_ASC_V1.cfx")
SCF = os.path.join(REPO, "example/RA140423_1200_v02.scf")   # или "" -> орбита из cfx (ORB_FILE)

ariadna.init_ephemeris(EPH)             # ОДИН раз за процесс, до расчёта

# Полиномы задержки и координат (u,v,w) ВСЕХ станций сеанса В ПАМЯТИ, без файлов.
res = ariadna.compute_task_polys(CFX, SCF, eop=EOP, block=60, degree=5, sample=6.0, tropo=True)

print("Станций:", len(res.delay))
for sp, uv in zip(res.delay, res.uvw):
    print(f"\n{sp.telescope}: блоков задержки {len(sp.blocks)}, блоков uvw {len(uv.blocks)}")
    b = sp.blocks[0]
    print(f"  первый блок: источник {b.source}, {b.utc_start:.6f}..{b.utc_stop:.6f} сут")
    print(f"  задержка P0..P5 (с): " + ", ".join(f"{c:.4e}" for c in b.coef))
    ub = uv.blocks[0]
    print(f"  u,v,w P0 (с): {ub.u[0]:.6e}, {ub.v[0]:.6e}, {ub.w[0]:.6e}")

# Значение полинома в момент t (сек от начала блока): tau(t) = sum P_k * t^k.
def eval_poly(coef, t):
    return sum(c * t ** k for k, c in enumerate(coef))

sp0 = res.delay[0]; blk = sp0.blocks[0]
print()
for t in (0.0, 30.0, 59.0):
    print(f"{sp0.telescope}: задержка в +{t:>4.0f} с блока = {eval_poly(blk.coef, t):.9e} с")

# Коэффициенты — обычные списки float; при желании -> numpy: import numpy as np; np.array(blk.coef)
assert len(res.delay) == len(res.uvw) and len(res.delay) > 0, "нет результатов"
print("\nOK: модуль ariadna работает.")
