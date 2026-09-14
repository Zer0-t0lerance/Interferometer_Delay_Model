#!/usr/bin/env python
# uvwlib.py - общий код для вычисления проекции базы из полиномов uvw.
# Используется и baseline.py (сравнение по базам), и ra_scan.py (счёт по сканам).
import os, re, glob, datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # Correlator_Tests/new
H    = os.path.join(ROOT, "analysis2")
PO   = os.path.join(H, "poly")
C    = 2.99792458e8
DE   = 12742000.0        # диаметр Земли, м

re_st = re.compile(r'start\s*=\s*(\d+)/(\d+)/(\d+)\s+(\d+)h(\d+)m(\d+)s')
re_sp = re.compile(r'stop\s*=\s*(\d+)/(\d+)/(\d+)\s+(\d+)h(\d+)m(\d+)s')
re_dt = re.compile(r'(\d+)d(\d+)m(\d+)y(\d+)h(\d+)m(\d+)s')


def with_band(fn, band):
    b, e = os.path.splitext(fn)
    return (b if b.endswith("_" + band) else b + "_" + band) + e


def parse_tfull(s):
    g = list(map(int, re_dt.match(s).groups()))
    return datetime.datetime(g[2], g[1], g[0], g[3], g[4], g[5])


def cfx_map(path):
    """iam_name -> базовое имя POLY_FILE"""
    m = {}; iam = poly = None
    for ln in open(path, encoding="utf-8", errors="replace"):
        s = ln.strip()
        if s.lower().startswith("name"):
            if iam and poly: m[iam] = poly
            iam = poly = None
        elif s.lower().startswith("iam_name"):
            iam = s.split("=", 1)[1].strip()
        elif s.lower().startswith("poly_file"):
            poly = s.split("=", 1)[1].strip().split(":")[-1].strip()
    if iam and poly: m[iam] = poly
    return m


def load_uvw(path):
    """-> [(start_dt, stop_dt, [коэф. u], [коэф. v])]"""
    blocks = []; st = sp = None; u = []; v = []
    for ln in open(path, encoding="utf-8", errors="replace"):
        a = re_st.search(ln)
        if a:
            if st and u: blocks.append((st, sp, u, v))
            g = list(map(int, a.groups()))
            st = datetime.datetime(g[2], g[1], g[0], g[3], g[4], g[5]); u = []; v = []
            continue
        b = re_sp.search(ln)
        if b:
            g = list(map(int, b.groups()))
            sp = datetime.datetime(g[2], g[1], g[0], g[3], g[4], g[5]); continue
        if ln.startswith("P"):
            p = [float(x) for x in ln.split("=")[1].split(",")]
            u.append(p[0]); v.append(p[1])
    if st and u: blocks.append((st, sp, u, v))
    return blocks


def uv_at(blocks, t):
    for st, sp, u, v in blocks:
        if st <= t < (sp if sp and sp > st else st + datetime.timedelta(seconds=60)):
            dt = (t - st).total_seconds()
            return (sum(c * dt ** k for k, c in enumerate(u)),
                    sum(c * dt ** k for k, c in enumerate(v)))
    return None


def build_index():
    """(exp, band, iam_name) -> блоки полинома uvw"""
    UV = {}
    for d in sorted(glob.glob(os.path.join(ROOT, "*"))):
        if not os.path.isdir(d) or os.path.basename(d).startswith("analysis"): continue
        exp = os.path.basename(d)
        for cfx in glob.glob(os.path.join(d, "*.cfx")):
            if cfx.endswith("_p.cfx"): continue
            mb = re.search(r'_([CLKSXQ])_', os.path.basename(cfx))
            if not mb: continue
            band = mb.group(1)
            for iam, pf in cfx_map(cfx).items():
                base, ext = os.path.splitext(with_band(pf, band))
                f = os.path.join(PO, exp, base + "_uvw" + ext)
                if os.path.exists(f):
                    UV[(exp, band, iam)] = load_uvw(f)
    return UV


def bproj(UV, exp, band, sta, ref, t):
    """Проекция базы sta-ref в диаметрах Земли, или None."""
    a = UV.get((exp, band, sta)); b = UV.get((exp, band, ref))
    if not (a and b): return None
    pa, pb = uv_at(a, t), uv_at(b, t)
    if not (pa and pb): return None
    return ((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2) ** .5 * C / DE
