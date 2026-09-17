#!/usr/bin/env python
# verify_batch.py - проверка результата пакетного прогона:
#   1) у каждого исходного задания есть <имя>_p.cfx рядом с ним;
#   2) метки времени в TIMEOFS нового задания СОВПАДАЮТ с метками исходного;
#   3) строка вывода коррелятора переименована в *_p.uvx;
#   4) полиномы лежат в <эксперимент>/new_poly и покрывают все станции задания.
#
# Запуск (из корня репозитория):
#   python tools/verify_batch.py [корень]       # по умолчанию корень = Correlator_Tests
import os, re, glob, sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else "Correlator_Tests"
if not os.path.isdir(ROOT):
    sys.exit("Нет такой папки: %s (запускать из корня репозитория)" % ROOT)

def blocks(path):
    out, cur = [], None
    for raw in open(path, encoding="utf-8", errors="replace"):
        t = raw.rstrip("\r\n").strip()
        if t.startswith("[$TLSC]"): cur = []; out.append(cur); continue
        if t.startswith("[$"): cur = None; continue
        if cur is not None: cur.append(t)
    return out

def space_block(path):
    for b in blocks(path):
        if any(x.startswith("ORB_FILE") for x in b): return b
    return []

def marks(path):
    m = {}
    for t in space_block(path):
        g = re.match(r'TIMEOFS(\d+)\s*=\s*([^,]+),\s*([-\d.eE+]+)', t)
        if g: m[g.group(1)] = (float(g.group(2)), float(g.group(3)))
    return m

def poly_names(path):
    names = []
    for b in blocks(path):
        for t in b:
            g = re.match(r'POLY_FILE\s*=\s*(.*)', t)
            if g: names.append(g.group(1).strip().split(":")[-1].strip())
    return names

def out_file(path):
    for raw in open(path, encoding="utf-8", errors="replace"):
        t = raw.rstrip("\r\n").strip()
        k = "".join(c for c in t[:10] if not c.isspace()).upper()
        if k.startswith("OUTFILE"): return t.split("=", 1)[1].strip()
    return ""

bad, note = [], []
n_task = n_mark = n_same = n_uvx = n_poly = n_polymiss = 0
for d in sorted(glob.glob(os.path.join(ROOT, "*"))):
    if not os.path.isdir(d) or os.path.basename(d) == "new": continue
    exp = os.path.basename(d)
    for src in sorted(glob.glob(os.path.join(d, "*.cfx"))):
        if src.endswith("_p.cfx"): continue
        n_task += 1
        dst = src[:-4] + "_p.cfx"
        if not os.path.exists(dst):
            bad.append("%s: нет %s" % (exp, os.path.basename(dst))); continue

        ms, md = marks(src), marks(dst)
        for idx, (dly, mk) in md.items():
            n_mark += 1
            if idx in ms:
                if abs(ms[idx][1] - mk) < 1e-9: n_same += 1
                else: bad.append("%s %s: метка TIMEOFS%s %.9f -> %.9f" %
                                 (exp, os.path.basename(dst), idx, ms[idx][1], mk))
            else:
                # Штатный случай: строки TIMEOFS в задании не было, метку взяли из имени
                # файла данных. Это не ошибка, но показать надо (см. tools/batch_cfx.cpp).
                note.append("%s %s: TIMEOFS%s добавлен — в исходном строки не было, "
                            "метка из имени файла" % (exp, os.path.basename(dst), idx))
        for idx in ms:
            if idx not in md:
                bad.append("%s %s: TIMEOFS%s пропал (в исходном был)" %
                           (exp, os.path.basename(dst), idx))

        o = out_file(dst)
        if o.lower().endswith("_p.uvx"): n_uvx += 1
        else: bad.append("%s %s: вывод коррелятора не переименован: %s" % (exp, os.path.basename(dst), o))

        npd = os.path.join(d, "new_poly")
        for nm in poly_names(dst):
            n_poly += 1
            f = os.path.join(npd, nm)
            if not os.path.exists(f):
                # регистр имени в cfx и на диске может отличаться
                alt = [x for x in os.listdir(npd) if x.lower() == nm.lower()] if os.path.isdir(npd) else []
                if not alt: bad.append("%s: нет полинома new_poly/%s" % (exp, nm)); n_polymiss += 1

print("заданий проверено: %d" % n_task)
print("строк TIMEOFS в новых заданиях: %d, метка совпала с исходной: %d" % (n_mark, n_same))
print("заданий с выводом *_p.uvx: %d" % n_uvx)
print("ссылок на полиномы: %d, отсутствует файлов: %d" % (n_poly, n_polymiss))
print()
if note:
    print("ШТАТНЫЕ ОТЛИЧИЯ (%d) — не ошибки:" % len(note))
    for b in note: print("  " + b)
    print()
if bad:
    print("ПРОБЛЕМЫ (%d):" % len(bad))
    for b in bad: print("  " + b)
else:
    print("проблем нет")
