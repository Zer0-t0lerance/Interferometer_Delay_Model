#!/usr/bin/env python
# check_marks.py - разведка перед пакетным прогоном: для каждого входного cfx смотрим
# строки FILExx космической станции и метки времени в строках TIMEOFSxx.
# Вопросы: у всех ли FILExx есть своя TIMEOFS, и совпадает ли метка с тем, что даёт имя файла.
#
# Запуск (из корня репозитория):
#   python tools/check_marks.py [корень]        # по умолчанию корень = Correlator_Tests
import os, re, glob, sys, datetime

ROOT = sys.argv[1] if len(sys.argv) > 1 else "Correlator_Tests"
if not os.path.isdir(ROOT):
    sys.exit("Нет такой папки: %s (запускать из корня репозитория)" % ROOT)

def ymd_doy_to_mjd(Y, DDD, hh, mm, ss):
    d = datetime.datetime(Y, 1, 1) + datetime.timedelta(days=DDD - 1,
                                                        hours=hh, minutes=mm, seconds=ss)
    return (d - datetime.datetime(1858, 11, 17)).total_seconds() / 86400.0

def mark_from_name(val):
    """Время из имени файла данных: YYYYDDDHHMMSS -> MJD UTC, как это делает модель."""
    for m in re.finditer(r'\d{13,}', val):
        d = m.group(0)[:13]
        return ymd_doy_to_mjd(int(d[:4]), int(d[4:7]), int(d[7:9]), int(d[9:11]), int(d[11:13]))
    return None

def scan(path):
    """-> (имя космической станции, {idx: значение FILExx}, {idx: метка из TIMEOFSxx})
    Блок станции разбирается ЦЕЛИКОМ: ORB_FILE стоит ПОСЛЕ строк FILExx, поэтому
    признак «космический блок» известен только к концу блока."""
    blocks, cur = [], None
    for raw in open(path, encoding="utf-8", errors="replace"):
        t = raw.strip()
        if t.startswith("[$TLSC]"):
            cur = []; blocks.append(cur); continue
        if t.startswith("[$") and cur is not None:
            cur = None
        if cur is not None:
            cur.append(t)
    for b in blocks:
        if not any(x.startswith("ORB_FILE") for x in b): continue
        station = ""
        files, marks = {}, {}
        for t in b:
            if t.lower().startswith("name") and not t.lower().startswith("name_"):
                station = t.split("=", 1)[1].strip()
            m = re.match(r'FILE(\d+)\s*=\s*(.*)', t)
            if m: files[m.group(1)] = m.group(2).strip()
            m = re.match(r'TIMEOFS(\d+)\s*=\s*([^,]+),\s*([-\d.eE+]+)', t)
            if m: marks[m.group(1)] = float(m.group(3))
        return station, files, marks
    return "", {}, {}

tot_files = tot_marks = miss = diff = 0
print("%-10s %-4s %-9s %5s %5s  %s" % ("сеанс", "диап", "станция", "FILE", "MARK", "замечания"))
for d in sorted(glob.glob(os.path.join(ROOT, "*"))):
    if not os.path.isdir(d) or os.path.basename(d) == "new": continue
    exp = os.path.basename(d)
    for f in sorted(glob.glob(os.path.join(d, "*.cfx"))):
        if f.endswith("_p.cfx"): continue
        band = (re.search(r'_([CLKSXQ])_', os.path.basename(f)) or [None, "?"])[1]
        sta, files, marks = scan(f)
        notes = []
        for idx in sorted(files):
            if idx not in marks:
                notes.append("нет TIMEOFS%s" % idx); miss += 1
            else:
                mn = mark_from_name(files[idx])
                if mn is None:
                    notes.append("имя FILE%s без времени" % idx)
                elif abs(mn - marks[idx]) > 1e-9:
                    notes.append("FILE%s: метка %.9f vs имя %.9f" % (idx, marks[idx], mn)); diff += 1
        extra = [i for i in marks if i not in files]
        if extra: notes.append("TIMEOFS без FILE: " + ",".join(sorted(extra)))
        tot_files += len(files); tot_marks += len(marks)
        print("%-10s %-4s %-9s %5d %5d  %s" % (exp, band, sta, len(files), len(marks),
                                               "; ".join(notes) if notes else "-"))
print()
print("ИТОГО: файлов данных космоса %d, меток в cfx %d, без метки %d, метка != имени файла %d"
      % (tot_files, tot_marks, miss, diff))
