import re, os, glob, csv
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT  = os.path.join(ROOT, "analysis2")

re_hdr = re.compile(r'#\s*(\w+)-POL\s+Time interval\s*=\s*\[(\d+:\d+:\d+)\s*-\s*(\d+:\d+:\d+)\[\s*reference antenna:\s*(\S+)')
re_clk = re.compile(r'^(#?)\s*CLOCK\s*=\s*(\S+?),\s*(\S+?),\s*([-\d.eE+]+),\s*([-\d.eE+]+),\s*([-\d.eE+]+)')
re_snr = re.compile(r'#\s*snr\s*=\s*([-\d.eE+]+)\s*,\s*IF\s*=\s*(-?\d+)')

def parse(path):
    """-> список записей. detected = строка CLOCK НЕ закомментирована '#'."""
    out=[]; pol=ref=t0=None
    for ln in open(path, encoding="utf-8", errors="replace").read().replace("\r\n","\n").split("\n"):
        h = re_hdr.search(ln)
        if h: pol, t0, ref = h.group(1), h.group(2), h.group(4); continue
        c = re_clk.match(ln)
        if c: out.append(dict(pol=pol, ref=ref, t=t0, tfull=c.group(3), sta=c.group(2),
                              detected=(c.group(1) != "#"),
                              delay=float(c.group(4)), rate=float(c.group(5)), acc=float(c.group(6)),
                              snr=None, iff=None, _raw=ln)); continue
        s = re_snr.search(ln)
        if s and out and out[-1]['snr'] is None:
            out[-1]['snr']=float(s.group(1)); out[-1]['iff']=int(s.group(2))
    return out

rows=[]; report=[]
for d in sorted(glob.glob(os.path.join(ROOT,"*"))):
    if not os.path.isdir(d) or os.path.basename(d).startswith("analysis"): continue
    exp=os.path.basename(d)
    for old in sorted(glob.glob(os.path.join(d,"*_ffr.txt"))):
        if old.endswith("_p_ffr.txt"): continue
        new = old[:-len("_ffr.txt")] + "_p_ffr.txt"
        if not os.path.exists(new): report.append("НЕТ ПАРЫ: "+os.path.basename(old)); continue
        band = (re.search(r'_([CLKSXQ])_', os.path.basename(old)) or [None,'?'])[1]
        O, N = parse(old), parse(new)
        # ключ: поляризация + опорная + станция + время интервала
        ko = {(r['pol'],r['ref'],r['sta'],r['t']): r for r in O}
        kn = {(r['pol'],r['ref'],r['sta'],r['t']): r for r in N}
        if len(ko)!=len(O) or len(kn)!=len(N): report.append("ДУБЛИ КЛЮЧЕЙ: %s %s"%(exp,band))
        only_o = sorted(set(ko)-set(kn)); only_n = sorted(set(kn)-set(ko))
        if only_o or only_n:
            report.append("РАСХОЖДЕНИЕ НАБОРА %s %s: только в старом %d, только в новом %d"%(exp,band,len(only_o),len(only_n)))
            for k in only_o[:3]: report.append("    только старый: %s"%(k,))
            for k in only_n[:3]: report.append("    только новый : %s"%(k,))
        for k in sorted(set(ko)&set(kn)):
            o,n = ko[k],kn[k]
            # запись опорной станции (сама с собой) - не база, пропускаем
            if o['sta']==o['ref'] or o['iff']==-100: continue
            rows.append(dict(exp=exp,band=band,pol=k[0],ref=k[1],sta=k[2],t=k[3],tfull=o['tfull'],
                             det_old=int(o['detected']), det_new=int(n['detected']),
                             snr_old=o['snr'], snr_new=n['snr'],
                             delay_old=o['delay'], delay_new=n['delay'],
                             rate_old=o['rate'], rate_new=n['rate'],
                             acc_old=o['acc'], acc_new=n['acc']))
cols=list(rows[0].keys())
with open(os.path.join(OUT,"compare.csv"),"w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=cols); w.writeheader(); w.writerows(rows)
print("сравнимых записей (одна база в один момент): %d"%len(rows))
print("проблемы разбора:", "нет" if not report else "")
for x in report: print("   "+x)
