import csv, os
from collections import defaultdict
H=os.path.dirname(os.path.abspath(__file__))
R=list(csv.DictReader(open(os.path.join(H,"compare.csv"),encoding="utf-8")))
for r in R:
    r['det_old']=int(r['det_old']); r['det_new']=int(r['det_new'])
    for k in ('snr_old','snr_new','delay_old','delay_new','rate_old','rate_new'): r[k]=float(r[k])

L=[]; A=L.append
A("="*100)
A("ПРЯМОЕ СРАВНЕНИЕ: старый прогон (без _p)  ->  новый прогон (с _p)")
A("Корреляция = строка CLOCK НЕ закомментирована '#'. Сравниваются одни и те же базы в одни и те же моменты.")
A("="*100)
A("Сравнимых записей (база x момент): %d"%len(R))
bo=sum(r['det_old'] for r in R); bn=sum(r['det_new'] for r in R)
A("  корреляция БЫЛА  : %d"%bo)
A("  корреляция СТАЛА : %d"%bn)
A("  изменение        : %+d"%(bn-bo))
A("")
kept =[r for r in R if r['det_old'] and r['det_new']]
lost =[r for r in R if r['det_old'] and not r['det_new']]
gain =[r for r in R if not r['det_old'] and r['det_new']]
none_=[r for r in R if not r['det_old'] and not r['det_new']]
A("  корреляция в обеих   : %d"%len(kept))
A("  ПРОПАЛА              : %d"%len(lost))
A("  ПОЯВИЛАСЬ            : %d"%len(gain))
A("  нет ни в одной       : %d"%len(none_))
A("")
A("-"*100); A("ВСЕ ЗАПИСИ, ГДЕ КОРРЕЛЯЦИЯ ПРОПАЛА (%d)"%len(lost)); A("-"*100)
A("  %-9s %-2s %-3s %-4s %-4s %-9s | SNR было -> стало"%("сеанс","д","поl","оп","ст","время"))
for r in sorted(lost,key=lambda r:(r['exp'],r['band'],r['t'])):
    A("  %-9s %-2s %-4s %-4s %-4s %-9s | %8.2f -> %8.2f"%(
        r['exp'],r['band'],r['pol'],r['ref'],r['sta'],r['t'],r['snr_old'],r['snr_new']))
A("")
A("-"*100); A("ВСЕ ЗАПИСИ, ГДЕ КОРРЕЛЯЦИЯ ПОЯВИЛАСЬ (%d)"%len(gain)); A("-"*100)
A("  %-9s %-2s %-4s %-4s %-4s %-9s | SNR было -> стало"%("сеанс","д","пол","оп","ст","время"))
for r in sorted(gain,key=lambda r:(r['exp'],r['band'],r['t'])):
    A("  %-9s %-2s %-4s %-4s %-4s %-9s | %8.2f -> %8.2f"%(
        r['exp'],r['band'],r['pol'],r['ref'],r['sta'],r['t'],r['snr_old'],r['snr_new']))
A("")
A("-"*100); A("ПО СЕАНСАМ: сколько записей с корреляцией было и стало"%()); A("-"*100)
A("  %-9s %-2s | записей | корреляция была | стала | пропала | появилась"%("сеанс","д"))
g=defaultdict(list)
for r in R: g[(r['exp'],r['band'])].append(r)
tot_l=tot_g=0
for k in sorted(g):
    v=g[k]; o=sum(x['det_old'] for x in v); n=sum(x['det_new'] for x in v)
    l=sum(1 for x in v if x['det_old'] and not x['det_new']); gg=sum(1 for x in v if not x['det_old'] and x['det_new'])
    tot_l+=l; tot_g+=gg
    mark="" if (l==0 and gg==0) else ("  <-- изменения" )
    A("  %-9s %-2s | %7d | %15d | %5d | %7d | %9d%s"%(k[0],k[1],len(v),o,n,l,gg,mark))
A("  %-9s %-2s | %7d | %15d | %5d | %7d | %9d"%("ИТОГО","",len(R),bo,bn,tot_l,tot_g))
txt="\n".join(L)
open(os.path.join(H,"REPORT_detection.txt"),"w",encoding="utf-8").write(txt)
print(txt)
