import csv, os, statistics as st
from collections import defaultdict
H=os.path.dirname(os.path.abspath(__file__))
R=[]
for r in csv.DictReader(open(os.path.join(H,"compare.csv"),encoding="utf-8")):
    r['det_old']=int(r['det_old']); r['det_new']=int(r['det_new'])
    for k in ('snr_old','snr_new','delay_old','delay_new','rate_old','rate_new'): r[k]=float(r[k])
    r['bproj']=float(r['bproj']) if r['bproj'] else None
    R.append(r)
BUG=lambda r: r['exp']=='raks08ad' and r['band']=='C' and r['ref']=='Ef'
bug=[r for r in R if BUG(r)]; D=[r for r in R if not BUG(r)]
SPACE=lambda r: r['sta']=='RA' or r['ref']=='RA'

L=[];A=L.append
A("="*104)
A("СТАТИСТИКА ЗАНОВО. Прямое сравнение: старый прогон (без _p) -> новый (_p).")
A("Корреляция = строка CLOCK не закомментирована '#'. Одни и те же базы в одни и те же моменты.")
A("ИСКЛЮЧЕНО: raks08ad C, опорная Ef (%d записей) - там полином был затёрт соседним диапазоном."%len(bug))
A("="*104)
A("Записей в работе: %d   (наземные %d, с РадиоАстроном %d)"%(
    len(D),sum(1 for r in D if not SPACE(r)),sum(1 for r in D if SPACE(r))))
bo=sum(r['det_old'] for r in D); bn=sum(r['det_new'] for r in D)
kept=[r for r in D if r['det_old'] and r['det_new']]
lost=[r for r in D if r['det_old'] and not r['det_new']]
gain=[r for r in D if not r['det_old'] and r['det_new']]
A("")
A("1. КОРРЕЛЯЦИЯ ЕСТЬ / НЕТ")
A("   было %d  ->  стало %d   (%+d)"%(bo,bn,bn-bo))
A("   в обеих %d | ПРОПАЛА %d | ПОЯВИЛАСЬ %d | нет ни в одной %d"%(
    len(kept),len(lost),len(gain),len(D)-len(kept)-len(lost)-len(gain)))
A("")
A("2. ЗАВИСИМОСТЬ ОТ ПРОЕКЦИИ БАЗЫ (в диаметрах Земли)")
BINS=[(0,0.25,"0.00-0.25"),(0.25,0.5,"0.25-0.50"),(0.5,1.0,"0.50-1.00"),
      (1.0,3.0,"1.0-3.0"),(3.0,6.0,"3.0-6.0"),(6.0,99,"6.0-10")]
A("   %-11s %6s | корр. была  стала | пропало появилось | доля детекта было->стало"%("проекция","всего"))
for lo,hi,lab in BINS:
    s=[r for r in D if r['bproj'] is not None and lo<=r['bproj']<hi]
    if not s: continue
    o=sum(r['det_old'] for r in s); n=sum(r['det_new'] for r in s)
    l=sum(1 for r in s if r['det_old'] and not r['det_new']); g=sum(1 for r in s if not r['det_old'] and r['det_new'])
    A("   %-11s %6d | %10d %6d | %7d %10d | %5.0f%% -> %.0f%%"%(lab,len(s),o,n,l,g,100*o/len(s),100*n/len(s)))
A("")
A("   ПРОПАВШИЕ и ПОЯВИВШИЕСЯ по проекции базы:")
A("   %-9s %-2s %-4s %-4s %-9s | проекция, D | SNR было -> стало"%("сеанс","д","оп","ст","время"))
for lab,arr in [("ПРОПАЛА",lost),("ПОЯВИЛАСЬ",gain)]:
    A("   --- %s ---"%lab)
    for r in sorted(arr,key=lambda r:r['bproj'] or 0):
        A("   %-9s %-2s %-4s %-4s %-9s | %10.2f | %8.2f -> %8.2f"%(
            r['exp'],r['band'],r['ref'],r['sta'],r['t'],r['bproj'],r['snr_old'],r['snr_new']))
A("")
A("3. SNR, ЗАДЕРЖКА, FRINGE RATE — прямой счёт по %d записям с корреляцией в обеих"%len(kept))
def direct(nm,fo,fn,unit="",absv=True):
    o=[abs(fo(r)) if absv else fo(r) for r in kept]; n=[abs(fn(r)) if absv else fn(r) for r in kept]
    up=sum(1 for a,b in zip(o,n) if b>a); dn=sum(1 for a,b in zip(o,n) if b<a); eq=len(o)-up-dn
    A("   %-22s стало БОЛЬШЕ у %3d, МЕНЬШЕ у %3d, без изм. %3d %s"%(nm,up,dn,eq,unit))
    return o,n
direct("SNR",lambda r:r['snr_old'],lambda r:r['snr_new'],absv=False)
direct("|остаточная задержка|",lambda r:r['delay_old'],lambda r:r['delay_new'])
direct("|остаточный rate|",lambda r:r['rate_old'],lambda r:r['rate_new'])
A("")
A("   SNR по проекции базы (сколько записей выросло / упало):")
A("   %-11s %6s | SNR вырос  упал | сильных(>20): вырос упал"%("проекция","всего"))
for lo,hi,lab in BINS:
    s=[r for r in kept if r['bproj'] is not None and lo<=r['bproj']<hi]
    if not s: continue
    up=sum(1 for r in s if r['snr_new']>r['snr_old']); dn=len(s)-up
    ss=[r for r in s if r['snr_old']>20]
    su=sum(1 for r in ss if r['snr_new']>r['snr_old'])
    A("   %-11s %6d | %9d %5d | %12d %5d"%(lab,len(s),up,dn,su,len(ss)-su))
txt="\n".join(L)
open(os.path.join(H,"REPORT.txt"),"w",encoding="utf-8").write(txt)
print(txt)
