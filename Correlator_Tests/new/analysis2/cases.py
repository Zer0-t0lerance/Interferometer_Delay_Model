import csv, os
from collections import defaultdict
H=os.path.dirname(os.path.abspath(__file__))
R=[]
for r in csv.DictReader(open(os.path.join(H,"compare.csv"),encoding="utf-8")):
    r['det_old']=int(r['det_old']); r['det_new']=int(r['det_new'])
    for k in ('snr_old','snr_new'): r[k]=float(r[k])
    r['bproj']=float(r['bproj']) if r['bproj'] else None
    R.append(r)
BUG=lambda r: r['exp']=='raks08ad' and r['band']=='C' and r['ref']=='Ef'
D=[r for r in R if not BUG(r)]
RUN=defaultdict(list); BASE=defaultdict(list)
for r in D:
    RUN[(r['exp'],r['band'])].append(r)
    BASE[(r['exp'],r['band'],r['pol'],r['ref'],r['sta'])].append(r)
def tot(a): return sum(x['det_old'] for x in a), sum(x['det_new'] for x in a)

L=[];A=L.append
A("="*110)
A("РАЗБОР КАЖДОГО СЛУЧАЯ: пропала / появилась корреляция")
A("Вопрос: это весь прогон диапазона, вся база, или один интервал?")
A("="*110)
rows=[]
for kind,sel in [("ПРОПАЛА",lambda r: r['det_old'] and not r['det_new']),
                 ("ПОЯВИЛАСЬ",lambda r: not r['det_old'] and r['det_new'])]:
    A(""); A("### %s ###"%kind)
    for r in sorted([x for x in D if sel(x)],key=lambda x:(x['exp'],x['band'],x['t'])):
        rk=(r['exp'],r['band']); bk=(r['exp'],r['band'],r['pol'],r['ref'],r['sta'])
        ro,rn=tot(RUN[rk]); bo,bn=tot(BASE[bk]); nb=len(BASE[bk]); nr=len(RUN[rk])
        if rn==0 and ro>0:   verdict="ВЕСЬ ПРОГОН потерял корреляцию"
        elif bn==0 and bo>0: verdict="ВСЯ БАЗА потеряла корреляцию (прогон цел)"
        elif bo==0 and bn>0: verdict="ВСЯ БАЗА приобрела корреляцию"
        else:                verdict="ТОЛЬКО ЭТОТ ИНТЕРВАЛ (база и прогон целы)"
        A("")
        A("  %s %s  база %s-%s  %s   проекция %.2f диам.Земли   SNR %.2f -> %.2f"%(
            r['exp'],r['band'],r['sta'],r['ref'],r['t'],r['bproj'],r['snr_old'],r['snr_new']))
        A("     прогон %s %s : записей %d, корреляций %d -> %d"%(r['exp'],r['band'],nr,ro,rn))
        A("     база   %s-%s : интервалов %d, корреляций %d -> %d"%(r['sta'],r['ref'],nb,bo,bn))
        A("     ВЫВОД: %s"%verdict)
        rows.append((kind,r,verdict,nr,ro,rn,nb,bo,bn))
A(""); A("="*110); A("ИТОГ ПО ТИПАМ СЛУЧАЕВ"); A("="*110)
cnt=defaultdict(int)
for kind,r,v,*_ in rows: cnt[(kind,v)]+=1
for (kind,v),n in sorted(cnt.items()): A("  %-10s %-45s : %d"%(kind,v,n))
A("")
A("  Прогонов, потерявших корреляцию ЦЕЛИКОМ  : %d"%sum(1 for k,v in RUN.items() if tot(v)[0]>0 and tot(v)[1]==0))
A("  Прогонов, приобретших корреляцию ЦЕЛИКОМ : %d"%sum(1 for k,v in RUN.items() if tot(v)[0]==0 and tot(v)[1]>0))
A("  Баз, потерявших корреляцию целиком       : %d"%sum(1 for k,v in BASE.items() if tot(v)[0]>0 and tot(v)[1]==0))
A("  Баз, приобретших корреляцию целиком      : %d"%sum(1 for k,v in BASE.items() if tot(v)[0]==0 and tot(v)[1]>0))
txt="\n".join(L); open(os.path.join(H,"REPORT_cases.txt"),"w",encoding="utf-8").write(txt)
