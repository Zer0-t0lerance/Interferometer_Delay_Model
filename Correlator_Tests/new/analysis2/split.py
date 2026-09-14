import csv, os
from collections import defaultdict
H=os.path.dirname(os.path.abspath(__file__))
R=[]
for r in csv.DictReader(open(os.path.join(H,"compare.csv"),encoding="utf-8")):
    r['det_old']=int(r['det_old']); r['det_new']=int(r['det_new'])
    for k in ('snr_old','snr_new','delay_old','delay_new','rate_old','rate_new'): r[k]=float(r[k])
    r['bproj']=float(r['bproj']) if r['bproj'] else None
    R.append(r)
BUG=lambda r: r['exp']=='raks08ad' and r['band']=='C' and r['ref']=='Ef'
D=[r for r in R if not BUG(r)]
SPACE=lambda r: r['sta']=='RA' or r['ref']=='RA'
BK=lambda r:(r['exp'],r['band'],r['pol'],r['ref'],r['sta']); RK=lambda r:(r['exp'],r['band'])
L=[];A=L.append
A("="*104)
A("РАЗДЕЛЬНО: НАЗЕМНО-КОСМИЧЕСКИЕ базы (с РадиоАстроном) и НАЗЕМНЫЕ")
A("="*104)
for name,sel in [("НАЗЕМНО-КОСМИЧЕСКИЕ (РадиоАстрон)",SPACE),("НАЗЕМНЫЕ",lambda r: not SPACE(r))]:
    S=[r for r in D if sel(r)]
    RUN=defaultdict(list); BASE=defaultdict(list)
    for r in S: RUN[RK(r)].append(r); BASE[BK(r)].append(r)
    tot=lambda a:(sum(x['det_old'] for x in a),sum(x['det_new'] for x in a))
    o,n=tot(S)
    kept=[r for r in S if r['det_old'] and r['det_new']]
    lost=[r for r in S if r['det_old'] and not r['det_new']]
    gain=[r for r in S if not r['det_old'] and r['det_new']]
    A(""); A("### %s ###"%name)
    A("  записей (база x момент): %d   баз: %d   прогонов с такими базами: %d"%(len(S),len(BASE),len(RUN)))
    A("  корреляция: было %d -> стало %d   (%+d)"%(o,n,n-o))
    A("  в обеих %d | ПРОПАЛА %d | ПОЯВИЛАСЬ %d | нет ни в одной %d"%(
        len(kept),len(lost),len(gain),len(S)-len(kept)-len(lost)-len(gain)))
    bl=[k for k,v in BASE.items() if tot(v)[0]>0 and tot(v)[1]==0]
    bg=[k for k,v in BASE.items() if tot(v)[0]==0 and tot(v)[1]>0]
    rl=[k for k,v in RUN.items() if tot(v)[0]>0 and tot(v)[1]==0]
    rg=[k for k,v in RUN.items() if tot(v)[0]==0 and tot(v)[1]>0]
    A("  БАЗ потеряло корреляцию целиком: %d %s"%(len(bl),[ "%s %s %s-%s"%(k[0],k[1],k[4],k[3]) for k in bl]))
    A("  БАЗ приобрело корреляцию целиком: %d %s"%(len(bg),[ "%s %s %s-%s"%(k[0],k[1],k[4],k[3]) for k in bg]))
    A("  ПРОГОНОВ потеряло целиком: %d %s"%(len(rl),[ "%s %s"%k for k in rl]))
    A("  ПРОГОНОВ приобрело целиком: %d %s"%(len(rg),[ "%s %s"%k for k in rg]))
    if lost or gain:
        A("  список смен статуса (проекция базы в диаметрах Земли):")
        for lab,arr in [("ПРОПАЛА",lost),("ПОЯВИЛАСЬ",gain)]:
            for r in sorted(arr,key=lambda r:r['bproj'] or 0):
                bo,bn=tot(BASE[BK(r)]); ro,rn=tot(RUN[RK(r)])
                sc=("вся база" if (bn==0 and bo>0) or (bo==0 and bn>0) else "один интервал")
                A("    %-10s %-9s %-2s %-3s-%-3s %-9s B=%5.2f  SNR %6.2f->%6.2f  база %d->%d  прогон %d->%d  [%s]"%(
                    lab,r['exp'],r['band'],r['sta'],r['ref'],r['t'],r['bproj'],r['snr_old'],r['snr_new'],bo,bn,ro,rn,sc))
    up=sum(1 for r in kept if r['snr_new']>r['snr_old']); dn=len(kept)-up
    ru=sum(1 for r in kept if abs(r['rate_new'])>abs(r['rate_old'])); rd=sum(1 for r in kept if abs(r['rate_new'])<abs(r['rate_old']))
    du=sum(1 for r in kept if abs(r['delay_new'])>abs(r['delay_old'])); dd=sum(1 for r in kept if abs(r['delay_new'])<abs(r['delay_old']))
    A("  величины по %d записям с корреляцией в обеих:"%len(kept))
    A("     SNR             больше %3d / меньше %3d"%(up,dn))
    A("     |fringe rate|   больше %3d / меньше %3d"%(ru,rd))
    A("     |задержка|      больше %3d / меньше %3d"%(du,dd))
    if S and any(r['bproj'] for r in S):
        A("  проекция базы: %.2f .. %.2f диаметров Земли"%(
            min(r['bproj'] for r in S if r['bproj']), max(r['bproj'] for r in S if r['bproj'])))
txt="\n".join(L); open(os.path.join(H,"REPORT_split.txt"),"w",encoding="utf-8").write(txt)
