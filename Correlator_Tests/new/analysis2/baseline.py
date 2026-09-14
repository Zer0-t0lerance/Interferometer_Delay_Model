import csv, os, re, glob, datetime
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # Correlator_Tests/new
H=os.path.join(ROOT,"analysis2"); PO=os.path.join(H,"poly")
C=2.99792458e8; DE=12742000.0        # диаметр Земли, м

def with_band(fn,band):
    b,e=os.path.splitext(fn)
    return (b if b.endswith("_"+band) else b+"_"+band)+e

def cfx_map(path):
    """iam -> базовое имя POLY_FILE"""
    m={}; iam=poly=None
    for ln in open(path,encoding="utf-8",errors="replace"):
        s=ln.strip()
        if s.lower().startswith("name"):
            if iam and poly: m[iam]=poly
            iam=poly=None
        elif s.lower().startswith("iam_name"): iam=s.split("=",1)[1].strip()
        elif s.lower().startswith("poly_file"): poly=s.split("=",1)[1].strip().split(":")[-1].strip()
    if iam and poly: m[iam]=poly
    return m

re_st=re.compile(r'start\s*=\s*(\d+)/(\d+)/(\d+)\s+(\d+)h(\d+)m(\d+)s')
re_sp=re.compile(r'stop\s*=\s*(\d+)/(\d+)/(\d+)\s+(\d+)h(\d+)m(\d+)s')
def load_uvw(path):
    """-> [(start_dt, stop_dt, [u коэф], [v коэф])]"""
    blocks=[]; st=sp=None; u=[]; v=[]
    for ln in open(path,encoding="utf-8",errors="replace"):
        a=re_st.search(ln)
        if a:
            if st and u: blocks.append((st,sp,u,v))
            g=list(map(int,a.groups())); st=datetime.datetime(g[2],g[1],g[0],g[3],g[4],g[5]); u=[];v=[]; continue
        b=re_sp.search(ln)
        if b:
            g=list(map(int,b.groups())); sp=datetime.datetime(g[2],g[1],g[0],g[3],g[4],g[5]); continue
        if ln.startswith("P"):
            p=[float(x) for x in ln.split("=")[1].split(",")]
            u.append(p[0]); v.append(p[1])
    if st and u: blocks.append((st,sp,u,v))
    return blocks

def uv_at(blocks,t):
    for st,sp,u,v in blocks:
        if st<=t<(sp if sp and sp>st else st+datetime.timedelta(seconds=60)):
            dt=(t-st).total_seconds()
            return (sum(c*dt**k for k,c in enumerate(u)), sum(c*dt**k for k,c in enumerate(v)))
    return None

re_dt=re.compile(r'(\d+)d(\d+)m(\d+)y(\d+)h(\d+)m(\d+)s')
UV={}; MAP={}
for d in sorted(glob.glob(os.path.join(ROOT,"*"))):
    if not os.path.isdir(d) or os.path.basename(d).startswith("analysis"): continue
    exp=os.path.basename(d)
    for cfx in glob.glob(os.path.join(d,"*.cfx")):
        if cfx.endswith("_p.cfx"): continue
        mb=re.search(r'_([CLKSXQ])_',os.path.basename(cfx))
        if not mb: continue
        band=mb.group(1); MAP[(exp,band)]=cfx_map(cfx)
        for iam,pf in MAP[(exp,band)].items():
            base,ext=os.path.splitext(with_band(pf,band))
            f=os.path.join(PO,exp,base+"_uvw"+ext)
            if os.path.exists(f): UV[(exp,band,iam)]=load_uvw(f)

R=list(csv.DictReader(open(os.path.join(H,"compare.csv"),encoding="utf-8")))
miss=0
for r in R:
    g=list(map(int,re_dt.match(r['tfull']).groups()))
    t=datetime.datetime(g[2],g[1],g[0],g[3],g[4],g[5])
    a=UV.get((r['exp'],r['band'],r['sta'])); b=UV.get((r['exp'],r['band'],r['ref']))
    r['bproj']=""
    if a and b:
        pa,pb=uv_at(a,t),uv_at(b,t)
        if pa and pb:
            r['bproj']="%.4f"%(((pa[0]-pb[0])**2+(pa[1]-pb[1])**2)**.5*C/DE)
    if not r['bproj']: miss+=1
cols=list(R[0].keys())
with open(os.path.join(H,"compare.csv"),"w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=cols); w.writeheader(); w.writerows(R)
bp=[float(r['bproj']) for r in R if r['bproj']]
print("проекция базы посчитана для %d из %d записей (не удалось: %d)"%(len(bp),len(R),miss))
print("диапазон проекции базы: %.2f .. %.2f диаметров Земли"%(min(bp),max(bp)))
import statistics as st
print("наземные базы: до %.2f;  с РадиоАстроном: %.2f .. %.2f"%(
    max(float(r['bproj']) for r in R if r['bproj'] and r['sta']!='RA' and r['ref']!='RA'),
    min(float(r['bproj']) for r in R if r['bproj'] and (r['sta']=='RA' or r['ref']=='RA')),
    max(float(r['bproj']) for r in R if r['bproj'] and (r['sta']=='RA' or r['ref']=='RA'))))
