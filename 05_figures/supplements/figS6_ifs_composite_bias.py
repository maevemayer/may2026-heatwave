import numpy as np, pandas as pd, xarray as xr
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.dates as mdates, os, matplotlib.colors as mcolors
fs = 17
axfs = fs-2   # x/y tick label size (offset from fs, not hardcoded)
plt.rcParams.update({"font.size":fs,"xtick.labelsize":axfs,"ytick.labelsize":axfs})
def _ipcc(name,fallback):
    for p in ["${HOME}/colormaps_IPCC/continuous_colormaps_rgb_0-1/"+name, os.path.join(os.path.dirname(os.path.abspath(__file__)),name)]:
        if os.path.exists(p): return mcolors.LinearSegmentedColormap.from_list(name,np.loadtxt(p))
    return plt.get_cmap(fallback)
TEMPDIV=_ipcc("temp_div.txt","RdBu_r")
# IFS: the colleague's tanom_decomp of 2026-10-19 (corrected advective T')
BASE=Path(os.environ.get("IFS_BASE","${IFS_TANOM_DECOMP}"))
HOURS=["00","06","12","18"]; TIMES=pd.date_range("2026-05-18 00:00","2026-05-30 18:00",freq="6h"); TIDX={t:i for i,t in enumerate(TIMES)}
P_EDGES=np.arange(175.,1026.,50.); P_CEN=0.5*(P_EDGES[:-1]+P_EDGES[1:]); TERMS=["TA","advective","adiabatic","diabatic"]
ERA5VAR={"TA":"TA","advective":"advection","adiabatic":"adiabatic","diabatic":"diabatic"}
IFSCOLS={"TA":["tanom0"],"advective":["tanom_advective"],"adiabatic":["tanom_adiab1","tanom_adiab2"],"diabatic":["tanom_diabatic"]}
def hov(V,P,T):
    df=pd.DataFrame({"v":V,"p":P,"t":T}).dropna(); df=df[(df.p>=P_EDGES[0])&(df.p<P_EDGES[-1])]; df["ti"]=df.t.map(TIDX); df=df.dropna(subset=["ti"]); df["pb"]=np.digitize(df.p,P_EDGES)-1
    g=df.groupby(["pb","ti"])["v"].mean(); G=np.full((len(P_CEN),len(TIMES)),np.nan)
    for (pb,ti),val in g.items():
        if 0<=int(pb)<len(P_CEN): G[int(pb),int(ti)]=val
    return G
# NEW pipeline: single ERA5 event origin file (vertical Hovmoller by ARRIVAL pressure/time)
era5={}
# ERA5: windfix + genesis as in Hotz et al. (2024)
d=xr.open_dataset(os.environ.get("POOLED_FILE","${DERIVED_DATA}/origin_contributions_6h_clim2017_2025_windfix/origin_contributions.nc"))
ap=d["arrival_p"].values; tt=pd.DatetimeIndex(d["time"].values).values
for term in TERMS:
    val=d[ERA5VAR[term]].values
    era5[term]=hov(val.ravel(),ap.ravel(),np.repeat(tt[:,None],val.shape[1],axis=1).ravel())
d.close()
csv=pd.read_csv("${DERIVED_DATA}/ifs_cutoff_activity_with_T2m_Z500.csv"); csv["initdir"]=csv["init"].str.replace("-","")+"_00"; csv["mem"]=csv["member"].map(lambda m:f"{int(m):02d}")
topT=list(zip(csv.nlargest(5,"t2m_C").initdir,csv.nlargest(5,"t2m_C").mem)); topZ=list(zip(csv.nlargest(5,"z500_gh").initdir,csv.nlargest(5,"z500_gh").mem))
def ifs_member(initd,mem):
    rows=[]
    for D in pd.date_range("2026-05-18","2026-05-30",freq="D"):
        for h in HOURS:
            p=BASE/initd/mem/f"tanom_{D:%Y%m%d}_{h}"
            if p.is_file(): df=pd.read_csv(p,sep=r"\s+"); df["_t"]=D+pd.Timedelta(hours=int(h)); rows.append(df)
    if not rows: return None
    d=pd.concat(rows,ignore_index=True)
    for c in ["p0"]+sum(IFSCOLS.values(),[]): d[c]=pd.to_numeric(d[c],errors="coerce")
    return {term:hov(d[IFSCOLS[term]].sum(axis=1).values,d.p0.values,d._t.values) for term in TERMS}
def comp(members):
    acc={t:[] for t in TERMS}
    for i,m in members:
        r=ifs_member(i,m)
        if r:
            for t in TERMS: acc[t].append(r[t])
    return {t:np.nanmean(np.array(acc[t]),axis=0) for t in TERMS}
compT=comp(topT); compZ=comp(topZ)
biasT={t:compT[t]-era5[t] for t in TERMS}; biasZ={t:compZ[t]-era5[t] for t in TERMS}

def odd_levels(dicts):
    v=np.nanpercentile(np.abs(np.concatenate([d[t].ravel() for d in dicts for t in TERMS])),99); vi=max(int(np.ceil(v)),3)
    if vi%2==0: vi+=1
    return np.arange(-vi,vi+2,2).astype(float)
lev_abs=odd_levels([compT,compZ]); lev_bias=odd_levels([biasT,biasZ])
ROWS=[("Top-5 T2m",compT,lev_abs),("Top-5 Z500",compZ,lev_abs),("Top-5 T2m − ERA5",biasT,lev_bias),("Top-5 Z500 − ERA5",biasZ,lev_bias)]

tnum=mdates.date2num(TIMES.to_pydatetime()); fig,axes=plt.subplots(4,4,figsize=(21,17),sharex=True,sharey=True)
cf_abs=cf_bias=None
for ri,(rlab,rd,lev) in enumerate(ROWS):
    for ci,term in enumerate(TERMS):
        ax=axes[ri,ci]; cf=ax.contourf(tnum,P_CEN,rd[term],levels=lev,cmap="RdBu_r",extend="both")
        ax.set_yscale("log"); ax.set_ylim(1000,200); ax.set_yticks([1000,800,600,400,200]); ax.minorticks_off(); ax.set_yticklabels([1000,800,600,400,200])
        ax.axvline(mdates.date2num(pd.Timestamp("2026-05-25")),color="k",ls="--",lw=1)
        if ri==0: ax.set_title("TA" if term=="TA" else term.capitalize(),fontsize=fs)
        if ci==0: ax.set_ylabel(f"{rlab}\nPressure (hPa)",fontsize=fs)
        ax.text(-0.02,1.04,f"({chr(97+ri*4+ci)})",transform=ax.transAxes,fontsize=fs+2,fontweight="bold",va="bottom",ha="left")   # letter OUTSIDE panel
    if ri==0: cf_abs=cf
    if ri==2: cf_bias=cf
cb1=fig.colorbar(cf_abs,ax=axes[0:2,:].ravel().tolist(),location="right",shrink=0.55,pad=0.012,ticks=lev_abs); cb1.set_label("[K]",fontsize=fs-2); cb1.ax.tick_params(labelsize=axfs)
cb2=fig.colorbar(cf_bias,ax=axes[2:4,:].ravel().tolist(),location="right",shrink=0.55,pad=0.012,ticks=lev_bias); cb2.set_label("Δ [K]  (bias vs ERA5)",fontsize=fs-2); cb2.ax.tick_params(labelsize=axfs)
for ax in axes[-1,:]:
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2)); ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d")); ax.tick_params(axis="x",labelrotation=90,labelsize=axfs); ax.set_xlabel("Arrival date",fontsize=fs)
OUT="${HOME_WORK}/figure3_supplement_composites_bias_windfix_ifs20261019.png"; fig.savefig(OUT,dpi=150,bbox_inches="tight")
fig.savefig(OUT.replace(".png",".pdf"),bbox_inches="tight"); print("saved",OUT,"+ pdf")
