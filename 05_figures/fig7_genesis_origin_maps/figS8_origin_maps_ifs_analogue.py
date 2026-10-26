import numpy as np, pandas as pd, xarray as xr, glob, os
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.colors as mcolors
# --- rasterize pcolormesh/pcolor in vector output (PDF/SVG) so maps match the PNG (no inter-cell seams) ---
import matplotlib.axes as _mpl_axes_raster, matplotlib.pyplot as _plt_raster
for _pcm_name in ("pcolormesh", "pcolor"):
    _pcm_orig = getattr(_mpl_axes_raster.Axes, _pcm_name)
    if not getattr(_pcm_orig, "_rasterized_patch", False):
        def _pcm_make(_orig):
            def _wrapped(self, *a, **k):
                k.setdefault("rasterized", True)
                _m = _orig(self, *a, **k)
                try: _m.set_edgecolor("face")
                except Exception: pass
                return _m
            _wrapped._rasterized_patch = True
            return _wrapped
        setattr(_mpl_axes_raster.Axes, _pcm_name, _pcm_make(_pcm_orig))
_plt_raster.rcParams["savefig.dpi"] = 200   # crisp rasterized meshes in the PDF (PNG keeps its explicit dpi)
# --- end rasterize patch ---

import cartopy.crs as ccrs, cartopy.feature as cfeature
def _ipcc(n,fb):
    p=os.path.join("${HOME_WORK}",n); return mcolors.LinearSegmentedColormap.from_list(n,np.loadtxt(p)) if os.path.exists(p) else plt.get_cmap(fb)
TEMPDIV=_ipcc("temp_div.txt","RdBu_r")
MOR=(-20.,0.,25.,40.); DLON=0.5; NTOP=5; EXT=[-90,25,10,70]; BAND_LO=(850,1000); BAND_HI=(250,750)
def wrap(x): return ((np.asarray(x)+180)%360)-180

# ================= IFS =================
# IFS: the colleague's tanom_decomp of 2026-10-19 (corrected advective T')
BASE=Path(os.environ.get("IFS_BASE","${IFS_TANOM_DECOMP}"))
HOURS=["00","06","12","18"]; DATES=["20260524","20260525"]
lonb=np.arange(-130,80+DLON,DLON); latb=np.arange(15,80+DLON,DLON); lonc=0.5*(lonb[:-1]+lonb[1:]); latc=0.5*(latb[:-1]+latb[1:])
tz=pd.read_csv("${DERIVED_DATA}/ifs_cutoff_activity_with_T2m_Z500.csv"); tz["initdir"]=tz["init"].str.replace("-","")+"_00"
# IFS origin: advective genesis (as the analogues); IFS_GENESIS=ta for the TA genesis
IFS_GENESIS=os.environ.get("IFS_GENESIS","adv")
GLON,GLAT=("lon_adv_genesis","lat_adv_genesis") if IFS_GENESIS=="adv" else ("lon_genesis","lat_genesis")
def member_maps(initd,mem):
    NEED=[GLON,GLAT,"p0","tanom_advective"]
    fr=[pd.read_csv(BASE/initd/mem/f"tanom_{Dd}_{h}",sep=r"\s+") for Dd in DATES for h in HOURS if (BASE/initd/mem/f"tanom_{Dd}_{h}").is_file()]
    if not fr: return None
    d=pd.concat(fr,ignore_index=True)
    for c in NEED: d[c]=pd.to_numeric(d[c],errors="coerce")
    d=d.dropna(subset=NEED); out={}
    for bnd in (BAND_LO,BAND_HI):
        dd=d[(d.p0>=bnd[0])&(d.p0<=bnd[1])]
        if len(dd)==0: out[bnd]=None; continue
        g,_,_=np.histogram2d(wrap(dd[GLON].values),dd[GLAT].values,bins=[lonb,latb],weights=dd.tanom_advective.values); out[bnd]=(g.T)/len(dd)
    return out
mm={}
for r in tz.itertuples():
    m=member_maps(r.initdir,f"{int(r.member):02d}")
    if m is not None: mm[(r.initdir,f"{int(r.member):02d}")]=m
tz=tz[tz.apply(lambda r:(r.initdir,f"{int(r.member):02d}") in mm,axis=1)].reset_index(drop=True)
def comp(rows,bnd): return np.nanmean([mm[(r.initdir,f"{int(r.member):02d}")][bnd] for r in rows.itertuples() if mm[(r.initdir,f"{int(r.member):02d}")][bnd] is not None],axis=0)
tt=tz.nlargest(NTOP,"t2m_C"); tzz=tz.nlargest(NTOP,"z500_gh")
A=comp(tt,BAND_LO)-comp(tzz,BAND_LO); B=comp(tt,BAND_HI)-comp(tzz,BAND_HI); C=comp(tt,BAND_HI); D=comp(tzz,BAND_HI)
print(f"IFS members: {len(tz)}")

# ================= ANALOGUE (best-25 cut-off, 250-750) =================
# Best-25 cut-off analogues (ranks 01..25, same set as Fig 4).  Bottom row = the
# top-NTOP vs lowest-NTOP analogues by cut-off (COL) activity within the best-25.
# COL activity is reconstructed from the catalogue exactly as Fig 4 (5-day rolling Sum(So) of
# closed cut-offs ex==1 & So>=6).  Per-analogue /n normalization then equal-weight mean over
# the selected analogues = the Figure-4 origin-map normalization.
import re
ANA_DIR=os.environ.get("ANA_DIR","${DERIVED_DATA}/analogue_origin_6h_windfix"); N_BEST=25; DAYS_BEFORE=1   # windfix + new genesis, ranks 01..25
CATALOGUE="${DERIVED_DATA}/col_features_AMJJA_1940_2025_largebox.nc"
alonb=np.arange(-180,180+DLON,DLON); alatb=np.arange(0,90+DLON,DLON)
alonc=0.5*(alonb[:-1]+alonb[1:]); alatc=0.5*(alatb[:-1]+alatb[1:]); anlat=alatc.size; anlon=alonc.size
def agrid(path,day):
    ds=xr.open_dataset(path); ti=pd.DatetimeIndex(ds["time"].values).normalize()
    sel=np.where((ti>=day-pd.Timedelta(days=DAYS_BEFORE))&(ti<=day))[0]; dd=ds.isel(time=sel)
    lon=wrap(dd["advection_lon_genesis"].values.ravel()); lat=dd["advection_lat_genesis"].values.ravel()
    val=dd["advection"].values.ravel(); arr=dd["arrival_p"].values.astype(float).ravel(); ds.close()
    ok=np.isfinite(lon)&np.isfinite(lat)&np.isfinite(val)&(arr>=BAND_HI[0])&(arr<=BAND_HI[1]); n=int(ok.sum())
    if n==0: return None
    g,_,_=np.histogram2d(lon[ok],lat[ok],bins=[alonb,alatb],weights=val[ok]); return (g.T)/n
# --- COL activity per day from the catalogue (5-day Sum(So), ex==1 & So>=6) ---
_ct=xr.open_dataset(CATALOGUE); _cf=_ct[["time","So","ex"]].to_dataframe().reset_index(); _ct.close()
_cc=_cf[(_cf.ex==1)&(_cf.So>=6)].copy(); _cc["day"]=pd.DatetimeIndex(_cc.time).normalize()
_dd=_cc.groupby("day")["So"].sum(); _dd=_dd.reindex(pd.date_range(_dd.index.min(),_dd.index.max(),freq="D"),fill_value=0.0)
COL_ACT=_dd.rolling(5,min_periods=1).sum()
def col_of(d):
    ts=pd.Timestamp(d).normalize(); return float(COL_ACT.loc[ts]) if ts in COL_ACT.index else 0.0
folders=[]
for f in glob.glob(ANA_DIR+"/rank*_MJ_*"):
    mm=re.search(r"rank(\d+)_MJ_(\d{8})",os.path.basename(f))
    if mm and os.path.exists(f+"/origin_contributions.nc"):
        folders.append((int(mm.group(1)),mm.group(2),f+"/origin_contributions.nc"))
folders=sorted(folders)[:N_BEST]
grids=[]; dates_ok=[]; cols=[]
for rk,ymd,path in folders:
    g=agrid(path,pd.Timestamp(ymd).normalize())
    if g is None: print("empty",ymd); continue
    grids.append(g); dates_ok.append(ymd); cols.append(col_of(ymd))
grids=np.array(grids); cols=np.array(cols); order=np.argsort(cols)
top=order[::-1][:NTOP]; low=order[:NTOP]
E=np.nanmean(grids[top],0); F=np.nanmean(grids[low],0)
print(f"analogue best-{N_BEST}: top-{NTOP} COL {[dates_ok[i] for i in top]} = {cols[top].round(1)}")
print(f"                        lowest-{NTOP} COL {[dates_ok[i] for i in low]} = {cols[low].round(1)}")

# ================= FIGURE 3x2 =================
proj=ccrs.PlateCarree(); fig,axes=plt.subplots(3,2,figsize=(14,12),subplot_kw={"projection":proj},constrained_layout=True)
def boxln(ax): lo0,lo1,la0,la1=MOR; ax.plot([lo0,lo1,lo1,lo0,lo0],[la0,la0,la1,la1,la0],"k--",lw=1.6,transform=proj)
def box_val(G,LC,LTC):
    Xg,Yg=np.meshgrid(LC,LTC); m=(Xg>=MOR[0])&(Xg<=MOR[1])&(Yg>=MOR[2])&(Yg<=MOR[3]); return float(np.nansum(np.where(np.isfinite(G),G,0.0)[m]))
def pct(*a):                                   # vmax from SIGNAL cells only (drop near-zero background)
    x=np.concatenate([np.abs(g[np.isfinite(g)]) for g in a])
    if x.size==0: return 1.0
    x=x[x>0.02*x.max()]
    return float(np.nanpercentile(x,99.0)) if x.size else 1.0
def area_density(F,LTC):   # per 0.5deg source cell -> per 10^6 km^2 (grid-independent; box integral in K preserved)
    R=6371.0; coslat=np.cos(np.deg2rad(np.asarray(LTC,float)))[:,None]
    return F/((np.deg2rad(DLON)*R)*(np.deg2rad(DLON)*R*coslat))*1e6
def draw(ax,LC,LTC,G,v,title,lett):
    Gd=area_density(G,LTC)   # plot as density; box_val below stays native K
    ax.set_extent(EXT,crs=proj); pm=ax.pcolormesh(LC,LTC,np.where(np.abs(G)>1e-12,Gd,np.nan),cmap=TEMPDIV,vmin=-v,vmax=v,shading="auto",transform=proj)
    boxln(ax); ax.add_feature(cfeature.COASTLINE,lw=0.6); ax.add_feature(cfeature.BORDERS,lw=0.3)
    gl=ax.gridlines(draw_labels=True,ls="--",alpha=0.3); gl.top_labels=gl.right_labels=False; gl.xlabel_style=gl.ylabel_style={"size":12}
    ax.set_title(title.split(') ',1)[1] if title.startswith('(') else title,fontsize=15)   # drop duplicate (x) prefix
    ax.text(-0.02,1.04,f"({lett})",transform=ax.transAxes,fontsize=20,fontweight="bold",va="bottom",ha="left")   # letter OUTSIDE panel
    ax.text(0.97,0.05,f"box = {box_val(G,LC,LTC):+.2f} K",transform=ax.transAxes,ha="right",va="bottom",fontsize=14,fontweight="bold",bbox=dict(boxstyle="square,pad=0.3",fc="white",ec="0.5",alpha=0.9))
    return pm
v1=pct(area_density(A,latc),area_density(B,latc)); v2=pct(area_density(C,latc),area_density(D,latc)); v3=pct(area_density(E,alatc),area_density(F,alatc))
p1=draw(axes[0,0],lonc,latc,A,v1,f"(a) IFS  top-{NTOP} T2m − top-{NTOP} Z500,  850–1000 hPa","a")
draw(axes[0,1],lonc,latc,B,v1,f"(b) IFS  top-{NTOP} T2m − top-{NTOP} Z500,  250–750 hPa","b")
p2=draw(axes[1,0],lonc,latc,C,v2,f"(c) IFS  top-{NTOP} T2m origin,  250–750 hPa","c")
draw(axes[1,1],lonc,latc,D,v2,f"(d) IFS  top-{NTOP} Z500 origin,  250–750 hPa","d")
p3=draw(axes[2,0],alonc,alatc,E,v3,f"Analogue top-{NTOP} COL activity,  250–750 hPa","e")
draw(axes[2,1],alonc,alatc,F,v3,f"Analogue lowest-{NTOP} COL activity,  250–750 hPa","f")
for p,row,lab in [(p1,0,"adv-origin difference [K per 10$^6$ km$^2$]"),(p2,1,"advective origin [K per 10$^6$ km$^2$]"),(p3,2,"advective origin [K per 10$^6$ km$^2$]")]:
    cb=fig.colorbar(p,ax=list(axes[row,:]),location="right",shrink=0.9,aspect=22,pad=0.015); cb.set_label(lab,fontsize=15); cb.ax.tick_params(labelsize=12)
fig.savefig("${HOME_WORK}/supplement_advective_origin_2x3_windfix_ifs20261019" + ("" if IFS_GENESIS=="adv" else "_IFSgenTA") + ".png",dpi=150,bbox_inches="tight")
fig.savefig("${HOME_WORK}/supplement_advective_origin_2x3_windfix_ifs20261019" + ("" if IFS_GENESIS=="adv" else "_IFSgenTA") + ".pdf",bbox_inches="tight"); print("saved png/pdf")
