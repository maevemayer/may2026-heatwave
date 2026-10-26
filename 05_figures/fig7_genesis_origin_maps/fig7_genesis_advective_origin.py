import numpy as np, pandas as pd, xarray as xr
from pathlib import Path
from scipy import stats
from scipy.stats import rankdata
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
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

from matplotlib.collections import LineCollection
import cartopy.crs as ccrs, cartopy.feature as cfeature
import cmocean, os, matplotlib.colors as mcolors
def _ipcc(name,fallback):   # same loader as Fig 1 (Robin's load_colormap): IPCC file -> local copy -> matplotlib fallback
    for p in ["${HOME}/colormaps_IPCC/continuous_colormaps_rgb_0-1/"+name,
              os.path.join(os.path.dirname(os.path.abspath(__file__)),name)]:
        if os.path.exists(p): return mcolors.LinearSegmentedColormap.from_list(name,np.loadtxt(p))
    return plt.get_cmap(fallback)
TEMPDIV=_ipcc("temp_div.txt","RdBu_r")   # IPCC AR6 temperature-divergent (matches Fig 1)
BETACMAP=TEMPDIV   # IPCC temp_div for the β-fit rows too
# IFS: the colleague's tanom_decomp of 2026-10-19 (corrected advective T')
BASE=Path(os.environ.get("IFS_BASE","${IFS_TANOM_DECOMP}"))
HOURS=["00","06","12","18"]; DATES=["20260524","20260525"]; BAND=(250,750); MOR=(-20.,0.,25.,40.); DLON=0.5; NPERM=2000; EXT=[-130,30,5,75]; NF_IFS=15.0
lon_bins=np.arange(-130,80+DLON,DLON); lat_bins=np.arange(15,80+DLON,DLON); lon_c=0.5*(lon_bins[:-1]+lon_bins[1:]); lat_c=0.5*(lat_bins[:-1]+lat_bins[1:])
nlat,nlon=len(lat_c),len(lon_c); LON,LAT=np.meshgrid(lon_c,lat_c)
def outline(ax,mask,LC,LTC,**kw):
    dl=LC[1]-LC[0]; dla=LTC[1]-LTC[0]; lo=np.append(LC-dl/2,LC[-1]+dl/2); la=np.append(LTC-dla/2,LTC[-1]+dla/2); s=[]
    for j in range(mask.shape[0]):
        for i in range(mask.shape[1]):
            if not mask[j,i]: continue
            if i==0 or not mask[j,i-1]: s.append([(lo[i],la[j]),(lo[i],la[j+1])])
            if i==mask.shape[1]-1 or not mask[j,i+1]: s.append([(lo[i+1],la[j]),(lo[i+1],la[j+1])])
            if j==0 or not mask[j-1,i]: s.append([(lo[i],la[j]),(lo[i+1],la[j])])
            if j==mask.shape[0]-1 or not mask[j+1,i]: s.append([(lo[i],la[j+1]),(lo[i+1],la[j+1])])
    ax.add_collection(LineCollection(s,transform=ccrs.PlateCarree(),**kw))
# IFS origin: advective genesis (as ERA5 and the analogues); IFS_GENESIS=ta for the TA genesis
IFS_GENESIS=os.environ.get("IFS_GENESIS","adv")
GLON,GLAT=("lon_adv_genesis","lat_adv_genesis") if IFS_GENESIS=="adv" else ("lon_genesis","lat_genesis")
NEED=[GLON,GLAT,"p0","tanom_advective"]
csv=pd.read_csv("${DERIVED_DATA}/ifs_COL_window_sweep_DAILYMEAN-5.csv"); _tz=pd.read_csv("${DERIVED_DATA}/ifs_cutoff_activity_with_T2m_Z500.csv")[["init","member","t2m_C","z500_gh"]]
csv=csv.merge(_tz,on=["init","member"]); csv["initdir"]=csv["init"].str.replace("-","")+"_00"
def mgrid(initd,mem):
    fr=[pd.read_csv(BASE/initd/mem/f"tanom_{Dd}_{h}",sep=r"\s+") for Dd in DATES for h in HOURS if (BASE/initd/mem/f"tanom_{Dd}_{h}").is_file()]
    if not fr: return None
    d=pd.concat(fr,ignore_index=True)
    for c in NEED: d[c]=pd.to_numeric(d[c],errors="coerce")
    d=d.dropna(subset=NEED); f5=d[(d.p0>=BAND[0])&(d.p0<=BAND[1])]
    if len(f5)==0: return None
    lon=((f5[GLON]+180)%360)-180
    vs,_,_=np.histogram2d(lon,f5[GLAT].values,bins=[lon_bins,lat_bins],weights=f5.tanom_advective.values); return (vs.T)/len(f5)
S=[];COL=[];T=[];Z=[];initlab=[]
for r in csv.itertuples():
    g=mgrid(r.initdir,f"{int(r.member):02d}")
    if g is None: continue
    S.append(g.ravel());COL.append(r.COL_d4/NF_IFS);T.append(r.t2m_C);Z.append(r.z500_gh);initlab.append(r.initdir)   # COL normalised by NF_IFS=15 -> daily-mean COL activity
S=np.asarray(S);COL=np.asarray(COL,float);T=np.asarray(T,float);Z=np.asarray(Z,float);initlab=np.asarray(initlab);nn=len(S)
occ=(np.abs(S>0).sum(0)>=8); inbox=((LON>=MOR[0])&(LON<=MOR[1])&(LAT>=MOR[2])&(LAT<=MOR[3])).ravel(); bc=occ&inbox
idxs=[np.where(initlab==i)[0] for i in sorted(set(initlab))]; rng=np.random.default_rng(0)
def perm():
    o=np.arange(nn)
    for ix in idxs: o[ix]=rng.permutation(ix)
    return o
Xc=S-S.mean(0); vCOL=((COL-COL.mean())**2).sum()/(nn-1)
def bA(c): cc=c-c.mean(); return (Xc*cc[:,None]).sum(0)/(nn-1)/vCOL
gA=bA(COL); mn=np.array([np.nanmax(np.where(occ,np.abs(bA(COL[perm()])),np.nan)) for _ in range(NPERM)]); sigA=np.where(occ,np.abs(gA)>np.quantile(mn,.95),False)
import json as _json
RIDGE_BETA=_json.load(open("${DERIVED_DATA}/ridge_control_t2m_z500_MJJA.json"))["beta_K_per_gpm"]  # learned MJJA Z500->T2m [K/gpm]
A=np.column_stack([np.ones(nn),Z]);P=A@np.linalg.pinv(A); Xr=S-P@S           # branch field: within-sample z500 control
Tr=T-RIDGE_BETA*Z; Tr=Tr-Tr.mean(); Xrc=Xr-Xr.mean(0)                        # T2m ridge control = learned MJJA slope
Bb=Xr[:,bc].sum(1); vB=((Bb-Bb.mean())**2).sum()/(nn-1)
def bB(y): yc=y-y.mean(); return (Xrc*yc[:,None]).sum(0)/(nn-1)/vB
gB=bB(Tr); mn=np.array([np.nanmax(np.where(occ,np.abs(bB(Tr[perm()])),np.nan)) for _ in range(NPERM)]); sigB=np.where(occ,np.abs(gB)>np.quantile(mn,.95),False)
brIFS=S[:,bc].sum(1); rC,pC=stats.pearsonr(COL,brIFS); sC=stats.spearmanr(COL,brIFS)[0]
def resid(y,x): Am=np.column_stack([np.ones(len(y)),x]); return y-Am@np.linalg.lstsq(Am,y,rcond=None)[0]
_TrL=T-RIDGE_BETA*Z; _TrL=_TrL-_TrL.mean()   # T2m ridge-removed via learned MJJA slope (branch stays within-sample)
rD,pD=stats.pearsonr(resid(brIFS,Z),_TrL); sD=stats.spearmanr(rankdata(resid(brIFS,Z)),rankdata(_TrL))[0]
betaC=np.nansum(gA[bc]); betaD=np.nansum(gB[bc])
# ERA5 branch (24-25)
# ERA5 event: windfix + genesis as in Hotz et al. (2024); origin = advection genesis (as the analogues)
d5=xr.open_dataset(os.environ.get("POOLED_FILE","${DERIVED_DATA}/origin_contributions_6h_clim2017_2025_windfix/origin_contributions.nc"))
adv=d5["advection"].values; e_lon=((d5["advection_lon_genesis"].values+180)%360)-180; e_lat=d5["advection_lat_genesis"].values; e_ap=d5["arrival_p"].values; t=pd.DatetimeIndex(d5["time"].values)
k=np.isin(t.strftime("%Y-%m-%d %H"),[f"2026-05-{dd} {h}" for dd in ["24","25"] for h in ["00","06","12","18"]])
aa=adv[k].ravel(); lo=e_lon[k].ravel(); laa=e_lat[k].ravel(); ap=e_ap[k].ravel(); ok=np.isfinite(aa)&(ap>=BAND[0])&(ap<=BAND[1]); nb=ok.sum()
gE,_,_=np.histogram2d(lo[ok],laa[ok],bins=[lon_bins,lat_bins],weights=aa[ok]); gE=(gE.T)/nb
# ============================ ANALOGUE side — everything reconstructed from the input files ============================
# The ONLY hardcoding is the analogue day lists. Per-date drivers are rebuilt at runtime:
#   T2M_ANOM, Z500_IDX  <- present-climate anomaly of the target-region daily means (netCDF)
#   COL_ACT (d4)        <- 5-day Sum(So) of closed cut-offs from the COL catalogue
#   advective grids (X) <- the 6-hourly pooled_origin netCDFs
import glob as _glob
AN_FOLDER    = os.environ.get("ANA_DIR", "${DERIVED_DATA}/analogue_origin_6h_windfix")   # windfix + new genesis analogue origin files (per rank)
AN_MEANS     = "${DERIVED_DATA}/target_region_daily_means_1940_2026.nc"
AN_CATALOGUE = "${DERIVED_DATA}/col_features_AMJJA_1940_2025_largebox.nc"
AN_BAND=BAND; AN_DAYS_BEFORE=1; AN_DLON=AN_DLAT=0.5; AN_NF=15.0; AN_MINCOVER=0.30   # (no longer used: see AN_MIN_N)  NF=15 -> daily-mean COL activity (matches NF_IFS so panels c/d share units)

# --- analogue set: best-25 CUT-OFF analogues from the NEW pipeline (ranks 01..25) ---
import re as _re
N_BEST_ANA=int(os.environ.get("NANA","25"))
_folders=[]
for _d in _glob.glob(AN_FOLDER+"/rank*_MJ_*"):
    _m=_re.search(r"rank(\d+)_MJ_(\d{8})", os.path.basename(_d))
    if _m and os.path.exists(_d+"/origin_contributions.nc"):
        _y=_m.group(2); _folders.append((int(_m.group(1)), f"{_y[:4]}-{_y[4:6]}-{_y[6:8]}", _d+"/origin_contributions.nc"))
_folders=sorted(_folders)[:N_BEST_ANA]
AN_ALL=[dd for _,dd,_ in _folders]
AN_FILE={dd:fp for _,dd,fp in _folders}                        # date -> new-pipeline origin file
AN_GROUP={dd:"factual" for dd in AN_ALL}                       # all cut-off -> single block
AN_FACTUAL=list(AN_ALL); AN_COUNTER=[]
print(f"[analogue] best-{N_BEST_ANA} cut-off (new pipeline): n={len(AN_FACTUAL)}")

# --- T2M_ANOM & Z500_IDX: present-climate anomaly, IDENTICAL to Fig 1/5 ---
#   target-box t2m = FRACTIONAL comp_target_box_mean of the raw 0.5-deg field (not the precomputed t2m_target),
#   climatology = 1995-2025 +-3 d doy window (fixed_doy_climatology), then RMST-detrend + MJJA projection.
_md=xr.open_dataset(AN_MEANS); _mt=pd.DatetimeIndex(_md.time.values).normalize()
_RAW_T2M="${ERA5_DERIVED}/era5_t2m_daily_NA_1940_2026_FULL_0p5deg.nc"
def _box_mean(da,ln,lt,lo0,lo1,la0,la1,frac):   # frac: fractional comp_target_box_mean; else plain cos-lat
    lon,lat=da[ln],da[lt]
    if frac:
        dl=float(abs(lon.values[1]-lon.values[0])); dla=float(abs(lat.values[1]-lat.values[0]))
        fx=((np.minimum(lon+dl/2,lo1)-np.maximum(lon-dl/2,lo0))/dl).clip(0,1)
        fy=((np.minimum(lat+dla/2,la1)-np.maximum(lat-dla/2,la0))/dla).clip(0,1)
        w=fx*fy*np.cos(np.deg2rad(lat)); v=np.isfinite(da)
        return (da.where(v)*w).sum((ln,lt))/(w*v).sum((ln,lt))
    s=da.sel({ln:slice(lo0,lo1)}); s=s.sel({lt:slice(la0,la1)}) if float(s[lt][0])<float(s[lt][-1]) else s.sel({lt:slice(la1,la0)})
    return s.weighted(np.cos(np.deg2rad(s[lt]))).mean((lt,ln))
# Load the raw 0.5-deg t2m ONCE; derive BOTH the target-box series AND the RMST exactly as Fig 1/5.
_da=xr.open_dataarray(_RAW_T2M).squeeze(drop=True)
_ln="lon" if "lon" in _da.coords else "longitude"; _lt="lat" if "lat" in _da.coords else "latitude"
if float(_da[_ln].max())>180: _da=_da.assign_coords({_ln:(((_da[_ln]+180)%360)-180)}).sortby(_ln)
_da=_da.sortby(_lt)
if float(_da.isel(time=0).mean())>100: _da=_da-273.15
def _reidx(da):
    s=pd.Series(np.asarray(da.load().values,float),index=pd.DatetimeIndex(da["time"].values).normalize())
    return s[~s.index.duplicated()].reindex(_mt)
_tgt=_reidx(_box_mean(_da,_ln,_lt,-5.5,5.0,44.5,55.0,True)).to_numpy()           # target box, fractional (Fig 1/5)
# RMST: European box (35-70N,15W-30E) cos-lat mean, 15-yr centred, 30-yr linear edge extrapolation (load_rmst_full)
_rm=_reidx(_box_mean(_da,_ln,_lt,-15,30,35,70,False)).sort_index()
_win=int(round(15*365.25)); _sm=_rm.rolling(_win,center=True,min_periods=_win).mean(); _val=_sm.dropna(); _nf=int(round(30*365.25))
def _edge(fit,at,av,ts):
    x=(fit.index-fit.index[0]).days.values.astype(float); c=np.polyfit(x,fit.values,1)
    return np.polyval(c,(ts-fit.index[0]).days.values.astype(float))-(np.polyval(c,(at-fit.index[0]).days)-av)
_h=_sm.index[_sm.index<_val.index[0]]; _t=_sm.index[_sm.index>_val.index[-1]]
_sm.loc[_h]=_edge(_val.iloc[:_nf],_val.index[0],_val.iloc[0],_h); _sm.loc[_t]=_edge(_val.iloc[-_nf:],_val.index[-1],_val.iloc[-1],_t)
_RMST=_sm.reindex(_mt); _da.close()
_REV=float(_RMST.reindex([pd.Timestamp("2025-05-25")]).iloc[0])
def _fig5_doy_clim(v):
    bmask=(_mt.year>=1995)&(_mt.year<=2025); bdoy=_mt.dayofyear.values[bmask]; bval=v[bmask]
    clim=np.full(367,np.nan)
    for D in range(1,367):
        mm=(np.abs(((bdoy-D+182)%365)-182)<=3)
        if mm.any(): clim[D]=np.nanmean(bval[mm])
    return clim[_mt.dayofyear.values]
def _present_anom(v):
    v=np.asarray(v,float); raw=v-_fig5_doy_clim(v); rm=_RMST.values
    _s=(_mt.month>=5)&(_mt.month<=8)&(_mt.year<=2025); ok=_s&np.isfinite(raw)&np.isfinite(rm); sl=np.polyfit(rm[ok],raw[ok],1)[0]   # MJJA 1940-2025 (exclude 2026), as Fig 5
    return pd.Series(raw+sl*(_REV-rm),index=_mt)
T2M_ANOM=_present_anom(_tgt); Z500_IDX=_present_anom(_md.z500_target.values); _md.close()

# --- COL_ACT (d4): 5-day Sum(So) of ex==1 So>=6 cut-offs from the catalogue ---
_ct=xr.open_dataset(AN_CATALOGUE); _cf=_ct[["time","So","ex"]].to_dataframe().reset_index(); _ct.close()
_cc=_cf[(_cf.ex==1)&(_cf.So>=6)].copy(); _cc["day"]=pd.DatetimeIndex(_cc.time).normalize()
_dd=_cc.groupby("day")["So"].sum(); _dd=_dd.reindex(pd.date_range(_dd.index.min(),_dd.index.max(),freq="D"),fill_value=0.0)
COL_ACT=_dd.rolling(5,min_periods=1).sum()

def _S(series,d):
    v=series.reindex([pd.Timestamp(d)]); return float(v.iloc[0])

# --- advective-origin grid per analogue (global 0.5deg, arrivals [D-1,D], 250-750 hPa band) ---
alon_bins=np.arange(-180,180+AN_DLON,AN_DLON); alat_bins=np.arange(0,90+AN_DLAT,AN_DLAT)
aLC=0.5*(alon_bins[:-1]+alon_bins[1:]); aLTC=0.5*(alat_bins[:-1]+alat_bins[1:]); aNLAT=aLTC.size; aNLON=aLC.size
_ap=AN_FILE                                                     # date -> new-pipeline origin file
def _agrid(date):
    p=_ap.get(date)
    if p is None: return None
    dsx=xr.open_dataset(p); day=pd.Timestamp(date).normalize(); ti=pd.DatetimeIndex(dsx["time"].values).normalize()
    sel=np.where((ti>=day-pd.Timedelta(days=AN_DAYS_BEFORE))&(ti<=day))[0]; ddx=dsx.isel(time=sel)
    # advective-branch ORIGIN = advection genesis (matches old pooled lon/lat = advection genesis)
    lon=((ddx["advection_lon_genesis"].values.ravel()+180)%360)-180; lat=ddx["advection_lat_genesis"].values.ravel()
    val=ddx["advection"].values.ravel(); arr=ddx["arrival_p"].values.astype(float).ravel(); dsx.close()
    ok=np.isfinite(lon)&np.isfinite(lat)&np.isfinite(val)&np.isfinite(arr)&(arr>=AN_BAND[0])&(arr<=AN_BAND[1]); n=int(ok.sum())
    if n==0: return None
    g,_,_=np.histogram2d(lon[ok],lat[ok],bins=[alon_bins,alat_bins],weights=val[ok]); return (g.T)/n
_ag=[]; _ad=[]
for d in AN_ALL:
    g=_agrid(d)
    if g is not None and np.isfinite(_S(T2M_ANOM,d)) and np.isfinite(_S(Z500_IDX,d)):
        _ag.append(g.ravel()); _ad.append(d)
aX=np.asarray(_ag); agrp=np.array([AN_GROUP[d] for d in _ad]); N_ALL=len(_ad)
aT=np.array([_S(T2M_ANOM,d) for d in _ad]); aZ=np.array([_S(Z500_IDX,d) for d in _ad])
aC=np.array([(_S(COL_ACT,d) if pd.Timestamp(d) in COL_ACT.index else 0.0) for d in _ad])/AN_NF
print(f"[analogue] n={N_ALL} (factual {int((agrp=='factual').sum())}, counterfactual {int((agrp=='counterfactual').sum())}) from {os.path.basename(AN_FOLDER)}")

aLON,aLAT=np.meshgrid(aLC,aLTC)
# AN_MIN_N: minimum number of analogues with a positive contribution for a cell to be drawn in (d)/(f)
# and to count in the box (default 4 of 25; was 30 % = 7, which blanked e.g. the Strait of Gibraltar)
AN_MIN_N=int(os.environ.get("AN_MIN_N", 4)); print(f"[analogue] cell kept if >= {AN_MIN_N} of {N_ALL} analogues are positive there")
a_occ=(np.abs(aX>0).sum(0)>=AN_MIN_N)
ainbox=((aLON>=MOR[0])&(aLON<=MOR[1])&(aLAT>=MOR[2])&(aLAT<=MOR[3])).ravel(); a_bc=a_occ&ainbox
an_occ=a_occ.reshape(aNLAT,aNLON); an_grid_a_all=np.nanmean(aX,0).reshape(aNLAT,aNLON)
_aidx=[np.where(agrp==g)[0] for g in sorted(set(agrp))]; _arng=np.random.default_rng(0)
def _aperm():
    o=np.arange(N_ALL)
    for ix in _aidx: o[ix]=_arng.permutation(ix)
    return o
def _abeta(Xm,y): Xc=Xm-Xm.mean(0); yc=y-y.mean(); return (Xc*yc[:,None]).sum(0)/(yc**2).sum()
def _afwer(Xm,y):
    b=_abeta(Xm,y); mx=np.array([np.nanmax(np.where(a_occ,np.abs(_abeta(Xm,y[_aperm()])),np.nan)) for _ in range(NPERM)])
    return b,np.where(a_occ,np.abs(b)>np.quantile(mx,.95),False)
def _ares(y,x): Am=np.column_stack([np.ones(len(y)),x]); return np.asarray(y,float)-Am@np.linalg.lstsq(Am,np.asarray(y,float),rcond=None)[0]
_ge,_se=_afwer(aX,aC); an_grid_e=np.where(a_occ,_ge,np.nan).reshape(aNLAT,aNLON); an_sig_e=_se.reshape(aNLAT,aNLON)  # (d) COL->branch
abranch=aX[:,a_bc].sum(1)
ann_e=np.array([stats.pearsonr(aC,abranch)[0],stats.pearsonr(aC,abranch)[1],stats.spearmanr(aC,abranch)[0],np.nansum(_ge[a_bc])])
# (f) advective branch -> T2m  DIRECT (no Z500 control; IFS panel (e) stays ridge-controlled)
_gf,_sf=_afwer(aX,aT); an_grid_f=np.where(a_occ,_gf,np.nan).reshape(aNLAT,aNLON); an_sig_f=_sf.reshape(aNLAT,aNLON)
_vB=((abranch-abranch.mean())**2).sum()
_bf=((abranch-abranch.mean())*(aT-aT.mean())).sum()/_vB
ann_f=np.array([stats.pearsonr(abranch,aT)[0],stats.pearsonr(abranch,aT)[1],stats.spearmanr(rankdata(abranch),rankdata(aT))[0],_bf])
print(f"[analogue] (d) COL->branch r={ann_e[0]:+.2f} beta_box={ann_e[3]:+.3f} | (f) branch->T2m DIRECT r={ann_f[0]:+.2f} beta_box={ann_f[3]:+.3f}")
aLON,aLAT=np.meshgrid(aLC,aLTC); inbox_an=(aLON>=MOR[0])&(aLON<=MOR[1])&(aLAT>=MOR[2])&(aLAT<=MOR[3]); bc_an=an_occ&inbox_an
gf=an_grid_f; boxint=np.nansum(gf[bc_an]); grid_f_fwd=gf*(ann_f[3]/boxint if boxint else 1.0)   # rescale to forward
GA2=np.where(occ,gA,np.nan).reshape(nlat,nlon); GB2=np.where(occ,gB,np.nan).reshape(nlat,nlon)
fs = 17
axfs = fs-2   # x/y tick label size (offset from fs, not hardcoded)
legfs = fs-3  # in-panel stats/legend box — unchanged size, just not hardcoded
plt.rcParams.update({"font.size":fs,"axes.titlesize":fs+1,"axes.labelsize":fs,"xtick.labelsize":axfs,"ytick.labelsize":axfs})
proj=ccrs.PlateCarree(); fig=plt.figure(figsize=(15,14.5),constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=0.02,h_pad=0.02,wspace=0.01,hspace=0.035)   # a bit more row spacing
def base(pos):
    ax=fig.add_subplot(3,2,pos,projection=proj); ax.set_extent(EXT,crs=proj)
    ax.add_feature(cfeature.LAND,alpha=0.45); ax.add_feature(cfeature.OCEAN,alpha=0.2); ax.add_feature(cfeature.COASTLINE,lw=0.5); ax.add_feature(cfeature.BORDERS,lw=0.25)
    gl=ax.gridlines(draw_labels=True,ls="--",alpha=0.3); gl.top_labels=gl.right_labels=False; gl.xlabel_style=gl.ylabel_style={"size":axfs}; return ax
def drawp(pos,l,LC,LTC,F,sig,title,vmax,ann=None,cmap="RdBu_r",vmin=None):
    ax=base(pos); pm=ax.pcolormesh(LC,LTC,F,cmap=cmap,vmin=(-vmax if vmin is None else vmin),vmax=vmax,transform=proj)
    if sig is not None: outline(ax,sig,LC,LTC,colors="k",linewidths=0.5)
    lo0,lo1,la0,la1=MOR; ax.plot([lo0,lo1,lo1,lo0,lo0],[la0,la0,la1,la1,la0],"k--",lw=1.4,transform=proj)
    ax.text(-0.04,1.03,f"({l})",transform=ax.transAxes,fontsize=fs+3,fontweight="bold",va="bottom",ha="right")   # letter OUTSIDE panel (upper-left)
    ax.set_title(title,fontsize=fs)
    if ann is not None:
        ax.text(0.98,0.03,f"box  Pearson {ann[0]:+.2f} ({'p<0.001' if ann[1]<0.001 else 'p='+format(ann[1],'.3f')})\nSpearman {ann[2]:+.2f}\nβ_box = {ann[3]:+.3f}",transform=ax.transAxes,fontsize=legfs,va="bottom",ha="right",bbox=dict(boxstyle="round",fc="white",alpha=0.85))
    return ax,pm
def rowbar(pms,axs,label): cb=fig.colorbar(pms,ax=axs,location='bottom',shrink=0.6,aspect=40,pad=0.02); cb.set_label(label,fontsize=fs-1); cb.ax.tick_params(labelsize=axfs)
VMAX_PCTL=99.5   # colour vmax = this percentile of |field| over SIGNAL cells (background ~0 excluded)
def pct(*arrs):
    a=np.concatenate([np.abs(np.asarray(x)[np.isfinite(x)]) for x in arrs]); a=a[a>1e-9]
    return float(np.nanpercentile(a,VMAX_PCTL)) if a.size else 1.0
def area_density(F,LATg):   # per 0.5deg source cell -> per 10^6 km^2 (grid-independent; box integral preserved). Both grids are 0.5deg.
    R=6371.0; cell_km2=(np.deg2rad(DLON)*R)*(np.deg2rad(DLON)*R*np.cos(np.deg2rad(LATg)))
    return F/cell_km2*1e6
# Row 1: origin — SEPARATE colorbars (amp, positives only) | col1=ERA5, col2=analogue
grid_a_pos=area_density(np.where(an_grid_a_all>0,an_grid_a_all,np.nan),aLAT)   # entire composite; positives only; K per 10^6 km^2
gE_pos=area_density(np.where(gE>0,gE,np.nan),LAT)                              # positives only; K per 10^6 km^2
vb=pct(gE_pos)        # ERA5 advective-origin scale (signal percentile, positives only)
va=pct(grid_a_pos)    # analogue advective-origin scale (signal percentile)
a1,p1=drawp(1,"a",lon_c,lat_c,gE_pos,None,"ERA5 advective branch (24–25 May)",vb,cmap=cmocean.cm.amp,vmin=0)
cb1=fig.colorbar(p1,ax=a1,location='bottom',shrink=0.75,aspect=25,pad=0.02); cb1.set_label("advective origin [K per 10$^6$ km$^2$]",fontsize=fs-1); cb1.ax.tick_params(labelsize=axfs)
a2,p2=drawp(2,"b",aLC,aLTC,grid_a_pos,None,f"Analogue advective branch (best-{N_ALL} cut-off)",va,cmap=cmocean.cm.amp,vmin=0)
cb2=fig.colorbar(p2,ax=a2,location='bottom',shrink=0.75,aspect=25,pad=0.02); cb2.set_label("advective origin [K per 10$^6$ km$^2$]",fontsize=fs-1); cb2.ax.tick_params(labelsize=axfs)
# --- optional Z500 contour overlay on the ERA5 panel (a): activates when a gridded Z500 file is present ---
# daily-mean ERA5 Z500, 1 deg, 60W-40E x 15-75N, 1940-2026 (2026 checked: no offset); m^2/s^2
Z500_FILE="${ERA5_DERIVED}/era5_z500_daily_NA_1940_2026_FULL_1deg.nc"
# off by default: a 24-25 May snapshot is misleading next to origins of different times
Z500_CONTOURS=os.environ.get("Z500_CONTOURS","0")=="1"
if Z500_CONTOURS and os.path.exists(Z500_FILE):
    zz=xr.open_dataset(Z500_FILE)
    zv=next(v for v in ("z","z500","zg","Z","var129") if v in zz.variables)
    lonn=next(c for c in ("lon","longitude") if c in zz.coords); latn=next(c for c in ("lat","latitude") if c in zz.coords); tn=next(c for c in ("time","valid_time") if c in zz.coords)
    zt=pd.DatetimeIndex(zz[tn].values); ksel=np.isin(zt.strftime("%Y-%m-%d"),["2026-05-24","2026-05-25"])
    Zf=zz[zv].isel({tn:np.where(ksel)[0]}).mean(tn) if ksel.any() else zz[zv].mean(tn)
    zlon=((zz[lonn].values+180)%360)-180; order=np.argsort(zlon); zlat=zz[latn].values
    Zf=Zf.squeeze(drop=True)                                          # drop the size-1 plev axis
    Z=np.asarray(Zf.values); Z=Z/9.80665 if np.nanmean(Z)>2e4 else Z   # m^2/s^2 -> gpm
    cs=a1.contour(zlon[order],zlat,Z[:,order],levels=np.arange(5400,6000,60),colors="k",linewidths=0.7,transform=proj)
    a1.clabel(cs,fmt="%d",fontsize=8)
    print("overlaid Z500 contours on panel (a)")
else:
    print("[note] panel (a) without Z500 contours" + ("" if Z500_CONTOURS else " (Z500_CONTOURS=1 to draw them)"))
# Row 2: COL -> branch   (IFS | analogue) — SEPARATE per-panel colour scales
GA2d=area_density(GA2,LAT); an_grid_e_d=area_density(an_grid_e,aLAT)
a3,p3=drawp(3,"c",lon_c,lat_c,GA2d,sigA.reshape(nlat,nlon),"IFS: cut-off activity → advective branch",pct(GA2d),ann=[rC,pC,sC,betaC],cmap=BETACMAP)
cb3=fig.colorbar(p3,ax=a3,location="bottom",shrink=0.85,aspect=25,pad=0.02); cb3.set_label("K adv / mean daily COL activity / 10$^6$ km$^2$",fontsize=fs-2)
a4,p4=drawp(4,"d",aLC,aLTC,an_grid_e_d,an_sig_e,"Analogue: cut-off activity → advective branch",pct(an_grid_e_d),ann=ann_e,cmap=BETACMAP)
cb4=fig.colorbar(p4,ax=a4,location="bottom",shrink=0.85,aspect=25,pad=0.02); cb4.set_label("K adv / mean daily COL activity / 10$^6$ km$^2$",fontsize=fs-2)
# Row 3: branch -> T2m   (IFS | Z500  ;  analogue DIRECT) — SEPARATE per-panel colour scales
GB2d=area_density(GB2,LAT); grid_f_fwd_d=area_density(grid_f_fwd,aLAT)
a5,p5=drawp(5,"e",lon_c,lat_c,GB2d,sigB.reshape(nlat,nlon),"IFS: advective branch → T2m | Z500",pct(GB2d),ann=[rD,pD,sD,betaD],cmap=BETACMAP)
cb5=fig.colorbar(p5,ax=a5,location="bottom",shrink=0.85,aspect=25,pad=0.02); cb5.set_label("K T2m / K branch / 10$^6$ km$^2$",fontsize=fs-2)
a6,p6=drawp(6,"f",aLC,aLTC,grid_f_fwd_d,an_sig_f,"Analogue: advective branch → T2m  (direct)",pct(grid_f_fwd_d),ann=ann_f,cmap=BETACMAP)
cb6=fig.colorbar(p6,ax=a6,location="bottom",shrink=0.85,aspect=25,pad=0.02); cb6.set_label("K T2m / K branch / 10$^6$ km$^2$",fontsize=fs-2)
OUT=f"${HOME_WORK}/Fig4_6panel_250_750hPa_0p5deg_6h_cutoff_N{N_BEST_ANA}_windfix_ifs20261019" + ("" if IFS_GENESIS=="adv" else "_IFSgenTA") + (f"_minN{AN_MIN_N}" if "AN_MIN_N" in os.environ else "") + ".png"; fig.savefig(OUT,dpi=150,bbox_inches="tight"); fig.savefig(OUT.replace(".png",".pdf"),bbox_inches="tight"); print("saved",OUT)
print(f"IFS  COL->branch beta_box={betaC:+.3f} | branch->T2m beta_box={betaD:+.2f}")
print(f"ANA  COL->branch beta_box={ann_e[3]:+.3f} | branch->T2m beta_box={ann_f[3]:+.2f}")
