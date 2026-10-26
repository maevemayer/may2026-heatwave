import numpy as np, pandas as pd, xarray as xr, os
from scipy import stats
from scipy.interpolate import griddata
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.dates as mdates, matplotlib.colors as mcolors
import cartopy.crs as ccrs, cartopy.feature as cfeature
FD="${DERIVED_DATA}/aifs"
DL="${DERIVED_DATA}/"
fs = 13
axfs = fs-2   # x/y tick label size (offset from fs, not hardcoded)
red,blue="#d62728","#1f77b4"; proj=ccrs.PlateCarree(); plt.rcParams.update({"font.size":fs,"hatch.linewidth":0.35,"xtick.labelsize":axfs,"ytick.labelsize":axfs})
def _ipcc(name,fb):
    for p in ["${HOME_WORK}/"+name,
              "/path/to/colormaps/"+name]:
        if os.path.exists(p): return mcolors.LinearSegmentedColormap.from_list(name,np.loadtxt(p))
    return plt.get_cmap(fb)
TEMPDIV=_ipcc("temp_div.txt","RdBu_r")

def fit_sn(v):
    v=np.asarray(v,float); mu,s=v.mean(),v.std(); raw=np.mean((v-mu)**3)
    sk=np.clip(abs(raw/s**3),1e-3,0.99)*np.sign(raw)
    d=np.sqrt(np.pi/2*abs(sk)**(2/3)/(abs(sk)**(2/3)+((4-np.pi)/2)**(2/3)))*np.sign(sk)
    a=d/np.sqrt(1-d**2); sc=s/np.sqrt(max(1-2*d**2/np.pi,1e-6)); return a,mu-sc*d*np.sqrt(2/np.pi),sc
def lbox(ax,letter):
    ax.text(0.0,1.1,f"({letter})",transform=ax.transAxes,fontsize=fs+2,fontweight="bold",va="bottom",ha="right",zorder=12)   # letter OUTSIDE panel
def tag_txt(ax,txt,y=0.85):
    ax.text(0.035,y,txt,transform=ax.transAxes,va="top",ha="left",fontsize=fs,zorder=11,
            bbox=dict(boxstyle="round,pad=0.2",fc="white",ec="none",alpha=0.7))
def pdf_panel(ax,cut,noc,letter,xlab,tag,legloc="upper right",vline=None,vlabel=None):
    xs=np.linspace(min(cut.min(),noc.min())-1,max(cut.max(),noc.max())+1,300)
    for g,c,nm in [(cut,red,"COL"),(noc,blue,"no COL")]:
        ax.hist(g,bins=12,density=True,alpha=0.30,color=c)
        sh,lo,sc=fit_sn(g); ax.plot(xs,stats.skewnorm.pdf(xs,sh,loc=lo,scale=sc),color=c,lw=2.3,label=f"{nm} (n={len(g)})")
        ax.axvline(g.mean(),color=c,ls="--",lw=1.1)
    if vline is not None:
        ax.axvline(vline,color="k",ls="--",lw=1.8,zorder=6,label=vlabel)   # realized ERA5 2026 anomaly
    d=cut.mean()-noc.mean(); wp=stats.ttest_ind(cut,noc,equal_var=False).pvalue; ks=stats.ks_2samp(cut,noc).pvalue
    ax.set_xlabel(xlab); ax.set_ylabel("Density"); ax.grid(alpha=.3); ax.legend(fontsize=fs-3,loc=legloc,handlelength=0.9,handletextpad=0.5)
    ax.set_title(f"{tag},  Δ={d:+.2f} K",fontsize=fs); lbox(ax,letter)

fig=plt.figure(figsize=(12,15)); gs=fig.add_gridspec(3,2,hspace=0.36,wspace=0.22,height_ratios=[1,1.12,1.12])

# ---------- (a) IFS PDF: ex==1, Z500-controlled ----------
tz=pd.read_csv(f"{DL}/ifs_box_daily_t2m_z500_20260525.csv")[["init","member","t2m_C","z500_m"]]   # new fractional box-mean (0.5deg t2m, 1deg z500)
cls=pd.read_csv(f"{DL}/ifs_COL_window_sweep_DAILYMEAN-5.csv")[["init","member","cutoff25"]]
mi=cls.merge(tz,on=["init","member"]); Ti=mi.t2m_C.values; Zi=mi.z500_m.values
# Z500 control in the ERA5-CLIMATOLOGICAL frame: the SAME learned slope (beta_anom) AND the SAME
# intercept (alpha_era5) are used for the IFS members, the analogues, and the ERA5 event -> one physical
# zero ("what ERA5 climatology expects for this ridge"), so the panels and the event line are comparable.
# res_i = (Ti - Tclim) - alpha - beta*(Zi - Zclim):  present-climate T2m anomaly minus the ERA5 ridge fit.
import json as _json
_j=_json.load(open("${DERIVED_DATA}/ridge_control_t2m_z500_MJJA.json"))
RIDGE_BETA=_j["beta_anom_K_per_gpm"]; ALPHA=_j["intercept_anom_K"]
TCLIM_E=_j["t2m_clim_doy_event_C"]; ZCLIM_E=_j["z500_clim_doy_event_m"]   # 25-May 1995-2025 doy clims (frac box mean)
res=(Ti-TCLIM_E)-ALPHA-RIDGE_BETA*(Zi-ZCLIM_E)   # Z500-controlled residual, ERA5-climatology zero (no ensemble re-centering)
g=mi.cutoff25.astype(bool).values
# panel (a) is drawn below, once OBS + Z_OBS are available for the ERA5 2026 realized residual line

# ---------- (b) AIFS PDF (binary wide, ex==1) -> T2m ANOMALY ----------
# Anomaly = T2m - target-box 25-May climatology (area-weighted daily-mean t2m, UK-France box,
# doy 145, 1995-2025, +/-3 d).  Computed from the ERA5 0.5deg daily file if present, else 13.77 C.
# (Constant offset: cutoff/no-cutoff Delta and PDF shapes unchanged, only the axis shifts.)
CLIM_T2M_FILE="${ERA5_DERIVED}/era5_t2m_daily_NA_1940_2026_FULL_0p5deg.nc"
def target_box_clim_obs(box=(-5.5,5.0,44.5,55.0),doy=145,yr=(1995,2025),half=3,fb_clim=13.77,fb_obs=21.15):
    # Uses Figure 1's FRACTIONAL cos-lat box mean (comp_target_box_mean): each edge cell is
    # weighted by the fraction of its area inside the box -> grid-independent, so OBS/CLIM and
    # hence the realized anomaly (+7.39 C) match Figure 1 exactly (a plain .sel slice gave +7.23).
    try:
        da=xr.open_dataarray(CLIM_T2M_FILE).squeeze(drop=True)
        ln="lon" if "lon" in da.coords else "longitude"; lt="lat" if "lat" in da.coords else "latitude"
        if float(da[ln].max())>180: da=da.assign_coords({ln:(((da[ln]+180)%360)-180)}).sortby(ln)
        da=da.sortby(lt)
        if float(da.isel(time=0).mean())>100: da=da-273.15
        lon,lat=da[ln],da[lt]; dlon=float(abs(lon.values[1]-lon.values[0])); dlat=float(abs(lat.values[1]-lat.values[0]))
        fx=((np.minimum(lon+dlon/2,box[1])-np.maximum(lon-dlon/2,box[0]))/dlon).clip(0,1)
        fy=((np.minimum(lat+dlat/2,box[3])-np.maximum(lat-dlat/2,box[2]))/dlat).clip(0,1)
        W=fx*fy*np.cos(np.deg2rad(lat)); valid=np.isfinite(da)
        b=((da.where(valid)*W).sum((ln,lt))/(W*valid).sum((ln,lt))).load()
        t=pd.DatetimeIndex(b["time"].values); dist=np.abs(((t.dayofyear.values-doy+182)%365)-182)
        m=(dist<=half)&(t.year>=yr[0])&(t.year<=yr[1])
        return float(b.values[m].mean()), float(b.sel(time="2026-05-25",method="nearest").values)
    except Exception as e:
        print("[clim] using fallback",fb_clim,fb_obs,"(",e,")"); return fb_clim,fb_obs
Z500_OBS_FILE="${ERA5_DERIVED}/era5_z500_daily_NA_1940_2026_FULL_1deg.nc"
def target_box_z500_obs(box=(-5.5,5.0,44.5,55.0),fb=5916.9):
    # ERA5 realized 25-May-2026 Z500 box mean (1deg, fractional cos-lat; SAME box/method as the IFS z500_m)
    try:
        da=xr.open_dataarray(Z500_OBS_FILE).squeeze(drop=True)
        ln="lon" if "lon" in da.coords else "longitude"; lt="lat" if "lat" in da.coords else "latitude"
        if float(da[ln].max())>180: da=da.assign_coords({ln:(((da[ln]+180)%360)-180)}).sortby(ln)
        da=da.sortby(lt)
        lon,lat=da[ln],da[lt]; dlon=float(abs(lon.values[1]-lon.values[0])); dlat=float(abs(lat.values[1]-lat.values[0]))
        fx=((np.minimum(lon+dlon/2,box[1])-np.maximum(lon-dlon/2,box[0]))/dlon).clip(0,1)
        fy=((np.minimum(lat+dlat/2,box[3])-np.maximum(lat-dlat/2,box[2]))/dlat).clip(0,1)
        W=fx*fy*np.cos(np.deg2rad(lat)); valid=np.isfinite(da)
        b=((da.where(valid)*W).sum((ln,lt))/(W*valid).sum((ln,lt))).load()
        z=float(b.sel(time="2026-05-25",method="nearest").values); return z/9.80665 if z>2e4 else z
    except Exception as e:
        print("[z500obs] using fallback",fb,"(",e,")"); return fb
CLIM,OBS=target_box_clim_obs(); OBS_ANOM=OBS-CLIM
print(f"target-box 25-May climatology = {CLIM:.2f} C | observed = {OBS:.2f} C | realized anomaly = {OBS_ANOM:+.2f} C")
# ---------- (a) IFS PDF: Z500-controlled residual, with ERA5 2026 realized residual line ----------
Z_OBS=target_box_z500_obs()
res_era5=(OBS-TCLIM_E)-ALPHA-RIDGE_BETA*(Z_OBS-ZCLIM_E)   # ERA5 2026 event on the SAME ERA5-climatological zero
print(f"ERA5 2026 Z500-controlled residual = {res_era5:+.2f} K (Z_OBS={Z_OBS:.1f} m) | IFS ens-mean (cold bias vs ERA5) = {res.mean():+.2f} K, max member = {res.max():+.2f} K")
pdf_panel(fig.add_subplot(gs[0,0]),res[g],res[~g],"a","T2m residual | $Z_{500}$  [K]","IFS $Z_{500}$-controlled",
          legloc="upper left",vline=res_era5,vlabel=f"ERA5 2026 ({res_era5:+.1f} K)")
pdf=pd.read_csv(f"{FD}/pdf_data_with_threshold.csv")
ac=pdf["T2m_With_Cutoff"].dropna().values-CLIM; an=pdf["T2m_Without_Cutoff"].dropna().values-CLIM
pdf_panel(fig.add_subplot(gs[0,1]),ac,an,"b","2 m temperature anomaly  [°C]","AIFS $Z_{500}$-matched",legloc="upper left",
          vline=OBS_ANOM,vlabel=f"ERA5 2026 ({OBS_ANOM:+.1f} K)")

# ---------- significance maps (c,e) from precomputed nc + Z500 contours ----------
sig=xr.open_dataset(f"{FD}/significance_cutoff_vs_no_cutoff.nc"); LO=sig["longitude"].values; LA=sig["latitude"].values
gx=np.arange(np.floor(LO.min()),np.ceil(LO.max())+0.4,0.4); gy=np.arange(np.floor(LA.min()),np.ceil(LA.max())+0.4,0.4); GX,GY=np.meshgrid(gx,gy)
EXT=[-25,25,30,65]
amax=float(max(np.abs(sig.anomaly_17th_C).max(),np.abs(sig.anomaly_18th_C).max())); amax=np.ceil(amax*2)/2
lv=np.arange(-amax,amax+0.25,0.5)
def map_panel(pos,letter,ak,pk,zk,lab):
    ax=fig.add_subplot(pos,projection=proj); ax.set_extent(EXT,crs=proj)
    Amap=griddata((LO,LA),sig[ak].values,(GX,GY),method="linear")
    Pmap=griddata((LO,LA),sig[pk].values,(GX,GY),method="linear")
    Zmap=griddata((LO,LA),sig[zk].values,(GX,GY),method="linear")   # already gpm
    pm=ax.contourf(GX,GY,Amap,levels=lv,cmap=TEMPDIV,extend="both",transform=proj)
    ax.contourf(GX,GY,Pmap,levels=[0.05,1.0],hatches=[".."],colors="none",transform=proj)     # inverted: light stipple NON-significant p>0.05
    cs=ax.contour(GX,GY,Zmap,levels=np.arange(5400,5920,40),colors="k",linewidths=1.0,transform=proj)
    ax.clabel(cs,fmt="%d",fontsize=8)
    ax.add_feature(cfeature.COASTLINE,lw=0.6); ax.add_feature(cfeature.BORDERS,lw=0.3)
    gl=ax.gridlines(draw_labels=True,ls="--",alpha=0.25); gl.top_labels=gl.right_labels=False; gl.xlabel_style=gl.ylabel_style={"size":axfs}
    ax.set_title(f"{lab}, top-10 COL",fontsize=fs); lbox(ax,letter)
    return ax,pm
axc,pmc=map_panel(gs[1,0],"c","anomaly_17th_C","p_values_17th","z500_top10_17th","init 17 May")
axe,pme=map_panel(gs[2,0],"e","anomaly_18th_C","p_values_18th","z500_top10_18th","init 18 May")
cb=fig.colorbar(pmc,ax=[axc,axe],location="bottom",shrink=0.8,aspect=32,pad=0.06); cb.set_label("T2m anomaly (°C)")

# ---------- Hovmollers (d,f): reconstruct top10 vs exceeding, ΔT + significance ----------
GPM=5860.0;   # AIFS subselection threshold (was 57500/9.80665 = 5863.4 gpm)
RNG=np.random.default_rng(0)
def lists(init):
    rm=pd.read_csv(f"{FD}/regional_metrics_Init_{init}_20260525.csv")
    exc=set(rm[rm["Z500_Regional_Mean_Daily_gpm"]>GPM]["Member"])
    ci=pd.read_csv(f"{FD}/cols_index_Init_{init}.csv")
    with_lp=set(ci[(ci.Target_Date=="2026-05-25")&(ci.up_or_down=="up")&(ci.ex==1)&(ci.Member.isin(exc))]["Member"].unique())
    d4=pd.read_csv(f"{FD}/COL_d4_index_Init_{init}_UP_EX1.csv")
    m=(d4.Member.isin(exc))&(d4.Total_So>0)&(d4.Member.isin(with_lp))
    top10=d4[m].sort_values("Total_So",ascending=False)["Member"].unique()[:10]
    return list(top10),list(exc),list(exc-with_lp)   # top10, all-exceeding, no-cutoff-exceeding
X0,X1=mdates.date2num(pd.Timestamp("2026-05-20")),mdates.date2num(pd.Timestamp("2026-05-30"))
def hov_panel(pos,letter,init,lab):
    v=xr.open_dataset(f"{FD}/vertical_temperature_profiles_Init_{init}.nc")
    lev=v["level"].values; vt=pd.DatetimeIndex(v["valid_time"].values); num=v["number"].values
    top10,exc,without=lists(init); ci=np.isin(num,top10); ei=np.isin(num,exc)   # ref = all exceeding
    T=v["temperature"].values                                   # (level, number, valid_time)
    dT=np.nanmean(T[:,ci,:],axis=1)-np.nanmean(T[:,ei,:],axis=1)
    P=v["p_values"].values                                      # precomputed significance (all-exceeding reference)
    tnum=mdates.date2num(vt.to_pydatetime()); ax=fig.add_subplot(pos)
    cf=ax.contourf(tnum,lev,dT,levels=np.arange(-4,4.1,0.5),cmap="RdBu_r",extend="both")
    ax.contourf(tnum,lev,P,levels=[0.05,1.0],hatches=[".."],colors="none")   # inverted: light stipple NON-significant p>0.05
    yt=np.arange(1000,199,-200)
    ax.set_yscale("log"); ax.set_ylim(1000,200); ax.set_yticks(yt); ax.set_yticklabels(yt); ax.minorticks_off()
    ax.set_xlim(X0,X1); ax.axvline(mdates.date2num(pd.Timestamp("2026-05-25")),color="k",ls="--",lw=1)
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2)); ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d")); ax.tick_params(axis="x",labelrotation=0)
    ax.set_ylabel("Pressure (hPa)")
    ax.set_title(f"{lab}, ΔT top-10 COL",fontsize=fs); lbox(ax,letter)
    return ax,cf
axd,cfd=hov_panel(gs[1,1],"d","2026051700","init 17 May")
axf,cff=hov_panel(gs[2,1],"f","2026051800","init 18 May")
cb2=fig.colorbar(cff,ax=[axd,axf],location="bottom",shrink=0.8,aspect=32,pad=0.10); cb2.set_label("ΔT (°C)")

# ---- align rows: let cartopy settle the maps' aspect-correct positions, then snap each hovmoller to its map ----
for _ in range(2):                                     # two passes so the maps' final drawn boxes are stable
    fig.canvas.draw()
    for mapax,hovax in [(axc,axd),(axe,axf)]:
        mb=mapax.get_position(); hb=hovax.get_position()
        hovax.set_position([hb.x0,mb.y0,hb.width,mb.height])   # match hovmoller y-span to the map (maps untouched)

OUT="${HOME_WORK}/Fig2_z500anom.png"
fig.savefig(OUT,dpi=150,bbox_inches="tight")
fig.savefig("${HOME_WORK}/Fig2_z500anom.pdf",bbox_inches="tight"); print("saved",OUT,"+ Desktop PDF (z500-anomaly ridge link)")
for init in ["2026051700","2026051800"]:
    t,e,w=lists(init); print(init,"top10:",t,"| exceeding n=",len(e),"| no-cutoff n=",len(w))
