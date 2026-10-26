"""SUPPLEMENT (new): ERA5 2026 EVENT origin maps at 250-750 hPa for ALL processes + TA.
Panels:  (a) Total anomaly (TA)   (b) Advection   (c) Adiabatic   (d) Diabatic
Arrivals 24-25 May 2026, mid-tropospheric advective branch (250-750 hPa origin band, 0.5 deg).

Normalization = the Figure-4 origin-map normalization: mean process contribution PER in-band
trajectory by genesis location.  A SHARED denominator nb (trajectories finite across ALL four
terms AND inside the arrival band) is used for every panel, so the decomposition is additive:
TA  ~=  advection + adiabatic + diabatic.

Letters are OUTSIDE the panels (paper style).  Saves to Downloads + Desktop.
"""
import numpy as np, pandas as pd, xarray as xr, os
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

import matplotlib.colors as mcolors
import cartopy.crs as ccrs, cartopy.feature as cfeature

def _ipcc(name, fb):
    for p in ["${HOME_WORK}/"+name,
              "/path/to/colormaps/"+name]:
        if os.path.exists(p): return mcolors.LinearSegmentedColormap.from_list(name, np.loadtxt(p))
    return plt.get_cmap(fb)
TEMPDIV = _ipcc("temp_div.txt", "RdBu_r")

BAND = (250, 750); MOR = (-20., 0., 25., 40.); DLON = 0.5; EXT = [-130, 30, 5, 75]
DL = "${DERIVED_DATA}"
lon_bins = np.arange(-130, 80+DLON, DLON); lat_bins = np.arange(15, 80+DLON, DLON)
lon_c = 0.5*(lon_bins[:-1]+lon_bins[1:]); lat_c = 0.5*(lat_bins[:-1]+lat_bins[1:])
def wrap(x): return ((np.asarray(x)+180) % 360)-180
SEL = [f"2026-05-{dd} {h}" for dd in ["24", "25"] for h in ["00", "06", "12", "18"]]

# order requested: TA first, then the three processes
TERMS = [("TA",        "a", "Total anomaly (TA)"),
         ("advection", "b", "Advection"),
         ("adiabatic", "c", "Adiabatic"),
         ("diabatic",  "d", "Diabatic")]

# --- NEW pipeline: single ERA5 event origin file (all terms + per-process genesis) ---
# Each process is binned at ITS OWN genesis location (<proc>_lon/lat_genesis).
# windfix + genesis as in Hotz et al. (2024)
NEWF = os.environ.get("POOLED_FILE", "${DERIVED_DATA}/origin_contributions_6h_clim2017_2025_windfix/origin_contributions.nc")
d = xr.open_dataset(NEWF)
t = pd.DatetimeIndex(d["time"].values); k = np.isin(t.strftime("%Y-%m-%d %H"), SEL)
ap = d["arrival_p"].values[k].ravel()
raw = {var: dict(val=d[var].values[k].ravel(),
                 lon=wrap(d[f"{var}_lon_genesis"].values[k].ravel()),
                 lat=d[f"{var}_lat_genesis"].values[k].ravel()) for var, _, _ in TERMS}
d.close()

# shared denominator: trajectory finite in ALL four terms AND inside the arrival band
finite_all = np.ones_like(ap, bool)
for var, _, _ in TERMS: finite_all &= np.isfinite(raw[var]["val"])
ok = finite_all & (ap >= BAND[0]) & (ap <= BAND[1]); nb = int(ok.sum())
print(f"ERA5 2026 event: {nb} in-band trajectories (shared denominator across all terms)")

maps = {}; boxmean = {}
for var, _, _ in TERMS:
    lo, la = raw[var]["lon"], raw[var]["lat"]                     # THIS process's own genesis
    m = ok & np.isfinite(lo) & np.isfinite(la)
    g, _, _ = np.histogram2d(lo[m], la[m], bins=[lon_bins, lat_bins], weights=raw[var]["val"][m])
    maps[var] = (g.T)/nb                                          # Fig-4 normalization: K per in-band trajectory
    ib = m & (lo >= MOR[0]) & (lo <= MOR[1]) & (la >= MOR[2]) & (la <= MOR[3])
    boxmean[var] = np.nansum(raw[var]["val"][ib])/nb
# per-process genesis -> maps are NOT cell-additive (each panel is its own origin); report box means
print(f"box means [K/traj]: TA={boxmean['TA']:+.3f}  adv={boxmean['advection']:+.3f}  "
      f"adiab={boxmean['adiabatic']:+.3f}  diab={boxmean['diabatic']:+.3f}")

def area_density(F, LTC):   # per 0.5deg source cell -> per 10^6 km^2 (grid-independent; box integral in K preserved)
    R = 6371.0; coslat = np.cos(np.deg2rad(np.asarray(LTC, float)))[:, None]
    return F/((np.deg2rad(DLON)*R)*(np.deg2rad(DLON)*R*coslat))*1e6
dmaps = {v: area_density(maps[v], lat_c) for v, _, _ in TERMS}   # plotted as density; box means stay native K

def sigpct(arrs, p=99.5):
    a = np.concatenate([np.abs(np.asarray(x)[np.isfinite(x)]) for x in arrs]); a = a[a > 1e-9]
    return float(np.nanpercentile(a, p)) if a.size else 1.0
# --- plot (density=True -> K per 10^6 km^2 [main];  density=False -> raw K/traj [comparison]) ---
plt.rcParams.update({"font.size": 14, "axes.titlesize": 15, "xtick.labelsize": 11, "ytick.labelsize": 11})
def make_fig(density, tag, unit):
    pdata = dmaps if density else maps
    vmax = sigpct([pdata[v] for v, _, _ in TERMS])
    proj = ccrs.PlateCarree(); fig = plt.figure(figsize=(15, 9.6), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.02, h_pad=0.02, wspace=0.01, hspace=0.05)
    axs = []
    for i, (var, letter, title) in enumerate(TERMS, 1):
        ax = fig.add_subplot(2, 2, i, projection=proj); ax.set_extent(EXT, crs=proj)
        ax.add_feature(cfeature.LAND, alpha=0.45); ax.add_feature(cfeature.OCEAN, alpha=0.2)
        ax.add_feature(cfeature.COASTLINE, lw=0.5); ax.add_feature(cfeature.BORDERS, lw=0.25)
        gl = ax.gridlines(draw_labels=True, ls="--", alpha=0.3); gl.top_labels = gl.right_labels = False
        gl.xlabel_style = gl.ylabel_style = {"size": 11}
        pm = ax.pcolormesh(lon_c, lat_c, np.where(np.abs(maps[var]) > 1e-12, pdata[var], np.nan),
                           cmap=TEMPDIV, vmin=-vmax, vmax=vmax, transform=proj)
        lo0, lo1, la0, la1 = MOR; ax.plot([lo0, lo1, lo1, lo0, lo0], [la0, la0, la1, la1, la0], "k--", lw=1.4, transform=proj)
        ax.text(-0.02, 1.04, f"({letter})", transform=ax.transAxes, fontsize=20, fontweight="bold",
                va="bottom", ha="left")                              # letter OUTSIDE panel
        ax.set_title(title, fontsize=15)
        ax.text(0.98, 0.03, f"box mean {boxmean[var]:+.3f} K/traj", transform=ax.transAxes, fontsize=11,
                va="bottom", ha="right", bbox=dict(boxstyle="round", fc="white", alpha=0.85))
        axs.append(ax)
    cb = fig.colorbar(pm, ax=axs, location="bottom", shrink=0.5, aspect=45, pad=0.02)
    cb.set_label(f"mean process contribution per in-band trajectory by genesis location  [{unit}]", fontsize=13)
    cb.ax.tick_params(labelsize=11)
    base = f"${HOME_WORK}/Supplement_ERA5_2026_origin_allprocess{tag}_windfix"
    fig.savefig(base + ".png", dpi=150, bbox_inches="tight")
    fig.savefig(base + ".pdf", bbox_inches="tight"); print("saved", base + ".png/.pdf")
    plt.close(fig)
make_fig(True,  "",          "K per 10$^6$ km$^2$")   # main (area density)
make_fig(False, "_rawKtraj", "K")                      # comparison (raw per-trajectory)
