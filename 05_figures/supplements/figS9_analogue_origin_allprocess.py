"""SUPPLEMENT (new pipeline): BEST-25 cut-off analogue origin maps, ALL processes + TA.
Panels:  (a) Total anomaly (TA)   (b) Advection   (c) Adiabatic   (d) Diabatic

Analogues = ranks 01..25 of the new-pipeline cut-off set (all cut-off, ordered by
quality). Per-analogue arrivals in [D-1, D], 250-750 hPa origin band, 0.5 deg.
Each process is binned at ITS OWN genesis location (<proc>_lon/lat_genesis).
Trajectory data: analogue_origin_6h/rank<NN>_MJ_<date>/origin_contributions.nc.

Normalization (Fig-4 origin-map convention): for EACH analogue a SHARED denominator
nb (trajectories finite in ALL four terms AND in-band) normalizes every term; the
analogues are then averaged with equal weight. Saves PNG + PDF.
"""
import numpy as np, pandas as pd, xarray as xr, glob, os, re
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

BAND = (250, 750); MOR = (-20., 0., 25., 40.); DLON = 0.5; DAYS_BEFORE = 1; EXT = [-130, 30, 5, 75]
ANA_DIR = os.environ.get("ANA_DIR", "${DERIVED_DATA}/analogue_origin_6h_windfix"); N_BEST = 25   # windfix + new genesis
lon_bins = np.arange(-180, 180+DLON, DLON); lat_bins = np.arange(0, 90+DLON, DLON)
lon_c = 0.5*(lon_bins[:-1]+lon_bins[1:]); lat_c = 0.5*(lat_bins[:-1]+lat_bins[1:])
def wrap(x): return ((np.asarray(x)+180) % 360)-180

TERMS = [("TA",        "a", "Total anomaly (TA)"),
         ("advection", "b", "Advection"),
         ("adiabatic", "c", "Adiabatic"),
         ("diabatic",  "d", "Diabatic")]


def analogue_maps(folder, day):
    """Per-analogue term maps, each at its OWN genesis, shared (all-term & in-band) denominator."""
    ds = xr.open_dataset(folder + "/origin_contributions.nc")
    ti = pd.DatetimeIndex(ds["time"].values).normalize()
    sel = np.where((ti >= day-pd.Timedelta(days=DAYS_BEFORE)) & (ti <= day))[0]
    dd = ds.isel(time=sel)
    ap = dd["arrival_p"].values.astype(float).ravel()
    val = {v: dd[v].values.ravel() for v, _, _ in TERMS}
    gen = {v: (wrap(dd[f"{v}_lon_genesis"].values.ravel()), dd[f"{v}_lat_genesis"].values.ravel()) for v, _, _ in TERMS}
    ds.close()
    finite_all = np.ones_like(ap, bool)
    for v, _, _ in TERMS: finite_all &= np.isfinite(val[v])
    ok = finite_all & (ap >= BAND[0]) & (ap <= BAND[1]); n = int(ok.sum())
    if n == 0: return None
    out = {}
    for v, _, _ in TERMS:
        lo, la = gen[v]; m = ok & np.isfinite(lo) & np.isfinite(la)
        g, _, _ = np.histogram2d(lo[m], la[m], bins=[lon_bins, lat_bins], weights=val[v][m])
        out[v] = (g.T)/n                                         # Fig-4 normalization: K per in-band trajectory
    return out


folders = []
for f in glob.glob(ANA_DIR + "/rank*_MJ_*"):
    mm = re.search(r"rank(\d+)_MJ_(\d{8})", os.path.basename(f))
    if mm and os.path.exists(f + "/origin_contributions.nc"):
        folders.append((int(mm.group(1)), mm.group(2), f))
folders = sorted(folders)[:N_BEST]
stacks = {v: [] for v, _, _ in TERMS}; used = []
for rk, ymd, folder in folders:
    m = analogue_maps(folder, pd.Timestamp(ymd).normalize())
    if m is None: print("empty", ymd); continue
    for v, _, _ in TERMS: stacks[v].append(m[v])
    used.append(ymd)
maps = {v: np.nanmean(np.array(stacks[v]), axis=0) for v, _, _ in TERMS}
print(f"new-pipeline best-{N_BEST} cut-off analogues: {len(used)} used")

Xg, Yg = np.meshgrid(lon_c, lat_c)
inbox = (Xg >= MOR[0]) & (Xg <= MOR[1]) & (Yg >= MOR[2]) & (Yg <= MOR[3])
boxmean = {v: float(np.nansum(np.where(inbox & np.isfinite(maps[v]), maps[v], 0.0))) for v, _, _ in TERMS}
print(f"box means [K/traj]: " + "  ".join(f"{v}={boxmean[v]:+.3f}" for v, _, _ in TERMS))


def area_density(F, LTC):
    R = 6371.0; coslat = np.cos(np.deg2rad(np.asarray(LTC, float)))[:, None]
    return F/((np.deg2rad(DLON)*R)*(np.deg2rad(DLON)*R*coslat))*1e6
dmaps = {v: area_density(maps[v], lat_c) for v, _, _ in TERMS}

def sigpct(arrs, p=99.5):
    a = np.concatenate([np.abs(np.asarray(x)[np.isfinite(x)]) for x in arrs]); a = a[a > 1e-9]
    return float(np.nanpercentile(a, p)) if a.size else 1.0

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
        ax.text(-0.02, 1.04, f"({letter})", transform=ax.transAxes, fontsize=20, fontweight="bold", va="bottom", ha="left")
        ax.set_title(title, fontsize=15)
        ax.text(0.98, 0.03, f"box mean {boxmean[var]:+.3f} K/traj", transform=ax.transAxes, fontsize=11,
                va="bottom", ha="right", bbox=dict(boxstyle="round", fc="white", alpha=0.85))
        axs.append(ax)
    cb = fig.colorbar(pm, ax=axs, location="bottom", shrink=0.5, aspect=45, pad=0.02)
    cb.set_label(f"mean process contribution per in-band trajectory by genesis location  [{unit}]", fontsize=13)
    cb.ax.tick_params(labelsize=11)
    base = f"${HOME_WORK}/Supplement_analogue25_origin_allprocess{tag}_windfix"
    fig.savefig(base + ".png", dpi=150, bbox_inches="tight")
    fig.savefig(base + ".pdf", bbox_inches="tight"); print("saved", base + ".png/.pdf")
    plt.close(fig)
make_fig(True,  "",          "K per 10$^6$ km$^2$")
make_fig(False, "_rawKtraj", "K")
