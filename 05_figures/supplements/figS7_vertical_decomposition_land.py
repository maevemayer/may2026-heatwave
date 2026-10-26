"""Figure 3 (vertical process decomposition) — LAND-ARRIVAL only.

Land-on-arrival version of figure3_vertical_decomposition.py: only trajectories
whose ARRIVAL grid cell is over land are kept.

ERA5 (reality) row is read from the pooled origin-contribution file
(pooled_origin_contributions.nc), which stores per trajectory the EXACT arrival
point (lon0, lat0, p0/arrival_p), TA at arrival, and the cumulative process
contributions (advection/adiabatic/diabatic) — so the arrival land mask is a
straight is_land(lon0, lat0) filter. IFS row is the tanom_decomp composite,
land-masked on (lon0, lat0).

Set LAND_ONLY = False to reproduce the all-arrival baseline through this code.
Needs a kernel with cartopy (Natural Earth land polygons) + h5netcdf.
"""
import numpy as np, pandas as pd, xarray as xr
from pathlib import Path
import os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.dates as mdates, matplotlib.colors as mcolors

# ============================ CONFIG ============================
LAND_ONLY = True            # True = land arrivals only (this figure); False = all arrivals
LAND_RES  = "50m"           # Natural Earth land resolution ("50m" fast, "10m" precise)

POOLED_FILE = "${DERIVED_DATA}/pooled_origin_contributions.nc"   # ERA5 row
IFS_BASE = Path("${IFS_TANOM_DECOMP}")                   # IFS row
IFS_CSV  = "${DERIVED_DATA}/ifs_cutoff_activity_with_T2m_Z500.csv"
OUT = "${HOME_WORK}/Fig3_vertical_2row_ERA5_LAND.png"

fs = 17
axfs = fs - 2
plt.rcParams.update({"font.size": fs, "xtick.labelsize": axfs, "ytick.labelsize": axfs})

HOURS = ["00", "06", "12", "18"]
TIMES = pd.date_range("2026-05-18 00:00", "2026-05-30 18:00", freq="6h"); TIDX = {t: i for i, t in enumerate(TIMES)}
P_EDGES = np.arange(175., 1026., 50.); P_CEN = 0.5 * (P_EDGES[:-1] + P_EDGES[1:])
TERMS = ["TA", "advective", "adiabatic", "diabatic"]
E5VAR = {"TA": "TA", "advective": "advection", "adiabatic": "adiabatic", "diabatic": "diabatic"}   # pooled-file names
IFSCOLS = {"TA": ["tanom0"], "advective": ["tanom_advective"],
           "adiabatic": ["tanom_adiab1", "tanom_adiab2"], "diabatic": ["tanom_diabatic"]}

# ---------------------- land mask (arrival point) ----------------------
from cartopy.io import shapereader
from shapely.ops import unary_union
from shapely.prepared import prep
try:
    from shapely import contains_xy
except ImportError:
    contains_xy = None
    from shapely.geometry import Point
_LAND = unary_union(list(shapereader.Reader(
    shapereader.natural_earth(resolution=LAND_RES, category="physical", name="land")).geometries()))
_PREP = prep(_LAND)
def is_land(lon, lat):
    lon = ((np.asarray(lon, float) + 180) % 360) - 180; lat = np.asarray(lat, float)
    if contains_xy is not None:
        return contains_xy(_LAND, lon, lat)
    return np.array([_PREP.contains(Point(x, y)) for x, y in zip(lon, lat)])

# ---------------------- Hovmoller binning (verbatim) ----------------------
def hov(V, P, T):
    df = pd.DataFrame({"v": V, "p": P, "t": T}).dropna(); df = df[(df.p >= P_EDGES[0]) & (df.p < P_EDGES[-1])]
    df["ti"] = df.t.map(TIDX); df = df.dropna(subset=["ti"]); df["pb"] = np.digitize(df.p, P_EDGES) - 1
    g = df.groupby(["pb", "ti"])["v"].mean(); G = np.full((len(P_CEN), len(TIMES)), np.nan)
    for (pb, ti), val in g.items():
        if 0 <= int(pb) < len(P_CEN): G[int(pb), int(ti)] = val
    return G

# ---------------------- ERA5 row: from pooled origin file, land-masked arrivals ----------------------
def era5_hov():
    ds = xr.open_dataset(POOLED_FILE)
    ntr = ds.sizes["traj_id"]
    lon0 = (((ds["lon0"].values + 180) % 360) - 180).ravel()
    lat0 = ds["lat0"].values.ravel()
    ap = ds["arrival_p"].values.ravel()
    tt = pd.DatetimeIndex(ds["time"].values)
    tgrid = np.repeat(tt.values[:, None], ntr, axis=1).ravel()
    keep = np.isfinite(lon0) & np.isfinite(lat0)
    if LAND_ONLY:
        keep = keep & is_land(np.where(np.isfinite(lon0), lon0, 0.0), np.where(np.isfinite(lat0), lat0, 0.0))
    out = {}
    for term in TERMS:
        val = ds[E5VAR[term]].values.ravel()
        m = keep & np.isfinite(val)
        out[term] = hov(val[m], ap[m], tgrid[m])
    ds.close()
    return out

# ---------------------- IFS row: land-masked arrivals (tanom_decomp) ----------------------
csv = pd.read_csv(IFS_CSV); csv["initdir"] = csv["init"].str.replace("-", "") + "_00"; csv["mem"] = csv["member"].map(lambda m: f"{int(m):02d}")
topT = list(zip(csv.nlargest(5, "t2m_C").initdir, csv.nlargest(5, "t2m_C").mem))
topZ = list(zip(csv.nlargest(5, "z500_gh").initdir, csv.nlargest(5, "z500_gh").mem))

def ifs_member(initd, mem):
    rows = []
    for D in pd.date_range("2026-05-18", "2026-05-30", freq="D"):
        for h in HOURS:
            p = IFS_BASE / initd / mem / f"tanom_{D:%Y%m%d}_{h}"
            if p.is_file():
                df = pd.read_csv(p, sep=r"\s+"); df["_t"] = D + pd.Timedelta(hours=int(h)); rows.append(df)
    if not rows: return None
    d = pd.concat(rows, ignore_index=True)
    for c in ["p0", "lon0", "lat0"] + sum(IFSCOLS.values(), []):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    if LAND_ONLY:
        d = d[np.isfinite(d.lon0) & np.isfinite(d.lat0)]
        d = d[is_land(d.lon0.values, d.lat0.values)]
    return {term: hov(d[IFSCOLS[term]].sum(axis=1).values, d.p0.values, d._t.values) for term in TERMS}

def comp(members):
    acc = {t: [] for t in TERMS}
    for i, m in members:
        r = ifs_member(i, m)
        if r:
            for t in TERMS: acc[t].append(r[t])
    return {t: np.nanmean(np.array(acc[t]), axis=0) for t in TERMS}

# ---------------------- build + plot ----------------------
print(f"LAND_ONLY={LAND_ONLY}  (Natural Earth {LAND_RES})")
era5 = era5_hov(); compT = comp(topT); compZ = comp(topZ)
rows = [("ERA5 (reality)", era5), ("Top-5 (T2m − Z500)", {t: compT[t] - compZ[t] for t in TERMS})]

tnum = mdates.date2num(TIMES.to_pydatetime()); fig, axes = plt.subplots(2, 4, figsize=(21, 7.5), sharex=True, sharey=True); cf_d = None
for ri, (rlab, rd) in enumerate(rows):
    if ri == 0:
        vmax = np.nanpercentile(np.abs(np.concatenate([rd[t].ravel() for t in TERMS])), 99); vmi = max(int(np.ceil(vmax)), 3)
        if vmi % 2 == 0: vmi += 1
        levels = np.arange(-vmi, vmi + 2, 2).astype(float)
    else:
        levels = np.arange(-7., 8., 2.)
    for ci, term in enumerate(TERMS):
        ax = axes[ri, ci]; cf = ax.contourf(tnum, P_CEN, rd[term], levels=levels, cmap="RdBu_r", extend="both")
        yt = np.arange(1000, 199, -200)
        ax.set_yscale("log"); ax.set_ylim(1000, 200); ax.set_yticks(yt); ax.minorticks_off(); ax.set_yticklabels(yt)
        ax.axvline(mdates.date2num(pd.Timestamp("2026-05-25")), color="k", ls="--", lw=1)
        if ri == 0: ax.set_title("TA" if term == "TA" else term.capitalize(), fontsize=fs)
        if ci == 0: ax.set_ylabel(f"{rlab}\nPressure (hPa)", fontsize=fs)
        ax.text(-0.02, 1.04, f"({chr(97 + ri * 4 + ci)})", transform=ax.transAxes, fontsize=fs + 2, fontweight="bold", va="bottom", ha="left")
    if ri == 0: fig.colorbar(cf, ax=axes[0, :].tolist(), shrink=0.85 * 1.15, pad=0.01, ticks=levels, label="[K]")
    else: cf_d = cf
fig.colorbar(cf_d, ax=axes[1, :].tolist(), shrink=0.85, pad=0.01, ticks=np.arange(-7, 8, 2), label="Δ [K]")
for ax in axes[-1, :]:
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2)); ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.tick_params(axis="x", labelrotation=90); ax.set_xlabel("Arrival date")
fig.suptitle("Vertical process decomposition — LAND arrival trajectories only" if LAND_ONLY else
             "Vertical process decomposition — all arrivals", fontsize=fs + 1, y=1.02)
fig.savefig(OUT, dpi=150, bbox_inches="tight"); print("saved", OUT)
