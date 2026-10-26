"""Figure 3 (regression) — OLS slope of each process on cut-off activity.

Two rows x four columns (pressure x arrival-day Hovmoller):
  row 1 = IFS members (n~100)      : slope of process on COL_d4/NF
  row 2 = best-25 cut-off analogues : slope of process on COL_ACT/NF
Shading = OLS slope  beta = d(process)/d(COL)  [degC per mean daily COL activity];
black contour = |r| = HL (0.5); hatching = NOT significant (p > PSIG), Fig-2d style.
Whole arrival box (no land restriction). X-axis = arrival day rel. to the peak.

ANALOGUES now come from the NEW corrected-6h pipeline, one file per rank:
  analogue_origin_6h/rank<NN>_MJ_<date>/origin_contributions.nc
(each carries lon0/lat0, arrival_p, TA/advection/adiabatic/diabatic + geneses).
"best 25" = ranks 01-25 (the folders are ordered by analogue quality, all cut-off).

Output: ./Fig3_slope_COL_IFS_analogues_N25.png  (+ npz cache next to data)
"""
import numpy as np, pandas as pd, os, glob, re
from pathlib import Path
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

# ============================ CONFIG ============================
D = "${DERIVED_DATA}"
# IFS: the colleague's tanom_decomp of 2026-10-19 (corrected advective T')
IFS = Path(os.environ.get("IFS_BASE", "${IFS_TANOM_DECOMP}"))
# analogues: windfix + genesis as in Hotz et al. (2024), rsynced from Levante analogues_windfix
ANA_DIR = os.environ.get("ANA_DIR", "${DERIVED_DATA}/analogue_origin_6h_windfix")
# output + cache tag from the data sources, so a cache of other inputs is never reused
_m = re.search(r"tanom_decomp_(\d{8})", str(IFS))                # data set of the colleague (delivery date)
TAG = "_ifs" + (_m.group(1) if _m else "old") + "_ana" + ("windfix" if "windfix" in ANA_DIR else "OLD")
IFS_CSV = D + "/ifs_COL_window_sweep_DAILYMEAN-5.csv"
CATALOGUE = D + "/col_features_AMJJA_1940_2025_largebox.nc"
N_ANA = int(os.environ.get("NANA", "25"))                         # best-N cut-off analogues (ranks 01..N)
# EXCLUDE=YYYYMMDD[,YYYYMMDD]: leave these analogues out; the next ranks fill up, so n stays N_ANA
EXCLUDE = [x for x in os.environ.get("EXCLUDE", "").replace(" ", "").split(",") if x]
if EXCLUDE: TAG += "_without" + "_".join(EXCLUDE)
OUT = f"${HOME_WORK}/Fig3_slope_COL_IFS_analogues_N{N_ANA}{TAG}.png"
CACHE = D + f"/fig3_slope_COL_IFS_analogues_N{N_ANA}{TAG}_cache.npz"    # delete to recompute

NF = 15.0                 # COL -> mean daily activity (slope units: degC per mean daily COL)
HL = 0.5                  # highlight |r| contour
PSIG = 0.05               # significance level (hatch where p > PSIG)
MIN_N = 20                # min group members with data in a cell
# 0 = all arrivals (paper figure); 1 = land arrivals only (lon0/lat0), both rows - not used in the paper
LAND_ONLY = os.environ.get("LAND_ONLY", "0") == "1"
LAND_ONLY_IFS = LAND_ONLY
if LAND_ONLY:
    OUT = OUT.replace(".png", "_LAND.png"); CACHE = CACHE.replace(".npz", "_LAND.npz")
HOURS = ["00", "06", "12", "18"]
IFS_CENTER = pd.Timestamp("2026-05-25")

P_EDGES = np.arange(175., 1026., 50.); P_CEN = 0.5 * (P_EDGES[:-1] + P_EDGES[1:])
REL_EDGES = np.arange(-7.125, 5.126, 0.25); REL_CEN = 0.5 * (REL_EDGES[:-1] + REL_EDGES[1:])
ORDER = ["TA", "advective", "adiabatic", "diabatic"]
IFSCOLS = {"TA": ["tanom0"], "advective": ["tanom_advective"],
           "adiabatic": ["tanom_adiab1", "tanom_adiab2"], "diabatic": ["tanom_diabatic"]}
ANA_VARMAP = {"TA": "TA", "advective": "advection", "adiabatic": "adiabatic", "diabatic": "diabatic"}


def hov(V, P, R):
    df = pd.DataFrame({"v": V, "p": P, "r": R}).dropna()
    df = df[(df.p >= P_EDGES[0]) & (df.p < P_EDGES[-1]) & (df.r >= REL_EDGES[0]) & (df.r < REL_EDGES[-1])]
    df["pb"] = np.digitize(df.p, P_EDGES) - 1
    df["ri"] = np.digitize(df.r, REL_EDGES) - 1
    g = df.groupby(["pb", "ri"])["v"].mean()
    G = np.full((len(P_CEN), len(REL_CEN)), np.nan)
    for (pb, ri), val in g.items():
        if 0 <= int(pb) < len(P_CEN) and 0 <= int(ri) < len(REL_CEN): G[int(pb), int(ri)] = val
    return G


def col_act_series():
    import xarray as xr
    ct = xr.open_dataset(CATALOGUE); cf = ct[["time", "So", "ex"]].to_dataframe().reset_index(); ct.close()
    cc = cf[(cf.ex == 1) & (cf.So >= 6)].copy(); cc["day"] = pd.DatetimeIndex(cc.time).normalize()
    dd = cc.groupby("day")["So"].sum(); dd = dd.reindex(pd.date_range(dd.index.min(), dd.index.max(), freq="D"), fill_value=0.0)
    return dd.rolling(5, min_periods=1).sum()


def compute():
    import xarray as xr
    is_land = None
    if LAND_ONLY_IFS:
        from cartopy.io import shapereader
        from shapely.ops import unary_union
        from shapely import contains_xy
        _LAND = unary_union(list(shapereader.Reader(
            shapereader.natural_earth(resolution="50m", category="physical", name="land")).geometries()))
        def is_land(lon, lat):
            return contains_xy(_LAND, ((np.asarray(lon, float) + 180) % 360) - 180, np.asarray(lat, float))

    # ---- IFS ----
    col = pd.read_csv(IFS_CSV); col["initdir"] = col["init"].str.replace("-", "") + "_00"
    col["mem"] = col["member"].map(lambda m: f"{int(m):02d}")
    COL = {(r.initdir, r.mem): r.COL_d4 / NF for r in col.itertuples()}
    ifs_grids = {t: [] for t in ORDER}; ifs_col = []
    read_dates = pd.date_range("2026-05-18", "2026-05-30", freq="D")
    for initd in sorted({k[0] for k in COL}):
        for mem in [f"{i:02d}" for i in range(1, 51)]:
            if (initd, mem) not in COL: continue
            rows = []
            for Dd in read_dates:
                for h in HOURS:
                    p = IFS / initd / mem / f"tanom_{Dd:%Y%m%d}_{h}"
                    if p.is_file():
                        df = pd.read_csv(p, sep=r"\s+")
                        df["_rel"] = ((Dd + pd.Timedelta(hours=int(h))) - IFS_CENTER) / pd.Timedelta(days=1)
                        rows.append(df)
            if not rows: continue
            d = pd.concat(rows, ignore_index=True)
            for c in ["p0", "lon0", "lat0"] + sum(IFSCOLS.values(), []): d[c] = pd.to_numeric(d[c], errors="coerce")
            if LAND_ONLY_IFS:
                d = d[np.isfinite(d.lon0) & np.isfinite(d.lat0)]
                d = d[is_land(d.lon0.values, d.lat0.values)]
            for t in ORDER:
                ifs_grids[t].append(hov(d[IFSCOLS[t]].sum(axis=1).values, d.p0.values, d._rel.values))
            ifs_col.append(COL[(initd, mem)])
    ifs_col = np.array(ifs_col, float)

    # ---- analogues (NEW pipeline, best-N = ranks 01..N) ----
    COL_ACT = col_act_series()
    folders = []
    for f in glob.glob(ANA_DIR + "/rank*_MJ_*"):
        m = re.search(r"rank(\d+)_MJ_(\d{8})", os.path.basename(f))
        if m and os.path.exists(f + "/origin_contributions.nc"):
            folders.append((int(m.group(1)), m.group(2), f))
    folders = [x for x in sorted(folders) if x[1] not in EXCLUDE][:N_ANA]
    ana_grids = {t: [] for t in ORDER}; ana_col = []
    for rk, ymd, folder in folders:
        date = pd.Timestamp(ymd)
        cval = float(COL_ACT.reindex([date.normalize()]).iloc[0]) if date.normalize() in COL_ACT.index else np.nan
        if not np.isfinite(cval): continue
        ds = xr.open_dataset(folder + "/origin_contributions.nc")
        tt = pd.DatetimeIndex(ds["time"].values); ap = ds["arrival_p"].values
        rel = ((tt - date) / pd.Timedelta(days=1)).values
        relg = np.repeat(rel[:, None], ap.shape[1], axis=1)
        keep = np.ones(ap.size, bool)
        if LAND_ONLY:
            lo, la = ds["lon0"].values.ravel(), ds["lat0"].values.ravel()
            keep = np.isfinite(lo) & np.isfinite(la)
            keep[keep] = is_land(lo[keep], la[keep])
        tmp = {t: hov(ds[ANA_VARMAP[t]].values.ravel()[keep], ap.ravel()[keep], relg.ravel()[keep]) for t in ORDER}
        ds.close()
        for t in ORDER: ana_grids[t].append(tmp[t])
        ana_col.append(cval / NF)
    ana_col = np.array(ana_col, float)

    def slope_r(stack, cvec):
        S = np.array(stack); SL = np.full((len(P_CEN), len(REL_CEN)), np.nan); R = SL.copy(); P = SL.copy()
        for i in range(len(P_CEN)):
            for j in range(len(REL_CEN)):
                y = S[:, i, j]; m = np.isfinite(y)
                if m.sum() >= MIN_N and np.nanstd(y[m]) > 0 and np.nanstd(cvec[m]) > 0:
                    SL[i, j] = np.polyfit(cvec[m], y[m], 1)[0]
                    rr, pp = stats.pearsonr(cvec[m], y[m]); R[i, j] = rr; P[i, j] = pp
        return SL, R, P

    out = {"n_ifs": len(ifs_col), "n_ana": len(ana_col)}
    for t in ORDER:
        sl, r, p = slope_r(ifs_grids[t], ifs_col); out[f"IFS_SL_{t}"] = sl; out[f"IFS_R_{t}"] = r; out[f"IFS_P_{t}"] = p
        sl, r, p = slope_r(ana_grids[t], ana_col); out[f"ANA_SL_{t}"] = sl; out[f"ANA_R_{t}"] = r; out[f"ANA_P_{t}"] = p
    np.savez(CACHE, **out)
    return out


data = np.load(CACHE) if os.path.exists(CACHE) else compute()
n_ifs = int(data["n_ifs"]); n_ana = int(data["n_ana"])
print(f"IFS members={n_ifs}, analogues={n_ana}")

# symmetric colour scale (0.1 steps) from 98th pct of |slope| over |r|>0.2 cells
allsl = []
for row in ("IFS", "ANA"):
    for t in ORDER:
        sl = data[f"{row}_SL_{t}"]; r = data[f"{row}_R_{t}"]
        allsl.append(sl[np.isfinite(sl) & (np.abs(r) > 0.2)])
allsl = np.concatenate([a for a in allsl if a.size]) if allsl else np.array([0.5])
vmax = max(0.1, np.ceil(np.nanpercentile(np.abs(allsl), 98) / 0.1) * 0.1)
lev = np.arange(-vmax, vmax + 0.1, 0.1)

fs = 16
plt.rcParams.update({"font.size": fs, "xtick.labelsize": fs - 3, "ytick.labelsize": fs - 3,
                     "hatch.linewidth": 0.7})
XLIM = (-5, 4)
fig, axes = plt.subplots(2, 4, figsize=(15, 8.2), sharex=True, sharey=True, constrained_layout=True)
rowinfo = [("IFS", f"{'IFS, land' if LAND_ONLY_IFS else 'IFS'} (n={n_ifs})"),
           ("ANA", (f"Best-{n_ana} analogues\nwithout " + ", ".join(x[:4] for x in EXCLUDE))
                   if EXCLUDE else f"Best-{n_ana} cut-off analogues{', land' if LAND_ONLY else ''}")]
cf = None
for ri, (key, rlab) in enumerate(rowinfo):
    for ci, term in enumerate(ORDER):
        ax = axes[ri, ci]
        SL = data[f"{key}_SL_{term}"]; R = data[f"{key}_R_{term}"]; P = data[f"{key}_P_{term}"]
        cf = ax.contourf(REL_CEN, P_CEN, SL, levels=lev, cmap="RdBu_r", extend="both")
        ax.contour(REL_CEN, P_CEN, np.abs(R), levels=[HL], colors="k", linewidths=1.4)
        Pm = np.ma.masked_invalid(P)
        ax.contourf(REL_CEN, P_CEN, Pm, levels=[PSIG, 1.0], hatches=["...."], colors="none")
        ax.set_yscale("log"); ax.set_ylim(1000, 200)
        ax.set_yticks([1000, 850, 700, 500, 400, 300, 200]); ax.set_yticklabels([1000, 850, 700, 500, 400, 300, 200]); ax.minorticks_off()
        ax.axvline(0, color="k", ls="--", lw=1); ax.set_xlim(*XLIM)
        if ri == 0: ax.set_title("TA" if term == "TA" else term.capitalize(), fontsize=fs + 1)
        if ci == 0: ax.set_ylabel(f"{rlab}\nPressure (hPa)", fontsize=fs)
        if ri == 1: ax.set_xlabel("Arrival day (rel. to 25th May)")
        ax.text(-0.02, 1.04, f"({chr(97 + ri * 4 + ci)})", transform=ax.transAxes, fontsize=fs + 1, fontweight="bold", va="bottom", ha="left")
cb = fig.colorbar(cf, ax=axes, location="bottom", shrink=0.5, aspect=45, pad=0.06)
cb.set_label(r"$\beta$  [°C per mean daily COL activity]   (contour $|r|=%.1f$; hatch: $p>%.2f$)" % (HL, PSIG), fontsize=fs - 1)
fig.savefig(OUT, dpi=150, bbox_inches="tight"); fig.savefig(OUT.replace(".png", ".pdf"), bbox_inches="tight"); print("saved", OUT, "+ pdf")
