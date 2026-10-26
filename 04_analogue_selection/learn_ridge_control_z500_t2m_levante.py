"""SUPPLEMENT: the ERA5 ridge control  Z500 -> T2m  learned over ALL AMJJA days.

Builds the climatological Z500->T2m transfer function the colleague suggested, so the
SAME control can be applied uniformly to the AIFS ensemble and the analogues (instead of
re-fitting a per-sample z500 slope). Everything is built exactly as Figures 1 & 5:

  T2m : target-box (-5.5..5 E, 44.5..55 N) FRACTIONAL cos-lat box mean of daily-mean t2m
        (comp_target_box_mean), anomaly vs the fixed 1995-2025 doy climatology (+-3 d),
        then RMST-detrended and pushed to present (event-year) climate.
  Z500: ridge box (-20..25 E, 30..68 N) cos-lat box mean (gpm), RMST-detrended to the
        same present climate (notebook's Z500_DETREND_SLOPE).
  RMST : European box (35..70 N, -15..30 E) cos-lat mean of daily t2m, 15-yr centred
        running mean, edges extrapolated from a linear fit to the first/last 30 yr.

Fit: OLS  T2m_anom_detr ~ Z500_detr  over AMJJA 1940-2025 (2026 excluded).

Outputs:
  ./supplement_ridge_control_t2m_z500_AMJJA.png / .pdf   (density scatter)
  ${DERIVED_DATA}/ridge_control_t2m_z500_AMJJA.json       (coefficients for Fig 2/3)
"""
import numpy as np, pandas as pd, xarray as xr, json, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.colors as mcolors
from scipy import stats

# ---- LEVANTE paths (same ERA5 daily files the Fig 5 analogue notebook uses). ----
# Override with env vars RIDGE_BASE / RIDGE_OUTDIR if your tree differs.
BASE   = os.environ.get("RIDGE_BASE", "${ANALOGUE_BASE}")
OUTDIR = os.environ.get("RIDGE_OUTDIR", ".")
T2M_FILE  = f"{BASE}/rawdat_0p5/era5_t2m_daily_NA_1940_2026_FULL_0p5deg.nc"
Z500_FILE = f"{BASE}/rawdat/era5_z500_daily_NA_1940_2026_FULL_1deg.nc"
G0 = 9.80665
TARGET = dict(lon=(-5.5, 5.0), lat=(44.5, 55.0))     # t2m target region
RMST_BOX = dict(lon=(-15, 30), lat=(35, 70))         # European RMST box
# Ridge Z500 = the LOCAL ridge over the target box (same region/def as Fig 2/4 z500_gh),
# NOT the large analogue-selection box, so the control matches the ridge-controlled analyses.
RIDGE_BOX = dict(lon=(-5.5, 5.0), lat=(44.5, 55.0))
CLIM_Y = (1995, 2025); CLIM_HALF = 3                 # fixed doy climatology
SEASON = os.environ.get("SEASON", "AMJJA")           # AMJJA / MJJA / JJA
MONTHS = {"AMJJA": (4, 8), "MJJA": (5, 8), "JJA": (6, 8)}[SEASON]
FIT_LAST = 2025                                       # slope fits exclude 2026
ROLL_YRS, EXTRAP_YRS = 15, 30
EVENT = pd.Timestamp("2026-05-25")
EVENT_REF = pd.Timestamp(year=min(EVENT.year, FIT_LAST), month=EVENT.month, day=EVENT.day)  # cap at 2025


def _open(path):
    da = xr.open_dataarray(path).squeeze(drop=True)
    ln = "lon" if "lon" in da.coords else "longitude"; lt = "lat" if "lat" in da.coords else "latitude"
    if float(da[ln].max()) > 180:
        da = da.assign_coords({ln: (((da[ln] + 180) % 360) - 180)}).sortby(ln)
    da = da.sortby(lt)
    return da, ln, lt


def coslat_box_mean(da, ln, lt, box):
    s = da.sel({ln: slice(*box["lon"])})
    s = s.sel({lt: slice(*box["lat"])}) if float(s[lt][0]) < float(s[lt][-1]) else s.sel({lt: slice(box["lat"][1], box["lat"][0])})
    w = np.cos(np.deg2rad(s[lt]))
    return s.weighted(w).mean(dim=(lt, ln))


def frac_box_mean(da, ln, lt, box):   # comp_target_box_mean (fractional cos-lat), Fig 1/5
    lon, lat = da[ln], da[lt]
    dlon = float(abs(lon.values[1] - lon.values[0])); dlat = float(abs(lat.values[1] - lat.values[0]))
    fx = ((np.minimum(lon + dlon/2, box["lon"][1]) - np.maximum(lon - dlon/2, box["lon"][0])) / dlon).clip(0, 1)
    fy = ((np.minimum(lat + dlat/2, box["lat"][1]) - np.maximum(lat - dlat/2, box["lat"][0])) / dlat).clip(0, 1)
    w = fx * fy * np.cos(np.deg2rad(lat)); v = np.isfinite(da)
    return (da.where(v) * w).sum((ln, lt)) / (w * v).sum((ln, lt))


def to_series(da):
    da = da.load()
    if float(da.mean()) > 150:   # K -> C for t2m
        da = da - 273.15
    s = pd.Series(np.asarray(da.values, float), index=pd.DatetimeIndex(da["time"].values).normalize())
    return s[~s.index.duplicated()]


# ---------------- RMST (identical to load_rmst_full / _analogue_rmst) ----------------
print("[1/5] RMST ...")
t2m_da, tln, tlt = _open(T2M_FILE)
rmst = to_series(coslat_box_mean(t2m_da, tln, tlt, RMST_BOX))
rmst = rmst.sort_index()
win = int(round(ROLL_YRS * 365.25))
sm = rmst.rolling(win, center=True, min_periods=win).mean()
valid = sm.dropna(); nfit = int(round(EXTRAP_YRS * 365.25))
def _edge(fit, anchor_t, anchor_v, times):
    x = (fit.index - fit.index[0]).days.values.astype(float); c = np.polyfit(x, fit.values, 1)
    xt = (times - fit.index[0]).days.values.astype(float); line = np.polyval(c, xt)
    xa = (anchor_t - fit.index[0]).days; return line - (np.polyval(c, xa) - anchor_v)
head = sm.index[sm.index < valid.index[0]]; tail = sm.index[sm.index > valid.index[-1]]
sm.loc[head] = _edge(valid.iloc[:nfit], valid.index[0], valid.iloc[0], head)
sm.loc[tail] = _edge(valid.iloc[-nfit:], valid.index[-1], valid.iloc[-1], tail)
RMST = sm
RMST_EVENT = float(RMST.reindex([EVENT_REF]).interpolate().iloc[0]) if EVENT_REF in RMST.index else float(RMST.asof(EVENT_REF))
print(f"      RMST_EVENT({EVENT_REF.date()}) = {RMST_EVENT:.3f} C")

# ---------------- T2m target-box anomaly (fixed 1995-2025 doy clim) + detrend ----------------
print("[2/5] T2m target anomaly ...")
tgt = to_series(frac_box_mean(t2m_da, tln, tlt, TARGET)).sort_index()
base = tgt[(tgt.index.year >= CLIM_Y[0]) & (tgt.index.year <= CLIM_Y[1])]
bdoy = base.index.dayofyear.values; bval = base.values
clim = np.full(367, np.nan)
for D in range(1, 367):
    d = np.abs(((bdoy - D + 182) % 365) - 182)
    m = d <= CLIM_HALF
    if m.any(): clim[D] = np.nanmean(bval[m])
anom = tgt.values - clim[tgt.index.dayofyear.values]
anom = pd.Series(anom, index=tgt.index)
rm_on_t = RMST.reindex(anom.index, method="nearest")
seas = (anom.index.month >= MONTHS[0]) & (anom.index.month <= MONTHS[1])
fit = seas & (anom.index.year <= FIT_LAST) & np.isfinite(anom.values) & np.isfinite(rm_on_t.values)
RMST_ANOM_SLOPE = float(np.polyfit(rm_on_t.values[fit], anom.values[fit], 1)[0])
anom_detr = anom + RMST_ANOM_SLOPE * (RMST_EVENT - rm_on_t.values)
print(f"      RMST_ANOMALY_SLOPE = {RMST_ANOM_SLOPE:+.3f} K/K ; event anom = {anom_detr.reindex([EVENT]).iloc[0]:+.2f} C")

# ---------------- Z500 ridge-box (gpm) + detrend to present climate ----------------
print("[3/5] Z500 ridge ...")
z_da, zln, zlt = _open(Z500_FILE)
zbox = frac_box_mean(z_da, zln, zlt, RIDGE_BOX).load()   # fractional cos-lat (= Fig 2/4 z500_gh)
zser = pd.Series(np.asarray(zbox.values, float), index=pd.DatetimeIndex(zbox["time"].values).normalize())
zser = zser[~zser.index.duplicated()].sort_index()
if zser.mean() > 2e4: zser = zser / G0     # geopotential -> gpm
rm_on_z = RMST.reindex(zser.index, method="nearest")
seaz = (zser.index.month >= MONTHS[0]) & (zser.index.month <= MONTHS[1])
fitz = seaz & (zser.index.year <= FIT_LAST) & np.isfinite(zser.values) & np.isfinite(rm_on_z.values)
Z500_SLOPE = float(np.polyfit(rm_on_z.values[fitz], zser.values[fitz], 1)[0])   # gpm per K
z_detr = zser + Z500_SLOPE * (RMST_EVENT - rm_on_z.values)                      # push to present climate
print(f"      Z500_DETREND_SLOPE = {Z500_SLOPE:+.2f} gpm/K")

# ---------------- learned ridge control: T2m_anom_detr ~ Z500_detr (AMJJA 1940-2025) ----------------
print("[4/5] regression ...")
df = pd.DataFrame({"t": anom_detr, "z": z_detr}).dropna()
df = df[(df.index.month >= MONTHS[0]) & (df.index.month <= MONTHS[1]) & (df.index.year <= FIT_LAST)]
BETA, INTERCEPT = np.polyfit(df.z.values, df.t.values, 1)
r, p = stats.pearsonr(df.z.values, df.t.values); R2 = r**2
n = len(df)
# event point (2026), detrended the same way
z_event = float(z_detr.reindex([EVENT]).iloc[0]) if EVENT in z_detr.index else np.nan
t_event = float(anom_detr.reindex([EVENT]).iloc[0]) if EVENT in anom_detr.index else np.nan
print(f"      BETA = {BETA:+.4f} K/gpm  intercept = {INTERCEPT:+.2f}  r = {r:+.3f} (R2={R2:.2f})  n = {n}")

out = dict(beta_K_per_gpm=BETA, intercept_K=INTERCEPT, r=r, r2=R2, n=int(n), p_value=float(p),
           rmst_anomaly_slope_K_per_K=RMST_ANOM_SLOPE, z500_detrend_slope_gpm_per_K=Z500_SLOPE,
           rmst_event_C=RMST_EVENT, event_ref=str(EVENT_REF.date()),
           target_box=TARGET, ridge_box=RIDGE_BOX, rmst_box=RMST_BOX,
           clim_years=list(CLIM_Y), clim_half_window_days=CLIM_HALF, season=SEASON, months=list(MONTHS),
           fit_last_year=FIT_LAST, reference="present (event-year) climate on both axes",
           event_point=dict(z500_detr_gpm=z_event, t2m_anom_detr_C=t_event),
           note=("Ridge-controlled T2m anomaly of a case = t2m_anom_detr - (beta*z500_detr + intercept). "
                 "Detrend a case's z500 to present climate with z500_detr = z500_gpm + z500_detrend_slope*(rmst_event - rmst_year)."))
os.makedirs(OUTDIR, exist_ok=True)
JOUT = f"{OUTDIR}/ridge_control_t2m_z500_{SEASON}.json"
json.dump(out, open(JOUT, "w"), indent=2); print("saved", JOUT)

# ---------------- density scatter ----------------
print("[5/5] plot ...")
plt.rcParams.update({"font.size": 13})
fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
hb = ax.hexbin(df.z.values, df.t.values, gridsize=70, mincnt=1, bins="log", cmap="magma_r")
xx = np.linspace(df.z.min(), df.z.max(), 100)
ax.plot(xx, BETA*xx + INTERCEPT, "k-", lw=2.2, label=f"OLS: {BETA:+.3f} K per gpm")
if np.isfinite(z_event):
    ax.scatter([z_event], [t_event], s=90, c="#1f77b4", edgecolor="k", zorder=6, label=f"25 May 2026 ({t_event:+.1f} °C)")
ax.axhline(0, color="0.6", lw=0.8, ls=":"); ax.axvline(float(z_detr.reindex(df.index).mean()), color="0.6", lw=0.8, ls=":")
ax.set_xlabel("Detrended ridge $Z_{500}$  [gpm, present climate]")
ax.set_ylabel("Detrended target-box T2m anomaly  [°C]")
ax.set_title(f"ERA5 ridge control ({SEASON} 1940–2025): $r$={r:+.2f}, $R^2$={R2:.2f}, $n$={n}", fontsize=13)
ax.legend(loc="upper left", frameon=True, fontsize=11)
cb = fig.colorbar(hb, ax=ax, label="days per bin (log)")
for ext in ("png", "pdf"):
    fig.savefig(f"{OUTDIR}/supplement_ridge_control_t2m_z500_{SEASON}.{ext}", dpi=200, bbox_inches="tight")
print(f"saved supplement_ridge_control_t2m_z500_{SEASON}.png/.pdf")
