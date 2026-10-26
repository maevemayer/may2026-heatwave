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

OC = "${ERA5_DERIVED}"
T2M_FILE  = f"{OC}/era5_t2m_daily_NA_1940_2026_FULL_0p5deg.nc"
Z500_FILE = f"{OC}/era5_z500_daily_NA_1940_2026_FULL_1deg.nc"
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
BOOT_B, BOOT_SEED = 1000, 0                           # bootstrap 95% CI (same scheme as Fig 5 temperature slope)


def _boot_slope_ci(x, y, B=BOOT_B, seed=BOOT_SEED):
    """95% CI of an OLS slope by B-sample nonparametric bootstrap (resample fit days w/ replacement).
    Returns (median, lo2.5, hi97.5) — identical scheme to the Fig 5 Robin/temperature-slope CIs."""
    x = np.asarray(x, float); y = np.asarray(y, float); n = x.size
    rng = np.random.default_rng(seed)
    idx = rng.choice(np.arange(n), size=(B, n), replace=True)
    sl = np.array([np.polyfit(x[i], y[i], 1)[0] for i in idx])
    return float(np.median(sl)), float(np.quantile(sl, 0.025)), float(np.quantile(sl, 0.975))


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
RMST_ANOM_SLOPE_CI = _boot_slope_ci(rm_on_t.values[fit], anom.values[fit])     # (median, lo, hi)
anom_detr = anom + RMST_ANOM_SLOPE * (RMST_EVENT - rm_on_t.values)
print(f"      RMST_ANOMALY_SLOPE = {RMST_ANOM_SLOPE:+.3f} K/K  95% CI [{RMST_ANOM_SLOPE_CI[1]:+.3f}, {RMST_ANOM_SLOPE_CI[2]:+.3f}] ; event anom = {anom_detr.reindex([EVENT]).iloc[0]:+.2f} C")

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
Z500_SLOPE_CI = _boot_slope_ci(rm_on_z.values[fitz], zser.values[fitz])         # (median, lo, hi)
z_detr = zser + Z500_SLOPE * (RMST_EVENT - rm_on_z.values)                      # push to present climate
print(f"      Z500_DETREND_SLOPE = {Z500_SLOPE:+.2f} gpm/K  95% CI [{Z500_SLOPE_CI[1]:+.2f}, {Z500_SLOPE_CI[2]:+.2f}]")

# ---------------- Z500 as a seasonal ANOMALY (remove doy climatology EXACTLY as for T2m), then detrend ----------------
# Same fixed 1995-2025 +/-3 d doy climatology used for T2m, applied to the ridge-box Z500 (gpm).
zbase = zser[(zser.index.year >= CLIM_Y[0]) & (zser.index.year <= CLIM_Y[1])]
zbdoy = zbase.index.dayofyear.values; zbval = zbase.values
zclim = np.full(367, np.nan)
for D in range(1, 367):
    d = np.abs(((zbdoy - D + 182) % 365) - 182); mm = d <= CLIM_HALF
    if mm.any(): zclim[D] = np.nanmean(zbval[mm])
z_anom = pd.Series(zser.values - zclim[zser.index.dayofyear.values], index=zser.index)   # Z500 anomaly [gpm]
fitza = seaz & (zser.index.year <= FIT_LAST) & np.isfinite(z_anom.values) & np.isfinite(rm_on_z.values)
Z500ANOM_SLOPE = float(np.polyfit(rm_on_z.values[fitza], z_anom.values[fitza], 1)[0])     # gpm per K (on the anomaly)
Z500ANOM_SLOPE_CI = _boot_slope_ci(rm_on_z.values[fitza], z_anom.values[fitza])
z_anom_detr = z_anom + Z500ANOM_SLOPE * (RMST_EVENT - rm_on_z.values)                     # push anomaly to present climate
print(f"      Z500ANOM_DETREND_SLOPE = {Z500ANOM_SLOPE:+.2f} gpm/K  95% CI [{Z500ANOM_SLOPE_CI[1]:+.2f}, {Z500ANOM_SLOPE_CI[2]:+.2f}]")

# ---------------- learned ridge control: T2m_anom_detr ~ Z500_detr (AMJJA 1940-2025) ----------------
print("[4/5] regression ...")
df = pd.DataFrame({"t": anom_detr, "z": z_detr}).dropna()
df = df[(df.index.month >= MONTHS[0]) & (df.index.month <= MONTHS[1]) & (df.index.year <= FIT_LAST)]
BETA, INTERCEPT = np.polyfit(df.z.values, df.t.values, 1)
BETA_CI = _boot_slope_ci(df.z.values, df.t.values)     # (median, lo, hi)
r, p = stats.pearsonr(df.z.values, df.t.values); R2 = r**2
n = len(df)
# event point (2026), detrended the same way
z_event = float(z_detr.reindex([EVENT]).iloc[0]) if EVENT in z_detr.index else np.nan
t_event = float(anom_detr.reindex([EVENT]).iloc[0]) if EVENT in anom_detr.index else np.nan
print(f"      BETA = {BETA:+.4f} K/gpm  95% CI [{BETA_CI[1]:+.4f}, {BETA_CI[2]:+.4f}]  intercept = {INTERCEPT:+.2f}  r = {r:+.3f} (R2={R2:.2f})  n = {n}")

# parallel link on the Z500 ANOMALY axis (seasonal cycle removed like T2m)
dfa = pd.DataFrame({"t": anom_detr, "z": z_anom_detr}).dropna()
dfa = dfa[(dfa.index.month >= MONTHS[0]) & (dfa.index.month <= MONTHS[1]) & (dfa.index.year <= FIT_LAST)]
BETA_A, INTERCEPT_A = np.polyfit(dfa.z.values, dfa.t.values, 1)
BETA_A_CI = _boot_slope_ci(dfa.z.values, dfa.t.values)
rA, pA = stats.pearsonr(dfa.z.values, dfa.t.values); R2A = rA**2
za_event = float(z_anom_detr.reindex([EVENT]).iloc[0]) if EVENT in z_anom_detr.index else np.nan
print(f"      BETA(anom) = {BETA_A:+.4f} K/gpm  95% CI [{BETA_A_CI[1]:+.4f}, {BETA_A_CI[2]:+.4f}]  r = {rA:+.3f} (R2={R2A:.2f})  event Z500' = {za_event:+.0f} gpm")

out = dict(beta_K_per_gpm=BETA, beta_ci95_K_per_gpm=[BETA_CI[1], BETA_CI[2]], intercept_K=INTERCEPT,
           r=r, r2=R2, n=int(n), p_value=float(p),
           rmst_anomaly_slope_K_per_K=RMST_ANOM_SLOPE, rmst_anomaly_slope_ci95_K_per_K=[RMST_ANOM_SLOPE_CI[1], RMST_ANOM_SLOPE_CI[2]],
           z500_detrend_slope_gpm_per_K=Z500_SLOPE, z500_detrend_slope_ci95_gpm_per_K=[Z500_SLOPE_CI[1], Z500_SLOPE_CI[2]],
           z500anom_detrend_slope_gpm_per_K=Z500ANOM_SLOPE, z500anom_detrend_slope_ci95_gpm_per_K=[Z500ANOM_SLOPE_CI[1], Z500ANOM_SLOPE_CI[2]],
           beta_anom_K_per_gpm=BETA_A, beta_anom_ci95_K_per_gpm=[BETA_A_CI[1], BETA_A_CI[2]], r_anom=rA, r2_anom=R2A,
           intercept_anom_K=float(INTERCEPT_A),                                      # ERA5 climatological alpha (anomaly fit): use the SAME alpha everywhere
           t2m_clim_doy_event_C=float(clim[EVENT.dayofyear]),                        # 25-May 1995-2025 doy clim (frac box mean)
           z500_clim_doy_event_m=float(zclim[EVENT.dayofyear]),                      # 25-May 1995-2025 z500 doy clim (frac box mean)
           event_point_anom=dict(z500_anom_detr_gpm=za_event, t2m_anom_detr_C=t_event),
           bootstrap=dict(B=BOOT_B, seed=BOOT_SEED, method="nonparametric resample of fit days; 2.5/97.5 pct"),
           rmst_event_C=RMST_EVENT, event_ref=str(EVENT_REF.date()),
           target_box=TARGET, ridge_box=RIDGE_BOX, rmst_box=RMST_BOX,
           clim_years=list(CLIM_Y), clim_half_window_days=CLIM_HALF, season=SEASON, months=list(MONTHS),
           fit_last_year=FIT_LAST, reference="present (event-year) climate on both axes",
           event_point=dict(z500_detr_gpm=z_event, t2m_anom_detr_C=t_event),
           note=("Ridge-controlled T2m anomaly of a case = t2m_anom_detr - (beta*z500_detr + intercept). "
                 "Detrend a case's z500 to present climate with z500_detr = z500_gpm + z500_detrend_slope*(rmst_event - rmst_year)."))
os.makedirs("${DERIVED_DATA}", exist_ok=True)
JOUT = f"${DERIVED_DATA}/ridge_control_t2m_z500_{SEASON}.json"
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
    fig.savefig(f"${HOME_WORK}/supplement_ridge_control_t2m_z500_{SEASON}.{ext}", dpi=200, bbox_inches="tight")
print(f"saved supplement_ridge_control_t2m_z500_{SEASON}.png/.pdf")

# ---------------- density scatter: Z500 ANOMALY axis (seasonal cycle removed like T2m) ----------------
fig2, ax2 = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
hb2 = ax2.hexbin(dfa.z.values, dfa.t.values, gridsize=70, mincnt=1, bins="log", cmap="magma_r")
xx2 = np.linspace(dfa.z.min(), dfa.z.max(), 100)
ax2.plot(xx2, BETA_A*xx2 + INTERCEPT_A, "k-", lw=2.2, label=f"OLS: {BETA_A:+.3f} K per gpm")
if np.isfinite(za_event):
    ax2.scatter([za_event], [t_event], s=90, c="#1f77b4", edgecolor="k", zorder=6, label=f"25 May 2026 ({t_event:+.1f} °C)")
ax2.axhline(0, color="0.6", lw=0.8, ls=":"); ax2.axvline(0, color="0.6", lw=0.8, ls=":")
ax2.set_xlabel("Detrended ridge $Z_{500}$ anomaly  [gpm, present climate]")
ax2.set_ylabel("Detrended target-box T2m anomaly  [°C]")
ax2.set_title(f"ERA5 ridge control, $Z_{{500}}$ anomaly ({SEASON} 1940–2025): $r$={rA:+.2f}, $R^2$={R2A:.2f}, $n$={len(dfa)}", fontsize=13)
ax2.legend(loc="upper left", frameon=True, fontsize=11)
fig2.colorbar(hb2, ax=ax2, label="days per bin (log)")
for ext in ("png", "pdf"):
    fig2.savefig(f"${HOME_WORK}/supplement_ridge_control_t2m_z500anom_{SEASON}.{ext}", dpi=200, bbox_inches="tight")
print(f"saved supplement_ridge_control_t2m_z500anom_{SEASON}.png/.pdf")
