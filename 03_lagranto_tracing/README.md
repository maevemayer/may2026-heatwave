# Corrected 6-hourly LAGRANTO pipeline

SLURM-submittable version of the notebook workflow, with the **6-hourly
snapshot** correction. All the physics/math is carried over verbatim from the
notebooks; only the items below change.

## Genesis rule (Oct 2026): as in Hotz et al. (2024)

Stages 05 and 06 integrate the processes from the TA genesis to arrival with the
genesis of the colleague's reference decomposition: the LAST point back from
arrival that still has the arrival sign of TA. If TA never changes sign within
15 days, the whole trajectory is integrated (her age = 360 h). Before, the sum
ran one 6-h interval past the FIRST opposite-sign point (two steps beyond her
genesis), and trajectories without a sign change (~6 %) got only the last 6 h,
which dropped 15 days of diabatic cooling for those air masses (+0.7 K diabatic
bias at 2026-05-27 12 UTC). With the new rule, our formulas on her trajectories
reproduce her terms (r = 1.00), and our 27 May means match hers to 0.01 K
(diabatic), 0.1 K (adiabatic) and 0.25 K (advective, trajectory differences).
Genesis at arrival (TA flips within the first 6 h) gives 0 K. Only stages 05-06
change: rerun with `STAGES="05 06"`. Old code: archive/corrected_lagranto_pipeline_6h_genesis_firstopp.

## Wind fix (Oct 2026): U, V divided by cos(latitude) - `_windfix`

ERA5 spectral u/v (params 131/132 in the DKRZ pool) are stored as U*cos(lat),
V*cos(lat); `cdo sp2gp` does not divide that out, so all earlier P-files had
horizontal winds too weak by cos(lat) (0.71 at 45N, 0.42 at 65N). Found by
comparing P20260527_12 with an independent MARS 0.5 deg P-file (ratio = cos(lat)
at every latitude; T and OMEGA identical). Stage 01 now divides U, V by cos(lat)
after the remap (verified: corr 1.000, slope 0.999, RMS 0.06-0.12 m/s against
the MARS winds) and flags corrected files; uncorrected wind intermediates are
recomputed. All trajectories/process/origin results from before are affected.

Rerun = new P-files only: `STAGES="01 04 05 06"` (default) builds P with stage 01,
links the existing, correct S-files (no U/V in them) from `S_REUSE_DIR` /
`ANALOGUE_S_REUSE_BASE` (analogues found by date), and runs 04-06 into new
`_windfix` trees. Full run: `STAGES="01 02 03 04 05 06"`.

## Climatology of the 2026 event: fixed 2017-2025 (Oct 2026)

The automatic rule in stage 02 (event year + 8 before, or the previous 9 if the
event year is incomplete) kept 2026 in its own climatology once ERA5T for 2026
was complete: the old run `.../RERUN/corrected_6h/tbar` averages **2018-2026**
(checked from the cdo history of all 140 tbar files). `config.sh` now sets
stage 02 now applies the paper rule (centred 9-year window clamped to
1940-2025, so 2026 -> 2017-2025 as for 2021-2025; Hotz et al. 2024; `CLIM_YEARS`
in `config.sh` overrides, years after 2025 are refused) and writes to `.../RERUN/corrected_6h_tc10hpa_clim2017_2025`. The
climatology changes TCLIM/DTDT in the S-files, so this event run is a full
run (stages 01-06, `PS_DATA_DIR=${DATA_DIR}`). Analogues already used this rule.
Check any finished climatology with `python check_climatology_years.py <tbar dir>`.

## Method update (Oct 2026): d(T_clim)/dp along the trajectory

The vertical gradient of the climatology is no longer the precomputed S-file
field `DTDP` (central difference between neighbouring model levels). As in
Hotz et al. (2024) / Roethlisberger & Papritz (2023), LAGRANTO traces T_clim
10 hPa below and above the parcel (`TCLIM:+10HPA`, `TCLIM:-10HPA` in
`tracevars_tc10hpa.shellscript`, next to `config.sh`), and stage 05 computes
[T_clim(p+10 hPa) - T_clim(p-10 hPa)] / 20 hPa (one-sided near the ground). It
enters advection and adiabatic with opposite signs, so it moves warming between
the two without changing their sum or TA. S-files no longer contain DTDP (~1/3
smaller). Stage 05 stops if the gradient is not positive in 300-900 hPa.

**Reusing the existing P/S files (default).** The new tracevars only need THETA,
TCLIM and DTDT from the S-files, which the old S-files contain (their extra DTDP
field is simply not traced). So stage 04 READS P/S from the old trees
(`PS_DATA_DIR` in `config.sh`, `ANALOGUE_PS_BASE` in `config_analogues.sh`) and
only stages 04-06 run, writing trajectories, decomposition and origin files into
the new `_tc10hpa` trees (`START_STAGE` defaults to 4; the submit scripts refuse
`START_STAGE<4` while P/S are read from elsewhere). For a full new run set
`PS_DATA_DIR=${DATA_DIR}` / `ANALOGUE_PS_BASE=${ANALOGUE_FINAL_BASE}`; stage 03
then writes S-files without DTDP.

All output goes to NEW `_tc10hpa` folders so the old results stay untouched:
event `.../RERUN/corrected_6h_tc10hpa`, analogues `..._tc10hpa` for the
origin, P/S/trajectory and scratch trees. The shared
`${WORK}/01_LAGRANTO/tracevars.shellscript` is no longer used and
left as it was. The pre-change pipeline is in
`./scripts/archive/corrected_lagranto_pipeline_6h_DTDPfield/`.

## What changed (and only this)

1. **No 3-hourly means.** `cdo timselmean(3)` is removed from stages 01 and 02.
2. **6-hourly instantaneous snapshots** (00/06/12/18 UTC) are selected with
   `-selhour,0,6,12,18` **before** the spectral transform (`sp2gp`) — fewer
   timesteps through the expensive step, so it's also the main speed-up.
3. **No rename step.** The old `04_02_rename` corrected a +1 h offset that came
   from *centred* 3h-bin means. Snapshots are instantaneous → times are already
   exactly 00/06/12/18 → nothing to rename. (This is also why the separate
   TCLIM +1 h fix is not needed here.)
4. **DTDT stencil = 6 h.** In stage 03 the climatology time derivative uses the
   ±6 h neighbours (`dt_seconds = 2*6*3600`) because the climatology now lives
   on a 6 h grid. Everything else in 03 (THETA, DTDP, PS, assembly) is identical.
5. **Domain cut to the Northern Hemisphere from the start.** Stages 01 and 02
   apply `-sellonlatbox,$DOMAIN_LATBOX` (default `0,360,0,90`) **chained onto the
   remap**, so ERA5 is trimmed the instant it is regridded and *every* downstream
   file (p, dp, P-files, tbar) is NH-only — less scratch I/O and RAM, faster
   `sp2gp`/remap chains. Set `DOMAIN_LATBOX=""` for the old global behaviour, or
   widen the box (e.g. `0,360,-15,90`) for an equator buffer. 15-day back
   trajectories from W Europe stay well inside the NH; any parcel that leaves the
   box is flagged `-999` by LAGRANTO and masked in stage 06.

Nothing else — grid (`r720x360`), 9-year selection, 21-day running mean,
central differences, potential-temperature formula, variable renames — is
untouched.

## Arrival location for land masking

Stage 06 (event **and** analogue) already writes the **exact arrival point**
(`lon0`, `lat0`, `p0`/`arrival_p`) per trajectory alongside TA and every process
contribution, so a land/ocean mask on arrival is a straight
`is_land(lon0, lat0)` filter on `origin_contributions.nc` — for both the ERA5
event and each analogue rank. (The earlier analogue pooled files that lacked
arrival coordinates came from the *old* framework; files produced by this
pipeline carry them.)

## Files

| file | stage |
|---|---|
| `config.sh` | shared paths / cadence / `LEVELS` knob (source this) |
| `01_prepare_era5_6h.py` + `.slurm` | ERA5 → P-files |
| `02_t_climatology_6h.py` + `.slurm` | 9-year T climatology (tbar) |
| `03_build_s_files_6h.py` + `.slurm` | build S-files (hyam/hybm folded in) |
| `04_run_parallel.csh` + `04_run_lagranto_driver.py` + `.slurm` | LAGRANTO |
| `05_traj_postprocess_6h.py` + `.slurm` | trajectory process decomposition |
| `06_origin_contributions_6h.py` + `.slurm` | origin contributions (TA + per-process genesis + exact arrival) |
| `submit_all.sh` | chain 01+02 → 03 → 04 → 05 → 06 with dependencies |

There is **no stage 04**: the old `04_write_coords.sh` (`ncks -A -v hyam,hybm P S`)
is folded into stage 03 — the coefficients are written straight into each S file
in xarray (no second NCO pass, no `module load nco`) — and the `04_02` rename is
obsolete because the 6-hourly snapshots already carry their true 00/06/12/18
timestamps.

Outputs go to a **new** tree (`${LAGRANTO_BASE}/RERUN/corrected_6h/`) so the old
3-hourly data is never overwritten. Edit paths in `config.sh` if needed.

## Prerequisites (one-time, manual)

- Set `set interval = "6"` in
  `$LAGRANTO/goodies/datelist.sh` (was `"3"`) before stage 05.

## Run the ERA5T 2026 event

```bash
chmod +x *.sh
./submit_all.sh 2026-04-30 2026-06-01 2026-05-15 2026-06-01
#               DATA_START  DATA_END  TRAJ_START  TRAJ_END
```
`DATA_*` must cover `[TRAJ_START − 15 days, TRAJ_END]` (15-day back trajectories).

## Run the 30 best analogues

Each analogue has its own dates/climatology years. Loop over your rank table,
setting the four dates per rank and its own output subdir. Minimal pattern:

```bash
while read RANK DATA_START DATA_END TRAJ_START TRAJ_END; do
  export BASE_OUT="${LAGRANTO_BASE}/RERUN/corrected_6h/rank${RANK}"
  ./submit_all.sh "$DATA_START" "$DATA_END" "$TRAJ_START" "$TRAJ_END"
done < ranks.txt
```
(`config.sh` derives all subpaths from `BASE_OUT`, so exporting it per rank
keeps each analogue isolated. The climatology year auto-selection in stage 02 is
unchanged — it uses the analogue's own base year, or the previous 9 if the base
year is incomplete.)

## Optional vertical restriction — `LEVELS`

`config.sh` has `LEVELS=""` (all 137 model levels → math identical to the
notebooks). Setting e.g. `LEVELS="70/137"` keeps the lowest 68 levels
(~250 hPa → surface), ~50% less 3D data and a matching regrid speed-up.

**Caveat for this study:** `create_startf` seeds trajectories on
`@profile(1000,200,17)` — up to **200 hPa**, above level 70 (~250 hPa) — and
15-day back trajectories often ascend into the UTLS. `70/137` would clip the
200 hPa starts and any air above ~250 hPa. Enable a level cut only after
confirming your trajectories stay below the chosen top (and raise the top, e.g.
~level 55–60 ≈ 150–170 hPa, if not). When `LEVELS` is set it is applied to the
4-D fields only; the p/dp interface-coefficient handling should be validated.

## Rough cost (shared partition, DKRZ node-hours)

Anchored to your logs (`pfiles-rank17` 0.64; `tclim-rank20` 1.9 for 6/22 dates
→ ~3.0 fresh full; `lagranto-rank3` 0.22):

- **~1.1–1.4 node-h per analogue** with the hour-selection (this pipeline).
- **30 analogues + event ≈ ~35–45 node-h**; budget ~50 with margin.
- Without hour-selection (regrid all 24 h) it would be ~80–115 node-h.
