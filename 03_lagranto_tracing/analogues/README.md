# Analogue pipeline — corrected 6-hourly

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

Per-rank, Slurm-chained port of your analogue framework, with the **same
6-hourly correction** as the event pipeline one level up. One rank per job, in
an `afterok` dependency chain; validation, staging and scratch cleanup are all
kept from your originals. **The analogue dates live in one shared file**
(`analogue_dates.txt`) and **any number of analogues** is supported.

## Stages (run in this order)

| stage | python (6h) | slurm wrapper | submit (chain 1→N) |
|---|---|---|---|
| 01 P-files | `01_prepare_era5_slurm_6h.py` | `01_run_prepare_era5_analogue.slurm` | `01_submit_prepare_era5.sh` |
| 02 climatology | `02_t_clim_slurm_6h.py` | `02_run_t_clim_analogue.slurm` | `02_submit_t_clim.sh` |
| 03 S-files (coords folded in) | `03_s_files_slurm_6h.py` | `03_run_s_files_analogue.slurm` | `03_submit_s_files.sh` |
| 04 LAGRANTO | `04_lagranto_slurm_6h.py` | `04_run_lagranto_analogue.slurm` (+ `04_run_lagranto_analogue.sh`) | `04_submit_lagranto.sh` |
| 05 decomposition | `05_traj_postprocess_slurm_6h.py` | `05_run_traj_postprocess_analogue.slurm` | `05_submit_traj_postprocess.sh` |
| 06 origin contrib | `06_origin_contributions_slurm_6h.py` | `06_run_origin_contributions_analogue.slurm` | `06_submit_origin_contributions.sh` |

Shared, edit-once: `analogue_dates.txt` (the dates) + `analogue_dates.py` (loader).

Consecutively numbered (`01/02/03/04/05`) to match the 2026-event pipeline. The
old `04_write_coords.sh` (hyam/hybm append) and `04_02_rename` are folded into
stage 03 — hyam/hybm are written straight into each S file, and there is **no
hour shift** because the 6-hourly snapshots are true to their 00/06/12/18
timestamps. 01 and 02 are independent; 03 needs both; 04 needs 03; 05 needs 04;
06 needs 05.

## What changed vs your originals (only this)

1. **No 3-hourly means** (`cdo.timselmean(3)` removed in stages 01 and 02).
2. **6-hourly snapshots** selected with `-selhour,0,6,12,18` **before** `sp2gp`.
3. **SLOTS = 00/06/12/18** everywhere; intermediate/final labels `_3h`→`_6h`,
   climatology files are `t_YYYY-MM-DD_6hHH_tbar.nc`.
4. **Rename phase dropped** in stage 03 (`04_02_rename` / `rename_to_synoptic`):
   instantaneous snapshots already carry exact 00/06/12/18 times, so there is no
   +1 h offset and no dual naming schedule. `--skip-rename` is still accepted but
   is a no-op.
5. **DTDT stencil = ±6 h** in stage 03 (`dt_seconds = 2*6*3600`), forced by the
   6 h climatology cadence.

Everything else is byte-for-byte your logic: the **centered 9-year climatology**
(`analogue_year ± 4`, clamped to 1940–2025), the 21-day running mean, hybrid
p/dp, THETA, DTDP, the hyam/hybm fold-in (old `04_write_coords.sh`), the
edge-slot accounting (`total − 2`), and **step 6 is already 6-hourly**
(`delta_t = 6*3600`, loader asserts 61 steps).

## Configuration — `config_analogues.sh` (single source of truth)

All paths, the Slurm account/partition, the CDO/ERA5/LAGRANTO locations and the
default worker counts live in **`config_analogues.sh`** (the analogue twin of
`../config.sh`). Every submit script sources it and every Python stage reads the
same values from the environment, so you edit paths in **one** place.

Key entries (override by exporting before you submit, or editing the file):

| variable | what lands there | default (all NEW for this rerun) |
|---|---|---|
| `ANALOGUE_HOME` | scripts + Slurm logs + origin output (uc1275) | `${LAGRANTO_BASE}/RERUN/corrected_lagranto_pipeline_6h/analogues` |
| `ANALOGUE_ORIGIN_DIR` | **only** `rankNN_MJ_<date>/origin_contributions.nc` (stage 06) | `${ANALOGUE_HOME}` |
| `ANALOGUE_FINAL_BASE` | **P/S files, trajectories (`results_6h`, `traj_netcdfs`), `t_clim/`** | `${LAGRANTO_RERUN}/analogues` |
| `ANALOGUE_LOG_DIR` | Slurm logs | `${ANALOGUE_HOME}/slurm_logs` |
| `ANALOGUE_SCRATCH_BASE` | scratch root (per-stage subdir appended) | `${SCRATCH}/lagranto/RERUN_2026_MAY_6h` |
| `ANALOGUE_PYTHON` | the analogues-2026 kernel python — **edit if it moved** | `.../${ANALOGUE_BASE}/kernel/analogues-2026/bin/python` |
| `CDO_BIN_DIR`, `ERA5_*` | CDO + ERA5 pools (fixed DKRZ) | as event |
| `LAGRANTO_ROOT` | LAGRANTO install | as event |
| `DOMAIN_LATBOX`, `LEVELS`, `HOURS_CSV` | physics knobs (same meaning as event) | `0,360,0,90` / all / `0,6,12,18` |

**Storage split:** the big, persistent trajectory data (P/S, `results_6h`,
`traj_netcdfs`, `t_clim`) stays on **bu1561** at a NEW `RERUN/...` path (so it
never collides with your old 3-hourly run); scratch is a NEW `RERUN_...` path;
and **only** the small stage-06 `origin_contributions.nc` is written to the new
uc1275 data dir. Stage 06 reads its trajectories from `ANALOGUE_FINAL_BASE` and
writes the origin file to `ANALOGUE_ORIGIN_DIR`. The Python defaults equal these
values, so running a stage directly still lands in the right place without the
config.

## Prerequisites

- Point `ANALOGUE_PYTHON` in `config_analogues.sh` at your analogues-2026 kernel.
- Set `set interval = "6"` in `${LAGRANTO_ROOT}/goodies/datelist.sh` before stage 04/05.

## Run all ranks

**Global submit (recommended)** — the analogue twin of `../submit_all.sh`. For
every rank it chains all six stages with the same dependency graph as the event
(01 & 02 → 03 → 04 → 05 → 06):

```bash
cd analogues
./submit_all_analogues.sh                   # all ranks in analogue_dates.txt
START_RANK=1 END_RANK=5 ./submit_all_analogues.sh
SERIAL_RANKS=1 ./submit_all_analogues.sh    # chain rank N after rank N-1 (bounds active scratch)
```

**Per-stage submit (legacy / partial reruns)** — each submits ranks `1 → N` in
its own `afterok` chain across ranks for that single stage:

```bash
./01_submit_prepare_era5.sh      # 01
./02_submit_t_clim.sh            # 02  (independent of 01)
./03_submit_s_files.sh           # 03  (after 01 & 02 finish)
./04_submit_lagranto.sh          # 04  (after 03)
./05_submit_traj_postprocess.sh  # 05  (after 04)
./06_submit_origin_contributions.sh  # 06  (after 05)
```
Override ranks with `START_RANK=1 END_RANK=5 ./...`; `END_RANK` defaults to all.

## Arrival window around the analogue date

Trajectories are started from `date − ANALOGUE_DAYS_BEFORE` to
`date + ANALOGUE_DAYS_AFTER` (6-hourly). The window may be **asymmetric** — set
each side independently in `config_analogues.sh` (default **4 before / 3 after**):

```bash
export ANALOGUE_DAYS_BEFORE="4"
export ANALOGUE_DAYS_AFTER="3"
```

All stages (01–05) read these from the environment, so change them once.
`ANALOGUE_LEAD_DAYS` (default 15) is separate — it is the back-trajectory length
(days of P/S data prepared before the earliest arrival), not the arrival window.
If you change the window, re-run the whole 01→06 chain for the affected ranks.

## Domain cut (Northern Hemisphere)

Stages 01 and 02 cut ERA5 to `DOMAIN_LATBOX` (default `0,360,0,90` = whole NH)
with `-sellonlatbox` chained onto the remap, exactly like the event pipeline, so
every P/tbar file is NH-only. It is inherited from the environment — submit with
`--export=ALL,...,DOMAIN_LATBOX=0,360,-15,90` to widen it, or `DOMAIN_LATBOX=`
(empty) to keep the global grid. Stage 06 already writes `lon0/lat0/arrival_p`
per trajectory, so an arrival land mask works for the analogues too.

## Vertical restriction

`LEVELS` is read from the environment in stages 01 and 02 (default empty = all
137 levels, math identical). To test the colleague's cut, submit with
`--export=ALL,...,LEVELS=70/137` — but see the top-level README caveat: it clips
the 200 hPa trajectory starts and any air above ~250 hPa.

## Analogue dates — one shared file

All five stages read the analogue dates from **`analogue_dates.txt`** via the
`analogue_dates.py` loader (`from analogue_dates import ANALOGUE_DATES`). No
stage embeds its own list anymore — **edit `analogue_dates.txt` once** and every
step (and the submit scripts' rank count) picks it up.

- Format: one ISO date (`YYYY-MM-DD`) per line; blank lines and `#`-comments
  ignored. Rank `N` = the N-th non-comment line (1-based), matching `--rank N`.
- **Any number of analogues** is allowed — the Python range check uses
  `len(ANALOGUE_DATES)` and the submit scripts count the dated lines, so
  `END_RANK` defaults to however many you list.
- The active set is **factual AMJJA "only cut-off"**; the variants (deep-trough,
  counterfactual MJ, cut-offs+troughs, no-distinction) are kept as commented
  blocks in the same file — uncomment one to switch.
- Check what's active: `python3 analogue_dates.py` prints the numbered list and
  the total.
