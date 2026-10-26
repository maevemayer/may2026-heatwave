#!/bin/csh

# LAGRANTO 15-day backward trajectory for ONE start time of ONE analogue.
#
# Trajectory math identical to the notebook / event workflow: same start-file
# box/profile, same 360 h backward caltra, same trace/tracevars. Only the
# wrapper is parametrized so each analogue reads its own P/S files and writes
# its own results:
#
#   Usage: run_lagranto_analogue.sh YYYYMMDD_HH SOURCE_DIR RESULTS_DIR TMP_BASE
#
# NOTE (6-hourly correction): set the datelist interval to the 6-hourly data
# resolution once in the LAGRANTO install:
#   .../lagranto.era5.test/lagranto.shellscript/goodies/datelist.sh
#   set interval = "6"
# The S-file stage now produces P/S files directly on the 00,06,12,18 schedule
# (no rename), so no filename shifting is needed.

if ($#argv != 4) then
    echo "Usage: $0 YYYYMMDD_HH SOURCE_DIR RESULTS_DIR TMP_BASE"
    exit 1
endif

set datetime    = $argv[1]
set source_dir  = $argv[2]
set results_dir = $argv[3]
set tmp_base    = $argv[4]
echo "Processing $datetime ..."

set tmpdir = ${tmp_base}/tmp_$datetime
mkdir -p $tmpdir

\cp ${LAGRANTO_BASE}/create_startf.param $tmpdir/
if (! $?TRACEVARS_FILE) then
    echo "TRACEVARS_FILE not set (04_lagranto_slurm_6h.py sets it)"; exit 1
endif
\cp ${TRACEVARS_FILE} $tmpdir/tracevars.shellscript   # ../tracevars_tc10hpa.shellscript
cd $tmpdir

module purge
module load python3
module load intel-oneapi-compilers
module load netcdf-c/4.8.1-intel-oneapi-mpi-2021.5.0-intel-2021.5.0
module load netcdf-fortran/4.5.3-intel-oneapi-mpi-2021.5.0-intel-2021.5.0

if (! $?LD_LIBRARY_PATH) setenv LD_LIBRARY_PATH ""
setenv LD_LIBRARY_PATH /path/to/hpc-modules/netcdf-fortran-4.5.3-r5r3ev/lib:$LD_LIBRARY_PATH
echo "netcdf fortran library linked to $LD_LIBRARY_PATH"

setenv MODEL era5
setenv NETCDF_FORMAT CF
setenv LAGRANTO ${WORK}/01_LAGRANTO/lagranto.era5.test/lagranto.shellscript

\ln -sf tracevars.shellscript tracevars

set backward_start = `python3 -c "from datetime import datetime, timedelta; dt=datetime.strptime('${datetime}','%Y%m%d_%H'); print((dt - timedelta(days=15)).strftime('%Y%m%d_%H'))"`

foreach date (`${LAGRANTO}/goodies/datelist.sh stdout -create $backward_start $datetime`)
    ln -sf ${source_dir}/P${date} ${tmpdir}/P${date}
    ln -sf ${source_dir}/S${date} ${tmpdir}/S${date}
end

set start_file = "${datetime}_startf.2"
set lonlat_file = "start_lonlat_${datetime}.txt"

cat >! make_lonlat.py << EOF
import sys

lons = [354.5 + 0.5*i for i in range(11)] + [0.0 + 0.5*i for i in range(11)]
lats = [44.5 + 0.5*i for i in range(22)]

with open(sys.argv[1], "w") as f:
    for lat in lats:
        for lon in lons:
            f.write("{:.1f} {:.1f}\n".format(lon, lat))
EOF

python3 make_lonlat.py "$lonlat_file"
if ($status != 0) exit 1

${LAGRANTO}/create_startf/create_startf.sh $datetime $start_file \
"file(${lonlat_file}) @profile(1000,200,17) @hPa" -changet
if ($status != 0) exit 1

set traj_file = "trajectory_${datetime}"

${LAGRANTO}/caltra/caltra.sh $datetime $backward_start $start_file $traj_file -j
if ($status != 0) exit 1

${LAGRANTO}/trace/trace.sh $traj_file $traj_file
if ($status != 0) exit 1

mkdir -p ${results_dir}
mv -n $traj_file ${results_dir}/

exit 0
