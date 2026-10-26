#!/bin/csh
# ===========================================================================
# Stage 05 - single-datetime LAGRANTO backward run (15 days).
# Adapted from 05_run_parallel-7.sh. Math/trajectory setup UNCHANGED.
# Only the data source dir and the results dir point at the corrected 6h tree.
#
# REMINDER: set  set interval = "6"  in
#   $LAGRANTO/goodies/datelist.sh   (was "3")
# so LAGRANTO steps on the 6-hourly P/S files.
#
# Usage: csh 04_run_parallel.csh YYYYMMDD_HH
# ===========================================================================
if ($#argv != 1) then
    echo "Usage: $0 YYYYMMDD_HH"
    exit 1
endif

set datetime = $argv[1]
echo "Processing $datetime ..."

set source_dir = ${PS_DATA_DIR}   # P/S files read (config.sh; may be the old run's data)
set results_dir = ${RESULTS_DIR}
set tmpdir = ${DATA_DIR}/tmp/tmp_$datetime   # run folder in the NEW tree, never in PS_DATA_DIR
mkdir -p $tmpdir

\cp ${LAGRANTO_BASE}/create_startf.param $tmpdir/
\cp ${TRACEVARS_FILE} $tmpdir/tracevars.shellscript   # tracevars_tc10hpa.shellscript (config.sh)
cd $tmpdir

module purge
module load python3
module load intel-oneapi-compilers
module load netcdf-c/4.8.1-intel-oneapi-mpi-2021.5.0-intel-2021.5.0
module load netcdf-fortran/4.5.3-intel-oneapi-mpi-2021.5.0-intel-2021.5.0

if (! $?LD_LIBRARY_PATH) setenv LD_LIBRARY_PATH ""
setenv LD_LIBRARY_PATH /path/to/hpc-modules/netcdf-fortran-4.5.3-r5r3ev/lib:$LD_LIBRARY_PATH

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
${LAGRANTO}/trace/trace.sh $traj_file $traj_file

mkdir -p $results_dir
mv -f $traj_file $results_dir/
exit 0
