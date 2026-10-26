#!/usr/bin/env python3
"""
Stage 05 driver - launch LAGRANTO for every 6-hourly start time in a period.
Adapted from 05_run_lagranto_parallel-6.ipynb. Trajectory setup UNCHANGED;
only the cadence is fixed to 6-hourly (hours = [0,6,12,18]).

Usage:
    python 04_run_lagranto_driver.py START_DATE END_DATE [N_JOBS]
    e.g. python 04_run_lagranto_driver.py 2026-05-10 2026-06-01 84
"""
import os
import sys
import subprocess
from datetime import datetime, timedelta
from joblib import Parallel, delayed

HERE = os.path.dirname(os.path.abspath(__file__))
SHELL_SCRIPT = os.path.join(HERE, "04_run_parallel.csh")
WORKING_DIR = "${LAGRANTO_BASE}/"
HOURS = [0, 6, 12, 18]                      # 6-hourly


def run_hour(dt):
    datetime_str = dt.strftime("%Y%m%d_%H")
    with open(os.devnull, "w") as devnull:
        subprocess.run(["csh", SHELL_SCRIPT, datetime_str],
                       cwd=WORKING_DIR, stdout=devnull, stderr=devnull, check=True)


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: 04_run_lagranto_driver.py START_DATE END_DATE [N_JOBS]")
    start_date = datetime.strptime(sys.argv[1], "%Y-%m-%d")
    end_date = datetime.strptime(sys.argv[2], "%Y-%m-%d")
    n_jobs = int(sys.argv[3]) if len(sys.argv) > 3 else 84

    current = start_date
    datetimes = []
    while current <= end_date:
        for h in HOURS:
            datetimes.append(current.replace(hour=h))
        current += timedelta(days=1)

    # first/last trimmed (as in the notebook) so neighbours exist for tracing
    Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(run_hour)(dt) for dt in datetimes[1:-1]
    )


if __name__ == "__main__":
    main()
