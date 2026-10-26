#!/usr/bin/env python3
"""Shared analogue-date list for every pipeline stage.

Reads analogue_dates.txt (next to this file) once: one ISO date per line, with
blank lines and #-comments ignored. Every stage does

    from analogue_dates import ANALOGUE_DATES

so the set of analogues - and how many there are - is edited in ONE place.
"""
from pathlib import Path

DATES_FILE = Path(__file__).with_name("analogue_dates.txt")


def load_analogue_dates(path=DATES_FILE):
    dates = []
    with open(path) as handle:
        for line in handle:
            token = line.split("#", 1)[0].strip()
            if token:
                dates.append(token)
    if not dates:
        raise RuntimeError(f"No analogue dates found in {path}")
    return dates


ANALOGUE_DATES = load_analogue_dates()

if __name__ == "__main__":
    for i, d in enumerate(ANALOGUE_DATES, 1):
        print(f"{i:>3}  {d}")
    print(f"total: {len(ANALOGUE_DATES)}")
