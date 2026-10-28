import datetime
from collections import defaultdict

import numpy as np
import earthkit.data as ekd
import earthkit.regrid as ekr

from tqdm import tqdm
from ecmwf.opendata import Client as OpendataClient

import copy
import os
import pickle
import multiurl #in some versions of earthkit they forgot to import it -_-
import multiurl.base #in some versions of earthkit they forgot to import it -_-
import argparse

os.environ["CDSAPI_RC"] = os.path.expanduser("~/.cdsapirc")  # your CDS API credentials
PARAM_SFC = ["10u", "10v", "2d", "2t", "msl", "skt", "sp", "tcw", "lsm", "z", "slor", "sdor"]
PARAM_SOIL = ["stl1", "stl2"]
PARAM_PL = ["z", "t", "u", "v", "w", "q"]
LEVELS = [1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 50]

def get_open_data(target_date, param, levelist=[]):
    """Fetches ERA5 reanalysis data from CDS for a target date and 6 hours before.

    Retrieves data on a 0.25-degree grid and interpolates it to the N320 Gaussian
    grid. For each parameter, the two time steps (t-6h and t) are stacked into
    a single array.
    
    Args:
        target_date (datetime.datetime): The target date and time for which to
            retrieve data. Data is also fetched 6 hours before this date.
        param (list[str]): List of parameter short names to retrieve
            (e.g. ``["10u", "10v", "2t"]``).
        levelist (list[int], optional): Pressure levels in hPa to retrieve.
            If non-empty, fetches from the pressure-level dataset; otherwise
            fetches from the single-level dataset. Defaults to ``[]``.
    
    Returns:
        dict[str, np.ndarray]: A mapping from parameter name to a stacked NumPy
            array of shape ``(2, N320_lat, N320_lon)``. Pressure-level parameters
            are keyed as ``"<param>_<level>"`` (e.g. ``"z_500"``), while
            surface parameters use the short name directly (e.g. ``"2t"``).
    """
    fields = defaultdict(list)
    # CDS separates surface and pressure level datasets
    dataset = "reanalysis-era5-pressure-levels" if levelist else "reanalysis-era5-single-levels"
    # Get the data for the current date and the previous date
    for date in [target_date - datetime.timedelta(hours=6), target_date]:
        request = {
            "product_type": "reanalysis",
            "param": param,
            "date": date.strftime("%Y-%m-%d"),
            "time": date.strftime("%H:%M"),
            "grid": [0.25, 0.25],
        }
        
        if levelist:
            request["pressure_level"] = levelist

        data = ekd.from_source("cds", dataset, **request)
        
        for f in data:
            assert f.to_numpy().shape == (721, 1440)
            values = f.to_numpy()
            
            # Interpolate the data from 0.25 to N320
            values = ekr.interpolate(values, {"grid": (0.25, 0.25)}, {"grid": "N320"})
            
            # Add the values to the list
            name = (
                f"{f.metadata('param')}_{f.metadata('levelist')}"
                if levelist
                else f.metadata("param")
            )
            fields[name].append(values)

    # Create a single matrix for each parameter
    for param, values in fields.items():
        fields[param] = np.stack(values)

    return fields

def build_initial_state(target_date):
    """Builds a complete initial model state by fetching all required ERA5 fields.

    Combines surface, soil, and pressure-level parameters into a single state
    dictionary. The fields are fetched at ``target_date`` and 6 hours before,
    and interpolated to the N320 Gaussian grid.
    
    Args:
        target_date (datetime.datetime): The target date and time for the
            initial state.
    
    Returns:
        dict: A state dictionary with the following keys:
    
            - ``"date"`` (datetime.datetime): The target date.
            - ``"fields"`` (dict[str, np.ndarray]): All fetched fields, keyed
              by parameter name (surface/soil) or ``"<param>_<level>"``
              (pressure-level). Each value is an array of shape
              ``(2, N320_lat, N320_lon)``.
    """
    fields = {}
    
    fields.update(get_open_data(target_date, param=PARAM_SFC))
    fields.update(get_open_data(target_date, param=PARAM_SOIL))
    fields.update(get_open_data(target_date, param=PARAM_PL, levelist=LEVELS))

    return dict(date=target_date, fields=fields)

def build_and_save_initial_states(
    init_control_state, date_input, save_path, exp_name):
    """Builds and saves a an initial state from era5

    Fetches the ERA5 data for the day and exports it in the right format for AIFS-Ens v1.
    
    Output files are written to a subdirectory named: <date>-<exp_name>
    
    Args:
        date_input (str): Initial date and time as a string ("%Y-%m-%d %H:%M:%S").
        save_path (str): Root directory under which output files are saved.
        exp_name (str): Name of the experiment for folder and file naming.
    
    Returns:
        None: Exports .pkl files 
    """
    DATE = datetime.datetime.strptime(date_input, "%Y-%m-%d %H:%M:%S")
    date_str = DATE.strftime("%Y%m%d%H")
    output_dir = os.path.join(save_path, f"{date_str}-{exp_name}")
    os.makedirs(output_dir, exist_ok=True)

    #save init state
    with open(
        os.path.join(output_dir, f"init_state-{date_str}-{exp_name}-init_state.pkl"), "wb"
    ) as f:
        pickle.dump(init_control_state, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--date_input", type=str, required=True)
    parser.add_argument("--exp_name", type=str, required=True)
    parser.add_argument("--output_save", type=str, required=True)
    args = parser.parse_args()
        
    print(f"Downloading and interpolating control state for {args.date_input}...")
    DATE = datetime.datetime.strptime(args.date_input, "%Y-%m-%d %H:%M:%S")
    init_control_state = build_initial_state(DATE)
    
    build_and_save_initial_states(
        init_control_state, args.date_input, args.output_save, args.exp_name
    )
    
    print("Done")
