# 06_derived_data — downloaded separately from Zenodo

The derived data for this study (~11 GB: ERA5 daily fields, the IFS Lagrangian
temperature-decomposition trajectories, origin-contribution fields, COL
catalogues, AIFS/IFS CSVs and ridge-control JSONs) is **too large for GitHub**
and is archived on Zenodo:

    Zenodo data archive DOI: 10.5281/zenodo.23666704
    https://doi.org/10.5281/zenodo.23666704

Download the archive, unzip it, and place its contents here so this folder
matches the layout in `../FILE_TREE.txt`. Then `config.yaml` works unchanged
(`DERIVED_DATA` / the `${...}` keys resolve to `./06_derived_data`).
