# README: Snake-Swimming-Vortex-Code

This repository contains the MATLAB code to reproduce the vortex figures and calculations in Gregorio et al 2026 (https://doi.org/10.1103/c57h-kx57).

## Citation

If you use this code, please cite:

Gregorio, E., Godoy-Diana, R., & Herrel, A. (2026). Turning without fins: quantifying the distinct kinematics and vortex dynamics of maneuvering swimming snakes. Physical Review E. https://doi.org/10.1103/c57h-kx57

## Data Repository

You can find the data repository associated with this code here:

Gregorio, Elizabeth; Godoy-Diana, Ramiro; Herrel, Anthony, 2026, "Swimming kinematics and volumetric wake measurements for Natrix maura and Nerodia rhombifer", https://doi.org/10.48579/PRO/5Q27ST, data.InDoRES

## Contact

Elizabeth Gregorio: elizabeth.gregorio@espci.fr

## Requirements

- MATLAB (developed/tested in R20XX; requires R2016b or later for local
  functions in script files)
- Image Processing Toolbox (`imread`, image handling)
- `brewermap`, `cmocean`, 'fullfig' are online MATLAB resources and should 
be installed separately if you want to reproduce figures.

## Files

- **`piv_functions.m`** — A MATLAB `classdef` with static methods
  bundling helper functions used by the analysis script:
  - `HydroParam_imp` — computes vorticity, Q-criterion, and related
    flow quantities from the raw velocity field
  - `isosurfacecolor` — computes RGB coloring for vortex isosurface
    rendering

  Called from the analysis script using dot notation, e.g.
  `piv_functions.LoadDat()`.

- **`PIV_PLOTS.m`** — Main analysis script. Loads one
  trial's data (set the `trial` variable at the top to `'turn1'`,
  `'fwd2'`, or `'fwd6'`), computes vortex impulse over time, fits
  polynomials to each vortex event window (from `vortex_range`), 
  generates the impulse and vortex-force figures, and visualization 
  of the vortices.

## How to run

1. Place `piv_functions.m`, `PIV_PLOTS.m`, and both `.mat` data
   files in the same folder.
2. Open the analysis script and set `trial` to one of `'turn1'`,
   `'fwd2'`, or `'fwd6'`.
3. Update the file path at the top of the script to point to
   `ddptv_data.mat` and `piv_functions.m`
4. Run the script.

## Acknowledgements

This code was adapted from the script written by Vincent Stin for his thesis (10.70675/ebec3ce5z7e2dz4ad9z89c4z3f21fdee0281).

## Funding

This work is funded by the Agence Nationale de la Recherche (France) through project DRAGON2 (ANR-20-CE02-0010).

## License

Licensed under the MIT License.
