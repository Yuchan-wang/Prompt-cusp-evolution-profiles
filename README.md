# Prompt-cusp-evolution-profiles
Code and derived halo/peak density-profile data for prompt-cusp formation and evolution to NFW in 64 WDM zoom-in simulations.

Density profile data
This repository contains radial density profiles used to generate the figures in the paper. 

## File Name
- `LocalRho_L3Nh20_profiles.npy`

This is a NumPy `.npy` file that stores one Python dictionary, containing:
- `profiles`: a dictionary keyed by strings like `L3HaloID{}k{}Nh20S{}`. HaloID is the halo identifier, k is the free-streaming cut-off scale, and S is the snapshot number where a profile can be measured.

k list : [7, 10, 12, 20, 24, 30, 35, 50]

Halo list : [108, 150, 206, 225, 260, 260, 305, 337, 443]

S ranges from 0 up to maximum 140, depending on the halo and k. Redshift decreases with increasing S.

- `meta`: metadata such as the halo/k lists and notes

Each `profiles[key]` entry is a dictionary with:
- `logr_halo`: 1D numpy array (log radius grid for total density profile)
- `rho_halo`:  1D numpy array (log total density)
- `logr_semi`: 1D numpy array (log radius grid for peak density profile)
- `rho_semi`:  1D numpy array (log peak density)
- `z`:          float (redshift corresponding to snapshot `S`)
