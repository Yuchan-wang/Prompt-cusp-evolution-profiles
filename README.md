# Prompt-cusp-evolution-profiles
Code and derived halo/peak density-profile data for prompt-cusp formation and evolution to NFW in 64 WDM zoom-in simulations.

Density profile data
This repository contains radial density profiles and codes in the paper.

## Data
- `LocalRho_L3Nh20_profiles.npy`

This is a NumPy `.npy` file that stores one Python dictionary, containing:
- `profiles`: a dictionary keyed by strings like `L3HaloID{}k{}Nh20S{}`. HaloID is the halo identifier, k is the free-streaming cut-off scale, and S is the snapshot number where a profile can be measured.

Each `profiles[key]` entry is a dictionary with:
- `logr_halo`: 1D numpy array (log radius grid for total density profile)
- `rho_halo`:  1D numpy array (log total density)
- `logr_semi`: 1D numpy array (log radius grid for peak density profile)
- `rho_semi`:  1D numpy array (log peak density)
- `z`:          float (redshift corresponding to snapshot `S`)

k list : [7, 10, 12, 20, 24, 30, 35, 50]

Halo list : [108, 150, 206, 225, 260, 260, 305, 337, 443]

S ranges from 0 up to maximum 140, depending on the halo and k. Redshift decreases with increasing S.

- `meta`: metadata such as the halo/k lists and notes

## Warm Dark Matter Halo Analysis Code

This repository contains the Python analysis pipeline used to investigate the evolution of prompt cusps in Warm Dark Matter (WDM) halos.

The code processes cosmological simulation data (HDF5 format), identifies Lagrangian density peaks, traces particle histories, and fits density profiles to characterize the central structure of dark matter halos.

The analysis is split into several modules, each handling a specific stage of the pipeline.

### 1. `WarmHalo_Peak_Identifier.py`
* **Purpose:** Identifies local density maxima (peaks) in the linear density field of the initial conditions.
* **Method:** Scans the Lagrangian space of the simulation particles to find local maxima in the linear overdensity field ($\delta$).
* **Output:** Stores the locations and properties of these peaks.

### 2. `WarmHalo_PreAnalysis_Subhalos.py`
* **Purpose:** Pre-processes subhalo data and generates density profiles for progenitors.
* **Function:** Loads merger tree information and calculates the spherically averaged density profiles for halos and subhalos across snapshots.

### 3. `WarmHalo_Analysis_Halo_Updated.py`
* **Purpose:** The core analysis script for the main halo and its central density structure.
* **Peak Selection:** Selects "peak particles" using a power-law envelope fit to the density-radius relation in the initial conditions (Lagrangian space).
* **Profile Fitting:** Fits the evolved halo density profiles using multiple models (NFW, Power-law, Broken Power-law) to determine if a central cusp persists.
* **Scanning:** Iteratively scans the fitting radius to find the optimal region describing the central cusp.

### 4. `WarmHalo_Analysis_Subhalos.py`
* **Purpose:** Analyzes the relationship between subhalos and the initial density peaks.
* **Function:** Matches identified subhalos in the evolved simulation to the Lagrangian peaks identified by `WarmHalo_Peak_Identifier.py`.

### 5. `WarmHalo_Analysis_ParticleHistory.py`
* **Purpose:** Traces the accretion history of particles.
* **Function:** Tracks when specific particles (e.g., those belonging to a peak) first fell into the main halo or its progenitors.

### 6. `General_functions.py`
* **Purpose:** A utility library containing shared mathematical and physical functions used across the analysis pipeline.

### 7. `SubHalo_Prog.py`

* **Purpose:** Computes the maximum main-progenitor mass for subhaloes.
---

### Requirements

The code requires a standard scientific Python environment. Key dependencies include:

* `numpy`
* `scipy`
* `matplotlib`
* `h5py` (for reading simulation snapshots)
* `pandas`

---

### Usage

The scripts are designed to be run from the command line, accepting arguments to specify the simulation run and halo parameters.

```bash
python <script_name>.py <Level> <HaloID> <k_index> <T_min> <T_max> <Nh> <Nf> <Nfile>
```
Arguments:

`Level`: Simulation resolution level.

`HaloID`: ID of the target halo (e.g., 0).

`k_index`: Index selecting the free-streaming scale ($k_{fs}$) from the internal list (controls WDM temperature).

`T_min`, `T_max`: Range of snapshots (time steps) to analyze.

`Nh`: Halo particle number parameter (used in path selection).

`Nf`: File number parameter (used in path selection).

`Nfile`: Location specifier of the snapshot files.
