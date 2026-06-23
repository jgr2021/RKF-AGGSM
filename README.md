# Robust Kalman Filters for Additive Gaussian–GSM Distributions

This repository contains the MATLAB code, processed experiment inputs, saved outputs, and figure-generation scripts associated with the study **“Robust Kalman Filters for Additive Gaussian–GSM Distributions.”**

The repository evaluates robust Kalman filtering under additive hybrid measurement noise, where a persistent Gaussian background component coexists with impulsive or heavy-tailed disturbances. It contains two complementary experiment packages:

- `ToyModel/` — Monte Carlo simulations for Gaussian-mixture, Gaussian-plus-Laplacian-outlier, additive Gaussian–Student's t, and additive Gaussian–Slash measurement-noise settings.
- `UrbanNav_Experiment/` — A two-dimensional GNSS/UrbanNav trajectory-estimation experiment with included processed MATLAB data, filter implementations, saved outputs, and plotting scripts.



## Quick start

### Synthetic experiments

```matlab
cd ToyModel
run_GMM_simulation
```

Change `case_id` near the beginning of a driver script to select a scenario. The main drivers are:

- `run_GMM_simulation.m`
- `run_Gaussian_Outliers_simulation.m`
- `run_AGST_simulation.m`
- `run_AGSlash_simulation.m`

### UrbanNav experiment

```matlab
cd UrbanNav_Experiment
run_UrbanNav_full_comparison
```

The main UrbanNav driver loads the processed inputs from `data/`, adds `filters/` to the MATLAB path, computes RMSE/MAE trajectories, and saves the consolidated result file in `results/`.

## Repository structure

```text
.
├── ToyModel/              Synthetic experiments and saved outputs
└── UrbanNav_Experiment/       UrbanNav experiment
    ├── data/                  Processed MATLAB input files
    ├── filters/               Filter implementations
    └── results/               Saved UrbanNav outputs
```

## Data Notice

The urban navigation data used in this repository are derived from the **UrbanNav Dataset**. We gratefully acknowledge the UrbanNav authors and maintainers for publicly releasing this challenging multisensory benchmark for research on robust positioning and navigation in urban environments.

UrbanNav provides GNSS, IMU, LiDAR, camera, and ground-truth measurements collected in challenging urban canyons. The original dataset, documentation, and usage instructions are available from the official UrbanNav GitHub repository:

https://github.com/IPNL-POLYU/UrbanNavDataset

This repository does not claim ownership of the UrbanNav dataset. The included files are used only for reproducing the experiments reported in our work. Users should follow the original UrbanNav data-access conditions, citation requirements, and any applicable redistribution restrictions.


## Citation

Please cite the associated manuscript when using this code. Add the final bibliographic record and DOI here after publication.
