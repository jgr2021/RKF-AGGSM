# Robust Kalman Filters for Additive Gaussian–GSM Distributions

This repository contains the MATLAB code, processed experiment inputs, saved outputs, and figure-generation scripts associated with the study **“Robust Kalman Filters for Additive Gaussian–GSM Distributions.”**

The repository evaluates robust Kalman filtering under additive hybrid measurement noise, where a persistent Gaussian background component coexists with impulsive or heavy-tailed disturbances. It contains two complementary experiment packages:

- `ToyModel_GMM/` — Monte Carlo simulations for Gaussian-mixture, Gaussian-plus-Laplacian-outlier, additive Gaussian–Student's t, and additive Gaussian–Slash measurement-noise settings.
- `UrbanNav_Experiment/` — A two-dimensional GNSS/UrbanNav trajectory-estimation experiment with included processed MATLAB data, filter implementations, saved outputs, and plotting scripts.

## Requirements

- MATLAB. The archive does not record a tested release.
- Statistics and Machine Learning Toolbox, used by functions such as `mvnrnd`, `trnd`, and `randsample`.

Run scripts from their own experiment directories because several scripts use relative file paths.

## Quick start

### Synthetic experiments

```matlab
cd ToyModel_GMM
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
├── ToyModel_GMM/              Synthetic experiments and saved outputs
└── UrbanNav_Experiment/       UrbanNav experiment
    ├── data/                  Processed MATLAB input files
    ├── filters/               Filter implementations
    └── results/               Saved UrbanNav outputs
```

Each directory contains its own README with more detailed usage notes.

## Reproducibility notes

- Synthetic drivers set `rng(2025)` for repeatable Monte Carlo trials.
- Some driver scripts save or overwrite files such as `*_results_*.mat`, `tune_*.mat`, and `*_best_params*.mat` in the current working directory.
- The repository includes precomputed `.mat`, `.eps`, `.png`, `.csv`, and `.xlsx` artifacts used to produce the reported comparisons.
- Some legacy source-file names retain `AGSMG`; these names are preserved for backward compatibility with the original experiment code.

## Data notice

The `UrbanNav_Experiment/data/` directory contains the processed MATLAB inputs used by the included scripts. Before making the repository public, confirm that you are authorized to redistribute every data file and any derivative data under the applicable data-source terms.

## Citation

Please cite the associated manuscript when using this code. Add the final bibliographic record and DOI here after publication.
