# ToyModel

This directory contains MATLAB simulations for robust Kalman filtering under several synthetic measurement-noise models. The experiments use a two-dimensional constant-velocity tracking model and compare KF, particle-filter, Huber, Gaussian robust-filter, Student's t, Gaussian–Student's t mixture, divergence-based, and additive Gaussian–GSM variants.

## Main experiment drivers

Run these scripts from this directory.

| Script | Noise setting | Scenario selector |
|---|---|---|
| `run_GMM_simulation.m` | Two-component Gaussian mixture model (GMM) | `case_id = 1,2,3,4` |
| `run_Gaussian_Outliers_simulation.m` | Gaussian background plus intermittent Laplacian outliers | `case_id = 1,2,3,4` |
| `run_AGST_simulation.m` | Additive Gaussian plus Student's t noise | `case_id = 1,2,3,4` |
| `run_AGSlash_simulation.m` | Additive Gaussian plus Slash noise | `case_id = 1,2,3,4` |

Each driver sets `rng(2025)`, executes Monte Carlo trials, reports average runtime, saves a case-level `.mat` file, and plots RMSE/MAE trajectories. Change `case_id` near the top of the selected driver before running it.

Example:

```matlab
cd ToyModel
run_GMM_simulation
```

## Tuning scripts

The `tune_*` scripts perform grid searches and save selected parameter files consumed by the main drivers.

- `tune_filter_hyperparameters*.m` — tunes Huber, GSTM, and RKF-TAU baselines.
- `tune_AGGSM_parameters*.m` — tunes additive Gaussian–GSM / AGST parameters.
- `tune_AGSlash_parameters*.m` — tunes additive Gaussian–Slash parameters.

The repository already contains the resulting `.mat` parameter files. Re-running a tuning script can overwrite them.

## Supporting code

- `Traj_*.m` — trajectory and measurement-noise generators.
- `test_demo*.m` — individual filter-evaluation routines.
- `compute_CRLB_*.m` — CRLB calculations for GMM or equivalent Gaussian settings.
- `fit_Gaussian_StudentT.m`, `GMMrnd.m`, and `slashrnd.m` — distribution-fitting and sampling utilities.
- `calculate_loss.m` — RMSE/MAE calculation.

## Figures and saved outputs

The directory includes precomputed artifacts:

- `*_results_*.mat` — saved per-case metrics and runtime outputs.
- `tune_*.mat` and `*_best_params*.mat` — tuning grids and selected hyperparameters.
- `*_RMSE.eps`, `*_MAE_by_case.eps`, and `.png` files — publication-oriented figures.
- `RMSE_traces_Case*.csv` and `AGST_performance_summary.xlsx` — exported numerical summaries.

Use the figure scripts to regenerate summary plots:

- `generate_all_RMSE_figures.m`
- `generate_AGST_RMSE_figures.m`
- `generate_AGSlash_RMSE_figures.m`
- `plot_metric_by_case*.m`

## Notes

- MATLAB with the Statistics and Machine Learning Toolbox is required (`mvnrnd`, `trnd`, and `randsample` are used).
- Several source-file names use the legacy acronym `AGSMG`; they are retained to preserve the original function calls and result-file compatibility.
- Run scripts from this directory, not from the repository root, because parameter and result files are loaded using relative paths.
