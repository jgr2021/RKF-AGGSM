# Current manuscript experiments

This directory contains the MATLAB implementations and frozen evaluation inputs
for the current manuscript **Robust Kalman Filters for Additive Gaussian–GSM
Distributions**. The root `ToyModel/` and `UrbanNav_Experiment/` directories are
retained as earlier experiments; their saved results are not the current tables.

## Run

Tested with MATLAB R2024b on Windows. Install Statistics and Machine Learning
Toolbox for functions used by the filters and calibration routines. No vendor
code download, Python environment, or absolute local path is needed for replay.

From the repository root:

```matlab
cd journal
validate_release          % 3 trajectories x 32 conditions x 14 methods + UrbanNav
run_synthetic(1000)       % all 1,000 trajectories in every manuscript condition
run_synthetic(1000,17:24) % Gaussian-mixture conditions only
run_urban                % full Odaiba-Trimble experiment, 14 methods
```

Results go to `generated/`. These commands do not change the archived inputs.
The full simulation replay takes substantially longer than the validation run.
Filter failure raises an error; failed runs are not silently removed.
The validation compares replayed estimates/errors against archived results.
Small floating-point differences can occur across MATLAB versions/platforms.

## Experiments and parameters

| Condition IDs | Noise family | Conditions | Trajectories | Steps |
|---|---|---:|---:|---:|
| 1–8 | Additive Gaussian + Student's t | 8 | 1,000 each | 30 |
| 9–16 | Additive Gaussian + Slash | 8 | 1,000 each | 30 |
| 17–24 | Gaussian mixture | 8 | 1,000 each | 30 |
| 25–32 | Gaussian + Laplace outliers | 8 | 1,000 each | 30 |

`parameters.json` records the actual noise conditions, method parameters,
calibration estimates, seeds, archive IDs, and reported mean RMSEs.
`data/<family>_C<n>.mat` contains the exact shared observations/truth, process
model, calibration samples, method parameters, and random-number seeds used
for replay. The `*_expected.mat` files retain per-trajectory reference RMSEs.
`results/selected_rmse.csv` is the 32-condition manuscript summary.

The synthetic statistic is **the mean of 1,000 trajectory RMSEs**, each
computed over all 30 steps and both position coordinates. It is not a single
last-step error, a steady-state RMSE, or the square root of a pooled mean.
Synthetic units are normalized; UrbanNav errors are in metres.

The manuscript conditions were selected after exploratory experiments.
`exploration/` retains the broader Student's t/Slash tables, all 12 conditions
of the fresh Gaussian–Laplace evaluation, and the two-route UrbanNav summary.
These selected cases are not an unbiased estimate of performance across all
possible noise laws. Some heavy-tail stress cases have infinite population
second moments and fall outside the finite-second-moment theory assumptions.

For Laplace, the saved settings were chosen per method from four archived
candidate configurations using independent tuning trajectories before the
1,000-trajectory evaluation. Other families use archived scenario-specific
settings. The tuning budgets are not identical across all 14 methods.
Replay uses those saved settings; it does not rerun the historical parameter
searches or regenerate raw GNSS measurements.

## Included methods

KF, PF, Huber KF, RKF-Gaussian, RKF-ST, RKF-GSTM, RKF-TAU, GSM-filter,
GGM, URKF, RKF-SGalphaS, NBKF, **RKF-AGST**, and **RKF-AGSlash**.

The proposed synthetic implementations retain their historical filenames:
`synthetic/filters/test_demoRKFAGSMG.m` is RKF-AGST;
`synthetic/filters/test_demoRKFAGSlash.m` is RKF-AGSlash.
`synthetic/call_method.m` documents the parameter-to-method mapping.
All core filter files are copied from the versions used for the archived
experiments; this release reorganizes their dependencies without changing
the recursions. `source_manifest.json` records their source-relative paths
and SHA-256 hashes. Those paths are provenance, not runtime dependencies.

Notable implementation conventions:

- RKF-TAU uses the corrected covariance recursion and saved `(tau,c)` settings;
  `c=0` is the nominal KF limit.
- GSM-filter uses a residual window of 30 and an iteration cap of 20.
- GGM uses the paper-recursion implementation; there is one GGM result column.
- URKF uses its measurement-uncertainty formulation with `c=0.5`.
- RKF-SGalphaS uses ECF calibration and a gamma-series/importance-sampling
  precision update. NBKF uses the saved nominal-covariance calibration and
  the outlier-probability gate in `competing_options`/`call_method`.

## UrbanNav scope

The manuscript experiment is the **Odaiba–Trimble** sequence. The processed
inputs include prediction-only epochs without a GNSS solution; all epochs
after the first 50 are scored on the same evaluation interval.
The experiment uses reference-derived velocity controls and calibration from
the earlier processed Odaiba route. It is an offline filter comparison, not
an independent online GNSS/INS deployment test. Calibration is included in
`urban_configuration.mat`, and the exact scored estimates are in
`urban_expected.mat`. The documented preprocessing used RTKLIB single-point
GPS+GLONASS with a 15-degree elevation mask (`urban/spp.conf`).

The original dataset is provided by
[UrbanNav](https://github.com/IPNL-POLYU/UrbanNavDataset). Follow its data-use
terms; the processed inputs do not replace the original dataset. We retain
only the data needed for this offline replay, not the raw Tokyo archive.
The exploratory two-route summary additionally includes I-MORKF; it is not
one of the 14 methods in the current manuscript replay.

`results/urban_runtime.csv` and `results/synthetic_runtime.mat` retain the
separate timing measurements. Runtime from a fresh replay is machine-dependent.
The synthetic timing benchmark uses the archived GMM/Laplace settings, not
all 32 current accuracy conditions.

## Baseline references

- D. Li et al., “GSM Filter: A Learning Model of Measurement Noise
  Distributions for Bayesian Filtering,” IEEE TAES, 2023.
- H. Zhu et al., “An Adaptive Kalman Filter With Inaccurate Noise Covariances
  in the Presence of Outliers,” IEEE TAC, 2022.
- “An update-resilient Kalman filtering approach” (URKF).
- P. Hao, O. Karakus and A. Achim, “Robust Kalman Filters Based on the
  Sub-Gaussian alpha-stable Distribution,” Signal Processing 224, 109574, 2024.
- Jia, Huang and Leung, “A Novel Robust Kalman Filter Based on
  Normal–Bernoulli Distribution for Non-Stationary Heavy-Tailed Measurement
  Noise,” IEEE TSP 73, 2025.

The comparison implementations and their adaptations are supplied for this
experiment. The repository does not claim they are official author releases.
Existing attribution in source files is preserved; no third-party vendor
source archive is republished in this directory.

## Condition-axis figure

Optional Python plot of the archived manuscript RMSE values:

```sh
python -m pip install -r requirements-plot.txt
python plot_conditions.py
```

This retains the manuscript's plotting styles and writes the figure under
`generated/`. It plots archived means; it does not rerun the filters.
