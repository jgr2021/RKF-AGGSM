# Robust Kalman Filters for Additive Gaussian–GSM Distributions

This repository contains the MATLAB code, processed experiment inputs, saved outputs, and figure-generation scripts associated with the study **“Robust Kalman Filters for Additive Gaussian–GSM Distributions.”**

The repository evaluates robust Kalman filtering under additive hybrid measurement noise, where a persistent Gaussian background component coexists with impulsive or heavy-tailed disturbances. It contains two complementary experiment packages:

- `ToyModel/` — Monte Carlo simulations for Gaussian-mixture, Gaussian-plus-Laplacian-outlier, additive Gaussian–Student's t, and additive Gaussian–Slash measurement-noise settings.
- `UrbanNav_Experiment/` — A two-dimensional GNSS/UrbanNav trajectory-estimation experiment with included processed MATLAB data, filter implementations, saved outputs, and plotting scripts.

## Background and Research Significance

Real-world measurements are rarely corrupted by only one type of uncertainty. In applications such as urban GNSS navigation, a persistent Gaussian background error remains present while multipath, blockage, and other abnormal events introduce intermittent large deviations.

Rather than approximating the total error by a single Gaussian or heavy-tailed distribution, the proposed **additive Gaussian--Gaussian scale mixture (AGGSM)** framework explicitly models the two coexisting sources as

$$\mathbf{v}_k=\mathbf{v}_k^{\mathrm{B}}+\mathbf{v}_k^{\mathrm{I}},$$

where $\mathbf{v}_k^{\mathrm{B}}$ denotes the persistent Gaussian background component and $\mathbf{v}_k^{\mathrm{I}}$ denotes the impulsive GSM component.

### Why AGGSM?

- **Preserves the physical noise structure:** background uncertainty does not disappear when an impulsive disturbance occurs.
- **Avoids unrealistic covariance collapse:** the variance offset \(d>0\) prevents the latent impulsive scale from approaching zero and retains a nonzero Gaussian background contribution.
- **Enables practical robust filtering:** a deterministic MAP-based latent-variable update yields tractable recursive filters for additive hybrid-noise estimation.

This repository provides reproducible MATLAB implementations of [RKF-AGST](https://github.com/jgr2021/RKF-AGST) and RKF-AGSlash, together with Gaussian-mixture, Gaussian-plus-Laplacian and real-world UrbanNav experiments. The included examples demonstrate the value of explicitly modeling additive background-plus-impulse noise for robust state estimation under challenging measurement conditions.


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

UrbanNav provides GNSS, IMU, LiDAR, camera, and ground-truth measurements collected in challenging urban canyons. The original dataset, documentation, and usage instructions are available from the official [UrbanNav GitHub repository](https://github.com/IPNL-POLYU/UrbanNavDataset).

## Acknowledgement

I would like to express my sincere gratitude to my supervisor, Prof. Xiao-Ping Zhang, for his invaluable research guidance and insightful revision feedback on this manuscript. I would also like to thank Prof. Zhenyu Liu for his valuable discussions, constructive suggestions, and continued support throughout this work.

This work originated from a course project on Bayesian Learning, [RKF-AGST](https://github.com/jgr2021/RKF-AGST), and was subsequently extended into the present study. I also thank Prof. Ercan E. Kuruoğlu and Dr. Pengcheng Hao for their valuable guidance and assistance during and after the course project.



## Citation

Please cite the associated manuscript when using this code. Add the final bibliographic record and DOI here after publication.
