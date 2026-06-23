# Robust Kalman Filters for Additive Gaussian–GSM Distributions

This repository contains the MATLAB code, processed experiment inputs, saved outputs, and figure-generation scripts associated with the study **“Robust Kalman Filters for Additive Gaussian–GSM Distributions.”**

The repository evaluates robust Kalman filtering under additive hybrid measurement noise, where a persistent Gaussian background component coexists with impulsive or heavy-tailed disturbances. It contains two complementary experiment packages:

- `ToyModel/` — Monte Carlo simulations for Gaussian-mixture, Gaussian-plus-Laplacian-outlier, additive Gaussian–Student's t, and additive Gaussian–Slash measurement-noise settings.
- `UrbanNav_Experiment/` — A two-dimensional GNSS/UrbanNav trajectory-estimation experiment with included processed MATLAB data, filter implementations, saved outputs, and plotting scripts.

## Background and Research Significance

Robust state estimation becomes challenging when measurement errors contain two coexisting sources of uncertainty: a persistent Gaussian background component and intermittent impulsive disturbances. This situation arises in practical sensing systems such as urban GNSS navigation, where ordinary measurement uncertainty remains present while multipath, blockage, or other abnormal events introduce occasional large errors.

A conventional Kalman filter models the total measurement error by a single Gaussian distribution. Many robust filters instead replace the total error by one heavy-tailed distribution. Although such models can mitigate outliers, they do not explicitly preserve the physical additive structure in which Gaussian background uncertainty remains active during impulsive events.

The proposed additive Gaussian--Gaussian scale mixture (AGGSM) framework addresses this issue by modeling the measurement noise as

$$
\mathbf{v}_k
============

\mathbf{v}^{\mathrm{B}}_k
+
\mathbf{v}^{\mathrm{I}}_k,
$$

where (\mathbf{v}^{\mathrm{B}}_k) represents the persistent Gaussian background component and (\mathbf{v}^{\mathrm{I}}_k) represents the impulsive GSM component. This formulation preserves both uncertainty sources within one recursive filtering model rather than approximating their sum by a single distribution.

A key feature of AGGSM is the variance offset (d>0), which prevents the latent impulsive scale from approaching zero and retains a nonzero Gaussian background contribution during abnormal measurements. Combined with a deterministic MAP-based latent-variable update, the resulting robust Kalman filters provide a tractable recursive solution for additive hybrid-noise estimation.

The repository includes Gaussian-mixture, Gaussian-plus-Laplacian, additive Gaussian--Student's (t), additive Gaussian--Slash, and UrbanNav experiments to evaluate this modeling principle under both controlled and real-world navigation conditions.


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
