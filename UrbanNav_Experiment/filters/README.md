# UrbanNav filter implementations

This directory contains the filter functions called by `run_UrbanNav_full_comparison.m` and the UrbanNav tuning scripts. Each function takes aligned position measurements, ground truth, control inputs, sampling interval, and method-specific parameters, and returns time-indexed error metrics.

## Methods

| File | Method |
|---|---|
| `KF_UrbanNav.m` | Conventional Kalman filter baseline. |
| `PF_UrbanNav.m` | Particle-filter baseline. |
| `Huber_UrbanNav.m` | Huber robust-filter baseline. |
| `RKF_Gaussian_UrbanNav.m` | Gaussian robust Kalman-filter baseline. |
| `RKF_ST_UrbanNav.m` | Student's t robust Kalman filter. |
| `RKF_GSTM_UrbanNav.m` | Gaussian–Student's t mixture robust Kalman filter. |
| `RKF_TAU_UrbanNav.m` | Divergence-based RKF-TAU baseline. |
| `RKF_AGSMG_UrbanNav.m` | Additive Gaussian–GSM / AGST implementation; the legacy filename is retained for compatibility. |
| `RKF_Slash_UrbanNav.m` | Additive Gaussian–Slash implementation. |

The parent scripts add this directory with:

```matlab
addpath('filters')
```

Run the parent experiment from `UrbanNav_Experiment/` so this relative path resolves correctly.
