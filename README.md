# Robust Kalman Filters for Additive Gaussian–GSM Distributions

MATLAB implementations of **RKF-AGST** and **RKF-AGSlash**, with synthetic
measurement-noise benchmarks and an offline UrbanNav positioning experiment.
This repository accompanies the manuscript *Robust Kalman Filters for Additive
Gaussian–GSM Distributions*. The current experiment code is in [`journal/`](journal/).

## Model and methods

The additive Gaussian–Gaussian scale mixture (AGGSM) model represents measurement
noise as the sum of a Gaussian background component and a GSM disturbance.
Conditional on the mixing scale, its covariance is `(d + eta_k) * Sigma_k`.
For fixed positive `d` and positive-definite `Sigma_k`, the offset retains the
modeled background contribution `d * Sigma_k` as the impulsive scale decreases.
It does **not** prevent `eta_k` itself from approaching zero, nor does the
parameterization alone establish separate identification of physical noise sources.

A VB–MAP inference scheme gives a Kalman-type prediction–correction recursion,
with inverse-Gamma and Slash mixing laws yielding RKF-AGST and RKF-AGSlash.
The experiments examine conditions where each method is useful; they do not
claim one filter is best under every non-Gaussian noise distribution.

## Quick start: current manuscript

Tested with MATLAB R2024b and Statistics and Machine Learning Toolbox.

```matlab
cd journal
validate_release          % short synthetic replay + full UrbanNav replay
run_synthetic(1000)       % 32 conditions, 14 methods, 1,000 trajectories each
run_urban                % Odaiba-Trimble, with common scored epochs
```

The current comparison includes KF, PF, Huber KF, RKF-Gaussian, RKF-ST,
RKF-GSTM, RKF-TAU, GSM-filter, GGM, URKF, RKF-SGalphaS, NBKF,
RKF-AGST and RKF-AGSlash. See the [experiment guide](journal/README.md)
for parameter conventions, metrics, calibration and reproduction scope.

## Repository structure

| Directory | Contents |
|---|---|
| `journal/` | Current 32-condition experiments, corrected filter versions, frozen inputs, replay scripts and parameter records |
| `journal/exploration/` | Broader exploratory result summaries, including conditions outside the manuscript selection |
| `ToyModel/` | Earlier synthetic examples and outputs, preserved for history |
| `UrbanNav_Experiment/` | Earlier processed UrbanNav experiment, preserved for history |

The root `RKF-AGGSM.zip` is also a legacy archive, not the current release.
Use GitHub **Code → Download ZIP** or clone the repository for the current code.
Do not mix legacy filter paths with `journal/`; `setup_journal` selects the
current implementation paths.

## Experiment scope

The synthetic benchmark contains eight selected conditions for each of four
noise families: additive Gaussian–Student's t, additive Gaussian–Slash,
Gaussian mixture, and Gaussian plus Laplace outliers. Every trajectory has
30 steps; the reported statistic is the mean of trajectory RMSEs.
Selection followed exploratory sweeps, whose broader summaries are retained.
No universal performance ranking or population mean-square guarantee for
infinite-variance stress tests is claimed.

The manuscript's UrbanNav result uses Odaiba–Trimble with reference-derived
velocity controls and previously calibrated parameters. It is an offline
filter comparison, not a standalone online GNSS/INS system. Missing GNSS
observations use prediction-only steps on the common scored interval.

## Data and attribution

UrbanNav inputs are derived from the
[UrbanNav Dataset](https://github.com/IPNL-POLYU/UrbanNavDataset).
We thank its authors and maintainers. Consult the original dataset's
documentation and data-use terms. This repository includes processed replay
inputs, not the raw Tokyo sensor archive or satellite imagery.
Baseline references and adaptation conventions are listed in the
[experiment guide](journal/README.md#baseline-references); these implementations
are not presented as the baseline authors' official software.

## Acknowledgement

I would like to express my sincere gratitude to my supervisor, Prof. Xiao-Ping Zhang, for his invaluable research guidance and insightful revision feedback on this manuscript. I would also like to thank Prof. Zhenyu Liu for his valuable discussions, constructive suggestions, and continued support throughout this work.

This work originated from a course project on Bayesian Learning, [RKF-AGST](https://github.com/jgr2021/RKF-AGST), and was subsequently extended into the present study. I also thank Prof. Ercan E. Kuruoğlu and Dr. Pengcheng Hao for their valuable guidance and assistance during and after the course project.

## Citation

Please cite the associated manuscript when using this code:
Gurui Jin, Pengcheng Hao, Xiao-Ping Zhang, Zhenyu Liu and Ercan E. Kuruoğlu,
“Robust Kalman Filters for Additive Gaussian–GSM Distributions.”
The final journal bibliographic record and DOI will be added after publication.
