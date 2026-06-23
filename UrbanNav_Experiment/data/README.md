# UrbanNav input data

This directory stores the processed MATLAB data files used by the UrbanNav experiment.

## Files used directly by the main loader

`load_UrbanNav_data.m` reads and aligns the following files:

| File | Expected variable | Role in the experiment |
|---|---|---|
| `pos.mat` | `pos` | Ground-truth planar position; the first two columns are used when more are present. |
| `traj.mat` | `traj` | GNSS/measurement trajectory; the first two columns are used when more are present. |
| `data_vel.mat` | `data_vel` | Velocity input used to construct the control term. |

The loader sets `dt = 0.1` seconds and truncates the arrays to a common length before filtering.

## Additional files

- `GPS_Data.mat`
- `Ground_Truth.mat`

These files are retained as auxiliary/original-format inputs in the supplied archive. The current main loader does not directly read them.

## Reuse and distribution

These are processed experiment inputs packaged with the code. Verify the provenance, attribution requirements, and redistribution permissions of the original data before making this directory publicly accessible.
