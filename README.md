# code_openacc

This folder contains the runnable supplementary code for the proposed
chance-constrained Gaussian truncation method. Baseline truncation methods have
been removed; the retained `CC`, `CC1`, `CC2`, `CC3`, and `CC4` labels all refer
to the proposed method with different risk thresholds. The `EKF`, `UKF`, and
`GSF` labels are filter backends used together with the proposed truncation
method, not competing truncation methods.

## Contents

- `Truncation/`: proposed chance-constrained truncation routines.
- `batch reaction/`: gas-phase reaction simulation.
- `robot localization/`: landmark-based localization simulation.
- `road tracking/`: real-world roadway tracking experiment.
- `gaussKLD.m`: shared Gaussian KL-divergence helper.
- `init_code_openacc.m`: package path initializer.

Risk-threshold labels used by the retained examples:

- `CC1`: 0.05%
- `CC2`: 0.5%
- `CC3`: 5%
- `CC4`: 50%

## Dependencies

Install the dependencies before running the examples.

1. MATLAB.
   The code uses modern MATLAB syntax and plotting utilities such as
   `exportgraphics`. MATLAB R2021a or newer is recommended.

2. MOSEK for MATLAB.
   MOSEK is not bundled with this repository. Install MOSEK separately, obtain
   a valid MOSEK license, and run the MATLAB setup step provided by MOSEK so
   that the `mosekopt` function is visible on the MATLAB path. This package does
   not hard-code a local MOSEK installation path.

   After installation, verify MOSEK in MATLAB:

   ```matlab
   [rcode, res] = mosekopt('symbcon echo(0)');
   assert(rcode == 0)
   ```

   If MATLAB cannot find `mosekopt`, reopen MATLAB after running the MOSEK setup
   script, or add the MOSEK MATLAB toolbox path according to the MOSEK
   installation guide for your operating system.

3. MATLAB toolboxes.
   The examples call the following MATLAB functions:

   - `trackingEKF`, `trackingUKF`, `trackingGSF`: Sensor Fusion and Tracking
     Toolbox.
   - `mvnrnd`, `normrnd`, `norminv`: Statistics and Machine Learning Toolbox.
   - `wrapToPi`: commonly available through Mapping Toolbox or Robotics System
     Toolbox. If it is unavailable in your MATLAB installation, replace it with
     an equivalent angle-wrapping helper.
   - `exportgraphics`, `array2table`: MATLAB plotting and table utilities
     available in recent MATLAB releases.

## Setup

Start MATLAB, move to the package root, and add the package to the MATLAB path:

```matlab
cd('path/to/code_openacc')
addpath(genpath(pwd))
```

Each entry script also calls `init_code_openacc.m`, so the scripts can be run
from the package root or from their own experiment folders.

Some examples use precomputed MOSEK problem templates:

- `road tracking/mosek_data_list.mat`
- `robot localization/mosek_data_list2.mat`

These files are included with the release and are required by the corresponding
entry scripts.

## Experiment 1: Batch Reaction Simulation

Run from MATLAB:

```matlab
cd('path/to/code_openacc/batch reaction')
total_batch_sim
```

This script runs the gas-phase reaction simulation with linear state
constraints. It prints tables for:

- RMSE
- average computation time per step
- average violation probability
- maximum violation probability

Lower RMSE and lower computation time are better. The violation-probability
tables report how often Monte Carlo samples violate the retained state
constraints.

## Experiment 2: Robot Localization Simulation

Run the main localization simulation:

```matlab
cd('path/to/code_openacc/robot localization')
total_localize_sim
```

This script evaluates landmark-based localization with field-of-view-induced
constraints. It prints average computation time, RMSE, average violation
probability, and maximum violation probability for the retained CC settings.
The `CC1` through `CC4` columns correspond to the risk thresholds listed above.

## Experiment 3: Road Tracking Experiment

Plot a representative real-world tracking case:

```matlab
cd('path/to/code_openacc/road tracking')
plot_tracking
```

This script loads the included roadway data, runs the retained CC method, and
saves:

- `tracking_scenario.pdf`: road geometry and representative tracking result.
- `tracking_rmse.pdf`: RMSE curve over time.

Run the aggregate road-tracking evaluation:

```matlab
cd('path/to/code_openacc/road tracking')
total_trackingSim
```

This script runs the included trajectories `exp_data1.mat` through
`exp_data5.mat` and prints aggregate tables for computation time, RMSE, position
RMSE, heading RMSE, velocity RMSE, average violation probability, and maximum
violation probability.

## Quick-Run Defaults

The default Monte Carlo counts are intentionally small so the supplementary
examples are easy to validate:

- `batch reaction/total_batch_sim.m`: `Nmc = 10`
- `robot localization/total_localize_sim.m`: `Nmc = 5`
- `road tracking/total_trackingSim.m`: `Nmc = 2`

Increase these values if you want lower-variance statistics for paper-quality
tables.

## Expected Outputs

Console tables summarize the main numerical results. In all experiments:

- RMSE measures estimation accuracy; lower is better.
- Computation time is reported in milliseconds per filtering step; lower is
  faster.
- Violation probability estimates how often sampled states violate the
  constraints; lower is safer and should be interpreted relative to the selected
  CC risk threshold.

Generated PDFs and `.mat` result files are analysis outputs and can be safely
deleted and regenerated. Keep the bundled `mosek_data_list.mat` and
`mosek_data_list2.mat` files because they are required runtime templates.
