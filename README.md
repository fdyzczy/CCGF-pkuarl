# Chance-Constrained Gaussian Filtering Examples

This repository contains the MATLAB supplementary code for the proposed
chance-constrained Gaussian truncation method. Each experiment compares an
unconstrained filter baseline (`UN`) with the proposed constrained filter
(`CC`). The `EKF`, `UKF`, and `GSF` labels are filter backends used with these
two filtering modes, not competing truncation methods.

## Installation

1. Install MATLAB and required toolboxes.

   MATLAB R2021a or newer is recommended. The examples use the following MATLAB
   functions:

   - `trackingEKF`, `trackingUKF`, `trackingGSF`: Sensor Fusion and Tracking
     Toolbox.
   - `mvnrnd`, `normrnd`, `norminv`: Statistics and Machine Learning Toolbox.
   - `wrapToPi`: commonly available through Mapping Toolbox or Robotics System
     Toolbox.
   - `array2table`: MATLAB table utility in recent MATLAB releases.

2. Install MOSEK for MATLAB before running the examples.

   MOSEK is not bundled with this repository. Install MOSEK separately, obtain a
   valid MOSEK license, and make sure MATLAB can call `mosekopt`.

   This code uses the MOSEK Optimization Toolbox interface through `mosekopt`.
   MOSEK 10 and MOSEK 11 both support this interface, but in MOSEK 11 it is
   documented as the older/deprecated Optimization Toolbox interface. Make sure
   you install and add the `mosekopt` toolbox path, not only the newer MATLAB
   API.

   Verify MOSEK in MATLAB:

   ```matlab
   which mosekopt -all
   [rcode, ~] = mosekopt('symbcon echo(0)');
   assert(rcode == 0)
   ```

   On Windows, if MATLAB finds `mosekopt.mexw64` but reports that the MEX file
   is invalid because a module cannot be found, the MOSEK binary directory is
   usually missing from `PATH`. For example, if MOSEK is installed in
   `C:\Program Files\Mosek\11.0`, run the following commands in MATLAB:

   ```matlab
   setenv('PATH', [getenv('PATH') ';C:\Program Files\Mosek\11.0\tools\platform\win64x86\bin']);
   addpath('C:\Program Files\Mosek\11.0\toolbox\r2017a');
   ```

3. Enter the project directory in MATLAB.

   ```matlab
   cd('path/to/code_openacc')
   addpath(genpath(pwd))
   ```

   Each experiment entry script also calls `init_code_openacc.m`, so the scripts
   can be run from the package root or from their own experiment folders.

## Running Experiments

These MATLAB examples are run by entering the corresponding folder and executing
the entry script. To change the risk threshold for `CC` or the number of Monte
Carlo runs, edit `risk_threshold` or `Nmc` near the top of the entry script.

### Batch Reaction Simulation

Run the gas-phase reaction simulation:

```matlab
cd('path/to/code_openacc/batch reaction')
total_batch_sim
```

Default settings:

- `risk_threshold = 1e-2`
- `Nmc = 10`

The script prints RMSE, computation time, average violation probability, and
maximum violation probability tables for `UN` and `CC`.

### Robot Localization Simulation

Run the landmark-based localization simulation:

```matlab
cd('path/to/code_openacc/robot localization')
total_localize_sim
```

Default settings:

- `risk_threshold = 0.5`
- `Nmc = 100`
- `num_landmarks = 16`

The script prints average computation time, RMSE, average violation probability,
and maximum violation probability tables for `UN` and `CC`.

### Road Tracking Experiment

Run the aggregate road-tracking evaluation:

```matlab
cd('path/to/code_openacc/road tracking')
total_trackingSim
```

Default settings:

- `risk_threshold = 0.0005`
- `Nmc = 2`

The script runs `exp_data1.mat` through `exp_data5.mat` and prints aggregate
tables for computation time, position RMSE, heading RMSE, velocity RMSE, average
violation probability, and maximum violation probability for `UN` and `CC`.

## Experiment Data and MOSEK Templates

The road-tracking experiment uses the included real-world data files:

```text
road tracking/exp_data1.mat
road tracking/exp_data2.mat
road tracking/exp_data3.mat
road tracking/exp_data4.mat
road tracking/exp_data5.mat
road tracking/constr_data.mat
```

Some examples also use precomputed MOSEK problem templates:

```text
road tracking/mosek_data_list.mat
robot localization/mosek_data_list2.mat
```

These template files are included with the release and are required by the
corresponding entry scripts.

## Repository Structure

```text
code_openacc/
  Truncation/             Proposed chance-constrained truncation routines
  batch reaction/         Batch reaction simulation
  robot localization/     Landmark-based localization simulation
  road tracking/          Real-world road-tracking experiment
  gaussKLD.m              Gaussian KL-divergence helper
  init_code_openacc.m     Package path and MOSEK initializer
```

## Outputs

Console tables summarize the main numerical results:

- RMSE measures estimation accuracy; lower is better.
- Computation time is reported in milliseconds per filtering step; lower is
  faster.
- Violation probability estimates how often sampled states violate the
  constraints. For `CC`, it should be interpreted relative to `risk_threshold`;
  for `UN`, it shows the unconstrained filter's constraint violations.

Generated PDFs and result files can be deleted and regenerated. Keep the bundled
`mosek_data_list.mat` and `mosek_data_list2.mat` files because they are required
runtime templates.
