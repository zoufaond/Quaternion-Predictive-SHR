# Quaternion-based optimal control for scapulohumeral coordination

This repository contains Python-based optimization workflows for shoulder optimal control, focused on:

- Scapulohumeral rhythm (SHR) prediction.
- Muscle parameter calibration.
- Validation simulations on daily motions.
- Euler vs quaternion kinematic tracking.

For understanding the optimization problem setup used in the paper context, start with `Python/healthy_SHR_prediction_example.ipynb`.

## Repository Structure

### Python
- `Python/equations.py`: 
    - symbolic/numeric model equations, muscle-tendon polynomials, objective terms, and EoM builders.

- `Python/trajectory_lib.py`: 
    - trajectory interpolation, EMG import/interpolation, quaternion/Euler helpers, and export to `.mat`/`.mot`.

- `Python/healthy_SHR_prediction_example.ipynb`:
  - Primary notebook to understand the problem formulation and optimization setup.
  - Demonstrates healthy SHR prediction with calibrated parameters.
  - All data required to run this notebook is included in `Motions/par2/`.

- `Python/SHR_prediction.ipynb`:
  - Runs multiple simulations with varying GH stability weights and rotator cuff limitation levels.
  - Includes healthy, RC-limited, and non-calibrated parameter scenarios.
  - Produces `res_SHR_0` - `res_SHR_6` outputs.

- `Python/muscle_params_calibration.ipynb`:
  - Calibrates grouped muscle scalers (`fmax`, `lceopt`) for selected muscle groups.
  - Uses trajectory tracking, GH stability, numerical regularization, and EMG-informed objective terms.

- `Python/quat_tracking_validation.ipynb`:
  - Validates quaternion model on functional tasks (for example shelf reaching/lifting).
  - Compares calibrated vs non-calibrated parameter.

- `Python/quat_tracking_elevations.ipynb` and `Python/eul_tracking_elevations.ipynb`:
  - Compare formulation and behavior between quaternion and Euler tracking workflows for elevation-type motions.

### MATLAB files
- `create_files.m` covers the following steps for a single participant and motion:

    1. Load inverse kinematics data from an OpenSim `.mot` file.
    2. Generate simulation trajectory .mat file to be tracked in optimal control (with clavicle axial rotation to correspond to minimal conoid ligament length.).
    3. Augment training data by adding noise across viable shoulder configurations, with viability checked against OpenSim wrapping geometry.
    4. Fit muscle-tendon path polynomials and save it in the `OS_model.mat`.

- `plot_results.m`:
    - Generates all figures used in the paper.

- `das3_polynomials.m`
    - Function to calculate the polynomial approximation of muscle-tendon paths.
    - Adapted from the DAS3 project (https://github.com/dasproject/DAS3) and modified to support quaternion-based represenation.

- `das_readosim.m`
    - Function to create `OS_model.mat` MATLAB struct.
    - Adapted from the DAS3 project (https://github.com/dasproject/DAS3).

- `Matlab_functions/*`:
    - Contains helping MATLAB functions.

## Python Environment

Recommended Python version: 3.10+

Main packages used in notebooks/modules:

- `numpy`
- `scipy`
- `sympy`
- `matplotlib`
- `jupyter`

Install core packages:

```bash
pip install numpy scipy sympy matplotlib jupyter
```

For `opty`, please use the official project documentation and installation guidance, as setup can be environment-dependent:

- https://opty.readthedocs.io/
