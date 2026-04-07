# Shoulder Modelling (Python Workflows)

This repository contains Python-based optimization workflows for shoulder biomechanics, focused on:

- Scapulohumeral rhythm (SHR) prediction.
- Muscle parameter calibration.
- Validation simulations on daily motions.
- Euler vs quaternion kinematic tracking.

The core model equations and objective terms are implemented in `Python/equations.py`, and trajectory/data utilities are implemented in `Python/trajectory_lib.py`.

## Scope Of This README

This is a first-draft README focused only on the Python code and notebooks.
MATLAB scripts and functions are intentionally not documented here yet.

## Repository Structure (Python)

- `Python/equations.py`: symbolic/numeric model equations, muscle-tendon polynomials, objective terms, and EoM builders.
- `Python/trajectory_lib.py`: trajectory interpolation, EMG import/interpolation, quaternion/Euler helpers, and export to `.mat`/`.mot`.
- `Python/SHR_prediction.ipynb`: SHR sensitivity simulations (healthy vs RC-limited, GH stability weighting sweeps).
- `Python/muscle_params_calibration.ipynb`: optimization of muscle group scalers (`fmax`, `lceopt`) using trajectory + EMG-informed terms.
- `Python/quat_tracking_validation.ipynb`: validation on functional motions (for example shelf reaching, lifting).
- `Python/quat_tracking_elevations.ipynb`: quaternion-based tracking experiments.
- `Python/eul_tracking_elevations.ipynb`: Euler-based tracking experiments.
- `Python/gen_matlab_functions.ipynb`: symbolic-generation helper notebook (legacy/utility workflow).

## Python Environment

Recommended Python version: 3.10+

Main packages used in notebooks/modules:

- `numpy`
- `scipy`
- `sympy`
- `matplotlib`
- `opty`
- `jupyter`

Install example:

```bash
pip install numpy scipy sympy matplotlib opty jupyter
```

## Data And Inputs

The Python workflows expect participant-specific data under `Motions/` and EMG files referenced by notebook paths.

Common files used by notebooks:

- OpenSim-derived model structs (for example `OS_model.mat`, `OS_model_prediction.mat`).
- Calibrated parameter files (for example `calibrated_params.mat`).
- Motion trajectories in MATLAB structs used by interpolation utilities.
- EMG MATLAB files (for calibration and some tracking variants).

Most notebooks use relative paths such as `../Motions/<participant>/...`, so they are intended to be run from the `Python/` directory context.

## Core API In `equations.py`

### Objective Builders

- `objective_traj_quat(...)`: quaternion-based trajectory tracking objective and Jacobians.
- `objective_traj_eul(...)`: Euler-based trajectory tracking objective and Jacobians.
- `objective_min_activation(...)`: activation minimization with optional EMG matching.
- `objective_state_diff(...)`: temporal smoothness penalty on states/inputs.
- `objective_max_GH_stab(...)`: GH stability objective from normalized reaction force components.
- `objective_scalers_restriction(...)`: regularization for calibrated muscle scalers.

### Muscle And Passive Structure Terms

- `polynomials_quat(...)`: quaternion-based muscle-tendon force/torque evaluation from polynomial and analytic paths, optional parameter scaling, optional RC force limiting.
- `polynomials_euler(...)`: Euler-based analog of polynomial muscle torque generation.
- `muscle_force(...)`: Hill-type muscle model.
- `act_dynamics(...)`: first-order activation dynamics.
- `conoid_force(...)`: conoid ligament force law.

### Equations Of Motion

- `create_eoms_quat_w_RF(...)`: quaternion EoMs including GH reaction-force states.
- `create_eoms_quat_no_RF(...)`: quaternion EoMs without GH reaction-force states.
- `create_eoms_eul(...)`: Euler-angle EoMs.

### Kinematics/Math Helpers

Quaternion multiplication, inverse, rotation matrices, transform helpers, analytic muscle lengths, and joint spring terms are included to support symbolic construction and numeric evaluation.

## Utilities In `trajectory_lib.py`

- `exp_trajectory_quat(...)`, `exp_trajectory_eul(...)`: load/interpolate trajectories onto optimization nodes.
- `exp_trajectory_quat_myobj(...)`, `exp_trajectory_eul_myobj(...)`: convert raw trajectory states to objective-tracking representations.
- `exp_emg(...)`: interpolate EMG data to optimization nodes.
- `initial_guess_from_solution(...)`: warm-start from previous `.mat` solution.
- `sol2struct(...)`: export optimization results to MATLAB struct format.
- `sol2mot_quat(...)`, `sol2mot_eul(...)`: export `.mot` files for OpenSim-compatible playback.

## Workflow Summary

### 1) Build Dynamics + Muscle Terms

Typical pattern in notebooks:

1. Load participant model struct from `Motions/...`.
2. Build EoMs (`create_eoms_*`).
3. Build muscle and passive torques (`polynomials_*`).
4. Optionally append activation dynamics states/equations.

### 2) Build Objectives

Objectives are composed as weighted sums of:

- Trajectory tracking (segment orientation/position terms).
- Activation and optionally excitation effort.
- GH stability penalty.
- Smoothness penalties on speeds, reaction-force states, and/or excitations.
- Optional calibration parameter regularization.

### 3) Solve NLP

Notebooks use `opty.Problem(...)` with:

- Midpoint transcription.
- Bounds on states/inputs.
- Quaternion unit-norm constraints at selected nodes.
- Problem-specific initial guesses (experimental trajectory, previous solutions, or calibration outputs).

### 4) Save Results

Solutions are written to:

- `.mat` result structs via `sol2struct(...)`.
- `.mot` kinematic files via `sol2mot_quat(...)` or `sol2mot_eul(...)`.

## Notebook Roles

- `SHR_prediction.ipynb`:
  - Runs multiple simulations with varying GH stability weights and rotator cuff limitation levels.
  - Includes healthy, RC-limited, and non-calibrated parameter scenarios.
  - Produces `res_SHR_*` outputs.

- `muscle_params_calibration.ipynb`:
  - Calibrates grouped muscle scalers (`fmax`, `lceopt`) for selected muscle groups.
  - Uses trajectory tracking, GH stability, smoothness, and EMG-informed objective components.
  - Intended to produce calibrated parameter files used by downstream simulations.

- `quat_tracking_validation.ipynb`:
  - Validates quaternion model on functional tasks (for example shelf reaching/lifting).
  - Can compare calibrated vs non-calibrated parameter settings.

- `quat_tracking_elevations.ipynb` and `eul_tracking_elevations.ipynb`:
  - Compare formulation and behavior between quaternion and Euler tracking workflows for elevation-type motions.

## Current Notes

- Some notebooks include legacy cells/function names from older iterations; use current function signatures from `Python/equations.py` and `Python/trajectory_lib.py` as source of truth.
- File paths are participant/motion specific. If a notebook fails with missing file errors, verify the selected participant, motion folder names, and expected `.mat` inputs.

## Next Documentation Steps

Potential follow-ups after this first draft:

1. Add a reproducible "quick start" run for one participant/motion with exact expected outputs.
2. Add per-notebook input/output tables (required files, generated files, key parameters).
3. Add troubleshooting notes for common optimization failures (initial guess, bounds, infeasibility).
