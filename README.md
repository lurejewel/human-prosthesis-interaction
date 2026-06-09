# Predictive Forward Simulation of Human Walking
*A Reflex-Based Musculoskeletal Model Optimised via CMA-ES*

## Overview

This repository provides a predictive, experimental-data-free forward dynamic
simulation framework for human walking.  A lower-limb musculoskeletal model is
driven by a muscle-reflex controller whose parameters are tuned by the Covariance
Matrix Adaptation Evolution Strategy (CMA-ES).  No motion capture data or
predefined joint trajectories are required.

**New in v2.0** — the reflex controller can now be either:

- a **hand-crafted** physiological controller with 28 optimisable parameters
  (legacy, v1.x), or
- a **LASSO-informed sparse linear phase controller** where each gait phase
  uses a learned 22 × 7 coefficient matrix (v2.0).  The LASSO controller is
  initialised from experimental data and then refined by CMA-ES inside the
  forward-dynamics loop.

> *The hand-crafted controller is documented in:*
> Jin W, Liu J, Zhang Q, et al. *Forward dynamics simulation of a simplified
> neuromuscular-skeletal-exoskeletal model based on the CMA-ES optimization
> algorithm: framework and case studies.* Multibody System Dynamics, 2024,
> 62(4): 525–558.

---

## Key Components

| Module | Description |
|--------|-------------|
| **Musculoskeletal Model** | 16 DOF, 14 Hill-type muscles (lower limbs only). Forward-dynamic simulation via OpenSim. |
| **Gait Phase Detection** | Real-time per-leg phase identification (5 phases: Early Stance, Late Stance, Liftoff, Swing, Landing). |
| **Muscle Reflex Controller** | Two interchangeable backends: (a) hand-crafted CNS-inspired reflex laws, (b) LASSO sparse linear per-phase model. |
| **Optimisation (IPOP-CMA-ES)** | Restart strategy with sigma boosts; population size doubles on restart. Parallelised over 6 workers. |
| **Fitness Evaluation** | Multi-objective: gait completeness, velocity tracking, joint hyperextension, GRF ceiling, metabolic effort. |

---

## Code Structure

```
.
├── demo_predictiveForwardSimulation_humanModel.m   % Main entry script
├── test_single_parameter.m                         % Single-parameter test + visuals
├── assets/                                         % Marker sets, IK/ID/RRA data, muscle CSV
├── model/                                          % OpenSim .osim model + geometry
├── results/                                        % Optimisation output (.mat) and simulation logs
├── functions/
│   ├── core/          % Simulation pipeline, CMA-ES, controllers, fitness
│   └── utils/         % I/O helpers (sto/trc/mat), vector conversion
├── tests/              % Unit tests (pure MATLAB, no OpenSim required)
├── docs/               % Format specifications
├── tools/              % OpenSim pipeline scripts (scale, IK, ID, RRA, SO)
├── Sparse Group LASSO Validation/  % LASSO fitting experiments + .mat export
└── README.md
```

---

## Requirements

| Software | Version / Notes |
|----------|-----------------|
| **MATLAB** | R2024b recommended |
| Parallel Computing Toolbox | required for `parfor` |
| **OpenSim** | 4.0 or later, with the MATLAB API configured |

> **OpenSim-MATLAB Interface**
> Instructions for configuring the OpenSim-MATLAB interface can be found at:
> https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53089380/Scripting+with+Matlab
> After successful configuration, OpenSim classes (e.g. `org.opensim.modeling.*`)
> should be accessible directly from MATLAB.

---

## How to Run

### 1. Full optimisation (IPOP-CMA-ES)

```matlab
demo_predictiveForwardSimulation_humanModel
```

The script loads a controller definition, initialises CMA-ES, and iterates until
stall criteria are met.  Results are automatically saved to `results/opt_result_*.mat`.

### 2. Quick single-parameter test

```matlab
test_single_parameter
```

Loads a previously saved `result.bestPara`, runs **one** forward simulation,
prints the fitness breakdown, and displays a 3×3 diagnostic figure (joint
angles, pelvis trajectory, gait phases, muscle excitations, GRF, speed).
Optionally plays a 3D visualisation of the gait.

---

## Controller Modes

The controller type is detected automatically by `ModelInfo.read_muscleReflex_array`:

| Mode | How to activate | Optimised parameter count |
|------|----------------|---------------------------|
| **Legacy hand-crafted** | Omit `reflexParamMap` / `reflexTemplate` on `modelInfo` | 28 (fixed) |
| **LASSO sparse linear** | Set `reflexParamMap` and `reflexTemplate` (loaded from `results/lasso_controller_result.mat`) | Variable — sum over phases of `nnz(mask{p}) + 7` |

The LASSO path is active in the current `demo_predictiveForwardSimulation_humanModel.m`.
To revert to the legacy controller, comment out the LASSO block and uncomment the
old `initPara = load(...)` line.

### LASSO controller format

See [`docs/LASSO_CONTROLLER_FORMAT.md`](docs/LASSO_CONTROLLER_FORMAT.md) for the
`.mat` file specification.  A LASSO `.mat` can be generated from experimental
data using the scripts in `Sparse Group LASSO Validation/`.

---

## Simulation and Optimisation Settings (configurable in the demo)

| Parameter | Default | Notes |
|-----------|---------|-------|
| Target walking speed | 1.0 m/s | |
| Simulation duration | 10 s | |
| Integration time step | 0.005 s | |
| `sigma` (CMA-ES initial step) | 0.02 | |
| `softPatience` | 30 gens | stall before sigma boost (×2.0) |
| `patience` | 100 gens | stall before IPOP restart |
| `maxRestarts` | 3 | IPOP restarts before hard stop |
| `gMax` | 1000 | absolute generation cap |
| Parallel workers | 6 | adjust to your CPU core count |

The optimisation terminates when:
- No fitness improvement occurs for `patience` generations **and** all
  `maxRestarts` IPOP restarts have been exhausted, or
- `gMax` generations are reached.

---

## Output Files

| Directory | File pattern | Contents |
|-----------|-------------|----------|
| `results/` | `opt_result_yyyy-mm-dd_HH-MM-SS.mat` | `result` struct: `bestFit`, `bestPara`, `fitHistory`, `timeHistory`, controller metadata, `bestReflexParams` |
| `results/` | `lasso_controller_result.mat` | LASSO controller definition (input to the demo) |
| `results/` | `sim_result_*.sto` | State histories (optional, when `simConfig.saveSTO = 1`) |
| `results/` | `opensim*.log` | OpenSim runtime logs |

---

## Scope and Limitations

- Focused on **level-ground, steady-state walking** at 1.0 m/s
- Reflex parameters are **bilaterally symmetric** (left and right legs share
  the same parameters)
- Upper-body dynamics are not included
- The LASSO controller requires a pre-fitted `.mat` file; it does **not**
  perform LASSO regression inside the optimisation loop

---

## Planned Extensions

- **Independent reflex parameters for left and right muscles** *(next minor release)*
- **Extended musculoskeletal models** with additional muscle actuators
- **Human-prosthesis coupled models** *(next major release)*
- **Additional motor tasks** beyond walking (slope walking, running, stair
  climbing; timeline TBD)

---

## Citation

If you use this code in your research, please cite the corresponding publication:

> Jin W, Liu J, Zhang Q, et al. *Forward dynamics simulation of a simplified
> neuromuscular-skeletal-exoskeletal model based on the CMA-ES optimization
> algorithm: framework and case studies.* Multibody System Dynamics, 2024,
> 62(4): 525–558.

---

## Contact

Wei JIN — Peking University — wjin24@stu.pku.edu.cn