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
| **Musculoskeletal Model** | 9 DOF (lower limbs + pelvis), 14 Hill-type muscles (7 per leg). Forward-dynamic simulation via OpenSim. |
| **Static Optimisation** | Iterative QP solver determines physiologically accurate initial muscle activations from first-frame kinematics and joint moments, replacing the uniform 0.05 default. |
| **Gait Phase Detection** | Real-time per-leg phase identification (5 phases: Early Stance, Late Stance, Liftoff, Swing, Landing). |
| **Muscle Reflex Controller** | Two interchangeable backends: (a) hand-crafted CNS-inspired reflex laws (legacy), (b) LASSO sparse linear per-phase model (active). |
| **Optimisation (IPOP-CMA-ES)** | Restart strategy with sigma boosts; population size doubles on restart. Parallelised over all available cores. |
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
│   ├── core/          % Simulation pipeline, CMA-ES, controllers, SO, fitness
│   └── utils/         % I/O helpers (sto/trc/mat), checkpoint, vector conversion
├── tests/              % Unit + validation tests (static optimisation, LASSO, muscle moments)
├── docs/               % Format specifications
├── tools/              % OpenSim pipeline scripts (scale, IK, ID, RRA, SO)
├── Sparse Group LASSO Validation/  % LASSO fitting experiments + UN.sto reference data
└── README.md
```

---

## Requirements

| Software | Version / Notes |
|----------|-----------------|
| **MATLAB** | R2024b recommended |
| Parallel Computing Toolbox | required for `parfor` |
| Optimization Toolbox | required for `quadprog` (iterative static optimisation) |
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

On a fresh start, the script:
1. Reads initial kinematics and joint moments from `UN.sto` (first frame).
2. Runs **iterative static optimisation** to compute physiologically accurate
   initial muscle activations (`a_opt`) via convex QP.
3. Initialises the LASSO controller and CMA-ES.
4. Iterates until stall criteria are met, injecting `a_opt` at every particle
   reset.

Results are automatically checkpointed every 10 generations and saved to
`results/opt_result_*.mat` on completion.  Set `resumeFromCheckpoint = true`
to restart from the last saved snapshot.

### 2. Quick single-parameter test

```matlab
test_single_parameter
```

Loads a previously saved `result.bestPara`, runs **one** forward simulation
(with static-optimisation initial activations), prints the fitness breakdown,
and displays diagnostic figures (joint angles, pelvis trajectory, gait phases,
muscle excitations, GRF, joint torques, gait-cycle-normalised plots).
Optionally plays a 3D visualisation of the gait.

### 3. Static optimisation validation (standalone)

```matlab
test_static_optimization_validation
```

Compares QP-computed initial muscle activations against the UN.sto reference:
scatter plot, grouped bar chart, joint moment reconstruction, and convergence
history.  Useful for verifying the static optimisation setup independently of
the full simulation pipeline.

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
| Static optimisation (SO) max iterations | 10 | inner QP loop |
| SO convergence tolerance | 1e-4 | infinity-norm of activation change |
| SO activation bounds | [0.01, 1.0] | |
| `sigma` (CMA-ES initial step) | 0.02 | |
| `softPatience` | 30 gens | stall before sigma boost (×2.0) |
| `patience` | 100 gens | stall before IPOP restart |
| `maxRestarts` | 3 | IPOP restarts before hard stop |
| `gMax` | 1000 | absolute generation cap |
| Parallel workers | `feature('numcores')` | auto-detected |

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
- Initial kinematics and joint moments are read from `UN.sto` (a pre-computed
  SCONE reference simulation); the static optimisation depends on this file
- Reflex parameters are **bilaterally symmetric** (left and right legs share
  the same parameters)
- Upper-body dynamics are not included
- The LASSO controller requires a pre-fitted `.mat` file; it does **not**
  perform LASSO regression inside the optimisation loop

---

## Recent Changes (v2.1)

- **Iterative static optimisation**: first-frame muscle activations are now
  computed via convex QP instead of using a uniform 0.05 default.  The solver
  accounts for muscle force-length-velocity relationships, passive forces, and
  knee coordinate-limit torques.
- **UN.sto-driven initialisation**: `initPose` is read from the SCONE reference
  simulation (`UN.sto`) rather than hardcoded.
- **Refactored entry point**: fresh-start setup is encapsulated in
  `prepare_fresh_start.m`; the iterative SO solver lives in
  `iterative_static_optimization.m`.
- **Checkpoint now stores `a_opt`**: resume is fully transparent.

## Planned Extensions

- **Independent reflex parameters for left and right muscles**
- **Extended musculoskeletal models** with additional muscle actuators
- **Human-prosthesis coupled models**
- **Additional motor tasks** beyond walking (slope walking, running, stair
  climbing)

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