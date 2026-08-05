# Predictive Forward Simulation of Human Walking
*A Reflex-Based Musculoskeletal Model Optimised via CMA-ES*

## Overview

This repository provides a predictive, experimental-data-free forward dynamic
simulation framework for human walking.  A lower-limb musculoskeletal model is
driven by a muscle-reflex controller whose parameters are tuned by the Covariance
Matrix Adaptation Evolution Strategy (CMA-ES).  No motion capture data or
predefined joint trajectories are required.

**Active controller (v2.1)** — the reflex controller is a **LASSO-informed
sparse linear phase controller** where each gait phase (or phase group) uses a
learned 26 × 9 coefficient matrix.  The LASSO controller is initialised from
experimental data and then refined by CMA-ES inside the forward-dynamics loop.

> **Legacy (v1.0):** a hand-crafted physiological controller with 28 optimisable
> parameters was used in earlier versions.  It has been **fully deprecated** —
> the parameter definitions file (`DEPRECATED_muscle_reflex_param_defs.m`) and
> the legacy excitation code in `cal_muscle_excitation.m` are retained for
> reference only and are no longer functional.

---

## Key Components

| Module | Description |
|--------|-------------|
| **Musculoskeletal Model** | 9 DOF (lower limbs + pelvis), 18 Hill-type muscles (9 per leg). Forward-dynamic simulation via OpenSim. Supported models: `human0918` (18 muscles, active). |
| **Static Optimisation** | Iterative QP solver determines physiologically accurate initial muscle activations from first-frame kinematics and joint moments, replacing the uniform 0.05 default. |
| **Gait Phase Detection** | Real-time per-leg phase identification (5 phases: Early Stance, Late Stance, Liftoff, Swing, Landing). |
| **Muscle Reflex Controller** | LASSO sparse linear per-phase model with optional **grouped format** (v2): phases sharing the same controller reference a single `beta`/`bias`/`mask`, eliminating parameter duplication. Legacy hand-crafted controller is deprecated. |
| **Optimisation (IPOP-CMA-ES)** | Restart strategy with soft sigma-boost before restart; population size doubles on restart. Parallelised over all available cores. Reproducible via `rngSeed` + `threefry`. |
| **Fitness Evaluation** | Multi-objective: gait completeness, velocity tracking, joint hyperextension (knee + ankle), GRF ceiling, metabolic effort. |
| **Checkpoint / Resume** | Full-state checkpoint saved every 10 generations. Resume from crash or intentional stop by setting `resumeFromCheckpoint = true`. |

---

## Code Structure

```
.
├── demo_predictiveForwardSimulation_humanModel.m   % Main entry script
├── assets/                                         % Marker sets, IK/ID/RRA data, muscle CSV
├── model/                                          % OpenSim .osim model + geometry
│   ├── human0714.osim                              %   14-muscle model (legacy)
│   ├── human0918.osim                              %   18-muscle model (active)
│   └── coupled_human-prosthesis_model.osim         %   prosthesis-coupled variant
├── results/                                        % Optimisation output (.mat), simulation logs, checkpoint
├── functions/
│   ├── core/          % Simulation pipeline, CMA-ES, controllers, SO, fitness
│   └── utils/         % I/O helpers (sto/trc/mat), checkpoint, vector conversion
├── tests/              % Unit + validation tests (SO, LASSO, muscle moments, inverse dynamics)
├── docs/               % Format specifications
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
to restart from the last saved snapshot (`results/checkpoint.mat`).
The checkpoint stores the full CMA-ES internal state, loop bookkeeping,
and all immutable configuration — resume is fully transparent.

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

### 4. Inverse dynamics & muscle moment validation

```matlab
test_inverse_dynamics_validation
test_muscle_moment_summation
```

Validates the OpenSim inverse dynamics pipeline and confirms that summed
individual muscle moments match net joint moments from the reference simulation.

### 5. LASSO controller round-trip test

```matlab
test_lasso_reflex_controller
```

Pure-MATLAB (no OpenSim) test of the LASSO controller loader/unpacker round-trip
for both legacy per-phase and grouped controller formats.

---

## Controller Modes

The controller type is detected automatically by `ModelInfo.read_muscleReflex_array`:

| Mode | How to activate | Optimised parameter count |
|------|----------------|---------------------------|
| **LASSO sparse linear** (active) | Set `reflexParamMap` and `reflexTemplate` (loaded from `results/lasso_controller_result.mat`) | Variable — sum over groups of `nnz(mask{g}) + 9` |
| **Legacy hand-crafted** (deprecated) | No longer functional; code retained for reference only | 28 (fixed) |

The LASSO path is active in the current `demo_predictiveForwardSimulation_humanModel.m`.

### LASSO controller: grouped vs. legacy format

The loader supports two LASSO struct formats:

| Format | Field | Description |
|--------|-------|-------------|
| **Grouped (v2)** ✅ recommended | `lasso.groups(g).beta/bias/mask/phases/label` | Controllers defined per **group**; phases sharing the same controller reference one parameter set — no duplication |
| **Legacy (v1)** | `lasso.beta{p}`, `lasso.bias{p}`, `lasso.mask{p}` | Per-phase cell arrays (`nPhases × 1`) |

In grouped mode, `initPara` contains only `nGroups` sets of parameters.
The unpacker expands them to per-phase `beta`/`bias` cells at runtime.
See [`docs/LASSO_CONTROLLER_FORMAT.md`](docs/LASSO_CONTROLLER_FORMAT.md) for the
complete `.mat` file specification.

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
| `rngSeed` | 2026 | `threefry` generator for reproducibility |
| `softPatience` | 30 gens | stall before soft sigma boost (×2.0) |
| `patience` | 100 gens | stall before IPOP restart |
| `maxRestarts` | 3 | IPOP restarts before hard stop |
| `minImprovement` | 1e-12 | minimum fitness delta to count as improvement |
| `gMax` | 2000 | absolute generation cap |
| Parallel workers | `feature('numcores')` | auto-detected |

The optimisation terminates when:
- No fitness improvement occurs for `softPatience` generations → a **soft
  sigma boost** (×`sigmaBoostFactor`) is applied once.
- If stalling continues for `patience` generations → an **IPOP restart** is
  triggered (population doubled, mean reset to best-so-far).
- After `maxRestarts` IPOP restarts with no further improvement, or when
  `gMax` generations are reached, optimisation stops.

---

## Output Files

| Directory | File pattern | Contents |
|-----------|-------------|----------|
| `results/` | `opt_result_yyyy-mm-dd_HH-MM-SS.mat` | `result` struct: `bestFit`, `bestPara`, `fitHistory`, `timeHistory`, controller metadata, `bestReflexParams` |
| `results/` | `checkpoint.mat` | Full CMA-ES state + loop bookkeeping for crash recovery (overwritten every 10 gens) |
| `results/` | `lasso_controller_result.mat` | LASSO controller definition (input to the demo) |
| `results/` | `sim_result_*.sto` | State histories (optional, when `simConfig.saveSTO = 1`) |
| `results/` | `opensim*.log` | OpenSim runtime logs |

---

## Scope and Limitations

- Focused on **level-ground, steady-state walking** at 1.0 m/s
- Initial kinematics and joint moments are read from `UN.sto` (a pre-computed
  SCONE reference simulation); the static optimisation depends on this file
- Reflex parameters are **bilaterally symmetric** (left and right legs share
  the same controller parameters)
- Upper-body dynamics are not included
- The LASSO controller requires a pre-fitted `.mat` file; it does **not**
  perform LASSO regression inside the optimisation loop
- The legacy hand-crafted reflex controller (v1.x) is **fully deprecated**
  and no longer functional
- The active model is `human0918` (18 muscles, 9 per leg); `human0714`
  (14 muscles, 7 per leg) is retained for legacy compatibility

---

## Recent Changes (v2.2)

- **Model upgrade**: switched from `human0714` (14 muscles, 7/leg) to
  `human0918` (18 muscles, 9/leg).  Muscle names were updated accordingly
  (hamstrings, bifemsh, glut_max, iliopsoas, rect_fem, vasti, gastroc,
  soleus, tib_ant).
- **Feature dimension update**: feature vector expanded from 22 to 26
  elements (9 fibre lengths + 9 MTU forces + 8 joint kinematics) to match
  the 9-muscle model.  Coefficient matrix is now 26 × 9.
- **Grouped LASSO format (v2)**: controllers can now be shared across phase
  groups (e.g., one stance controller for phases 0–1, one swing controller
  for phases 2–4).  Eliminates parameter duplication and reduces the CMA-ES
  search dimension.
- **Checkpoint / resume**: full CMA-ES state saved every 10 generations to
  `results/checkpoint.mat`.  Set `resumeFromCheckpoint = true` to resume
  from the last snapshot.  The checkpoint stores `a_opt` and `dofNames` for
  fully transparent recovery.
- **RNG reproducibility**: deterministic `threefry` generator with fixed
  `rngSeed = 2026`.
- **Staged stall handling**: soft sigma-boost (×2.0) is attempted before
  triggering an IPOP restart, giving the optimiser a gentler nudge.
- **Legacy controller deprecated**: the hand-crafted 28-parameter reflex
  controller is no longer functional.  The code is retained (commented out)
  in `cal_muscle_excitation.m` and `DEPRECATED_muscle_reflex_param_defs.m`.
- **New validation tests**: `test_inverse_dynamics_validation.m` and
  `test_muscle_moment_summation.m` for verifying the OpenSim pipeline.
- **Ankle limit measure**: added to the multi-objective fitness to penalise
  ankle hyperextension/hyperflexion.
- **Simulation cache**: `forward_simulation.m` pre-builds a cache of Java
  handles to avoid per-frame JNI overhead.

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