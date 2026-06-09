# LASSO Reflex Controller — `.mat` Result Format

This document describes the layout expected by `load_lasso_reflex_controller.m`
when loading a LASSO-informed sparse linear-phase controller.

---

## File name convention

Any name is fine (e.g., `lasso_controller_result.mat`).  
The file must be placed under `results/` and the path must match the variable
`lassoFile` set in `demo_predictiveForwardSimulation_humanModel.m`.

---

## Variable name

The `.mat` file must contain **exactly one** of the following:

| Case | What the loader does |
|------|----------------------|
| A variable **named `lasso`** | Uses it directly. |
| A **single struct variable** with another name | If the file contains exactly one variable and it is a struct, the loader uses it as the `lasso` struct. |
| Anything else | Throws a clear error. |

---

## Required fields

### `lasso.beta` — coefficient matrices

Cell array, size **`nPhases × 1`** (or `1 × nPhases`).

Each cell `lasso.beta{p}` must be a **22 × 7** `double` matrix.

| Dimension | Meaning |
|-----------|---------|
| 22 rows   | Features (see §Feature ordering below) |
| 7 columns | Canonical same-side muscles (see §Muscle ordering below) |

Entry `beta{p}(i, j)` is the linear coefficient from feature `i` to
canonical muscle `j` in phase `p`.

> Coefficients may be positive, negative, or zero.  
> Structural zeros (indicated by `lasso.mask`) must be exactly `0`.

### `lasso.bias` — bias vectors

Cell array, size **`nPhases × 1`** (or `1 × nPhases`).

Each cell `lasso.bias{p}` must have **7 elements** (convertible to `1 × 7`).

> **Alias**: if `lasso.bias` does not exist, the loader also checks for
> `lasso.b` and uses it as the bias field.

---

## Optional fields

| Field | Type | Default if missing |
|-------|------|---------------------|
| `lasso.mask` | Cell array of `22 × 7` `logical` | Inferred as `abs(beta{p}) > 0` |
| `lasso.nPhases` | Scalar integer | `numel(lasso.beta)` |
| `lasso.phaseIds` | Vector (length `nPhases`) | `0 : (nPhases-1)` |
| `lasso.featureNames` | Cell array `22 × 1` or `1 × 22` of `char` | *(none)* |
| `lasso.muscleNames` | Cell array `7 × 1` or `1 × 7` of `char` | *(none)* |

### `lasso.mask`

A `logical` matrix of the same size as each `beta{p}`.

- `mask{p}(i,j) == true`  → feature `i` → muscle `j` is a **free parameter**
  to be optimised by CMA-ES.
- `mask{p}(i,j) == false` → the coefficient is a **structural zero**; it will
  never be non-zero and never enter the CMA-ES vector.

> **Validation rule**: if `mask` is provided, every `beta{p}` entry where
> `mask{p}` is `false` must be **exactly zero**. The loader throws an
> error if this is violated.

### `lasso.phaseIds`

Gait-phase ID numbers.  Must be in the same order as the `beta`/`bias`/`mask`
cell arrays.

For the current 5-phase gait model, expected values are:

```
[0, 1, 2, 3, 4]
```

| Phase ID | Name |
|----------|------|
| 0 | Early Stance |
| 1 | Late Stance |
| 2 | Liftoff |
| 3 | Swing |
| 4 | Landing |

---

## Feature ordering — 22 elements (rows of `beta`)

| Index | Feature | Source |
|-------|---------|--------|
| 1–7 | Normalised fibre length of same-side muscles 1–7 | `modelInfo.dy.muscle.lCEN` |
| 8–14 | Normalised MTU force of same-side muscles 1–7 | `modelInfo.dy.muscle.fATN` |
| 15 | Pelvis tilt angle | `stateHistory('pelvis_tilt/value')` |
| 16 | Pelvis tilt angular velocity | `stateHistory('pelvis_tilt/speed')` |
| 17 | Same-side hip angle | `stateHistory('hip_flexion_{r/l}/value')` |
| 18 | Same-side hip angular velocity | `stateHistory('hip_flexion_{r/l}/speed')` |
| 19 | Same-side knee angle | `stateHistory('knee_flexion_{r/l}/value')` |
| 20 | Same-side knee angular velocity | `stateHistory('knee_flexion_{r/l}/speed')` |
| 21 | Same-side ankle angle | `stateHistory('ankle_dorsiflexion_{r/l}/value')` |
| 22 | Same-side ankle angular velocity | `stateHistory('ankle_dorsiflexion_{r/l}/speed')` |

> All joint kinematics use the same state-map indexing and sign conventions
> as the original hand-crafted controller.  Features 15–22 are read from
> `modelInfo.st.model.initPose` for `frameIndex == 1`, otherwise from
> `modelInfo.dy.stateHistory` at `frameIndex - 1`.

---

## Canonical muscle ordering — 7 per side (columns of `beta`)

The 7 same-side muscles in the canonical (right-leg) order, matching
`assets/coupled_human-prosthesis_model_muscleProp.CSV` rows 2–8:

| Column | Muscle name (canonical, without suffix) | Full name (OpenSim) |
|--------|----------------------------------------|---------------------|
| 1 | `hamstrings` | `hamstrings_r` / `hamstrings_l` |
| 2 | `glut_max` | `glut_max_r` / `glu_max_l` |
| 3 | `iliopsoas` | `iliopsoas_r` / `iliopsoas_l` |
| 4 | `vasti` | `vasti_r` / `vasti_l` |
| 5 | `gastroc` | `gastroc_r` / `gastroc_l` |
| 6 | `soleus` | `soleus_r` / `soleus_l` |
| 7 | `tibia` | `tibia_r` / `tibia_l` |

> **Bilateral symmetry**: the same `beta{p}` and `bias{p}` are used for
> both the right leg and the left leg.  The only difference between right
> and left excitation comes from their own phase and their own same-side
> feature vector.

---

## Complete MATLAB example

This snippet generates a valid `lasso_controller_result.mat` for the
5-phase, 22-feature, 7-muscle controller.  Customise `beta` and `bias` with
your own LASSO-fitted values.

```matlab
nPhases = 5;
nFeatures = 22;
nMuscles = 7;

lasso = struct();
lasso.nPhases = nPhases;
lasso.phaseIds = 0:(nPhases-1);          % [0 1 2 3 4]
lasso.beta = cell(nPhases, 1);
lasso.bias = cell(nPhases, 1);
lasso.mask = cell(nPhases, 1);

% Optional labels (display only, not used in computation)
lasso.featureNames = {
    'len_ham'; 'len_glu'; 'len_ili'; 'len_vas'; 'len_gas'; 'len_sol'; 'len_tib';
    'frc_ham'; 'frc_glu'; 'frc_ili'; 'frc_vas'; 'frc_gas'; 'frc_sol'; 'frc_tib';
    'pelvis_tilt'; 'pelvis_tilt_vel';
    'hip_ang'; 'hip_vel'; 'knee_ang'; 'knee_vel'; 'ankle_ang'; 'ankle_vel'
    };
lasso.muscleNames = {'hamstrings'; 'glut_max'; 'iliopsoas'; 'vasti'; ...
                      'gastroc'; 'soleus'; 'tibia'};

for p = 1:nPhases
    % Replace with your LASSO-fitted values
    lasso.beta{p} = your_fitted_beta_for_phase_p;   % 22 × 7 double
    lasso.bias{p} = your_fitted_bias_for_phase_p;   % 1 × 7 double

    % Either provide mask explicitly ...
    lasso.mask{p} = your_mask_for_phase_p;           % 22 × 7 logical

    % ... or let the loader infer it from abs(beta) > 0 (omit mask field).
end

save('results/lasso_controller_result.mat', 'lasso');
```

---

## Quick checklist before saving

- [ ] File saved under `results/` as a `.mat`.
- [ ] Contains either a variable named `lasso`, or exactly one struct variable.
- [ ] `lasso.beta` is `nPhases × 1` cell, each cell `22 × 7`.
- [ ] `lasso.bias` (or `lasso.b`) is `nPhases × 1` cell, each cell has 7 values.
- [ ] Every `beta{p}` entry outside `mask{p}` is exactly `0`.
- [ ] `lasso.phaseIds` matches the phase IDs output by `cal_gait_phase.m` (`0:4`).
- [ ] Feature rows 1–14 follow the canonical muscle order (ham, glu, ili, vas, gas, sol, tib).
