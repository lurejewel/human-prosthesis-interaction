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

## Format versions

The loader supports two formats, auto-detected:

| Format | Detected by | Description |
|--------|-------------|-------------|
| **Grouped (v2)** ✅ recommended | `isfield(lasso, 'groups')` | Controllers are defined per **group**; phases that share the same controller reference the same `beta`/`bias`/`mask`. No parameter duplication. |
| **Legacy (v1)** | otherwise | Per-phase `beta`/`bias`/`mask` cell arrays (`nPhases × 1`). |

---

## Grouped format (v2) — recommended

### Required fields

| Field | Type | Description |
|-------|------|-------------|
| `lasso.format` | `char` | `'grouped'` |
| `lasso.nPhases` | scalar | Total number of gait phases (e.g., `5`) |
| `lasso.phaseIds` | `1 × nPhases` | Phase ID numbers (e.g., `[0 1 2 3 4]`) |
| `lasso.nGroups` | scalar | Number of distinct controller groups (e.g., `2`) |
| `lasso.groups` | `1 × nGroups` struct array | Each element defines one controller group |

### `lasso.groups(g)` fields

| Field | Type | Description |
|-------|------|-------------|
| `.label` | `char` | Human-readable label (e.g., `'stance'`, `'swing'`) |
| `.phases` | `1 × k` double | Phase IDs that use this controller (e.g., `[0 1]`) |
| `.beta` | `22 × 7` double | Coefficient matrix (features × canonical muscles) |
| `.bias` | `1 × 7` double | Bias/intercept vector |
| `.mask` | `22 × 7` logical | Structural zero mask (`true` = optimizable entry) |

### Example grouped struct

```matlab
lasso = struct();
lasso.format = 'grouped';
lasso.nPhases = 5;
lasso.phaseIds = 0:4;
lasso.nGroups = 2;
lasso.groups = struct();

% Stance controller — used by phases 0 and 1
lasso.groups(1).label  = 'stance';
lasso.groups(1).phases = [0 1];
lasso.groups(1).beta   = betaStance;     % 22×7 double
lasso.groups(1).bias   = biasStance;     % 1×7 double
lasso.groups(1).mask   = maskStance;     % 22×7 logical

% Swing controller — used by phases 2, 3, and 4
lasso.groups(2).label  = 'swing';
lasso.groups(2).phases = [2 3 4];
lasso.groups(2).beta   = betaSwing;      % 22×7 double
lasso.groups(2).bias   = biasSwing;      % 1×7 double
lasso.groups(2).mask   = maskSwing;      % 22×7 logical

save('results/lasso_controller_result.mat', 'lasso');
```

### How the loader expands grouped format

1. **`reflexTemplate.mask`** — expanded to per-phase cells (phases in the same
   group share the same mask).
2. **`reflexTemplate.groupOfPhase`** — `1 × nPhases` mapping from phase index
   → group index.
3. **`reflexParamMap.groups(g)`** — stores the layout of each group in the
   flat parameter vector (no duplication).
4. **`initPara`** — contains only `nGroups` sets of parameters.  The unpacker
   expands them to per-phase `beta`/`bias` cells at runtime.

---

## Legacy format (v1)

### Required fields

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

## Optional fields (both formats)

| Field | Type | Default if missing |
|-------|------|---------------------|
| `lasso.nPhases` | Scalar integer | `numel(lasso.beta)` |
| `lasso.phaseIds` | Vector (length `nPhases`) | `0 : (nPhases-1)` |
| `lasso.featureNames` | Cell array `22 × 1` or `1 × 22` of `char` | *(none)* |
| `lasso.muscleNames` | Cell array `7 × 1` or `1 × 7` of `char` | *(none)* |

### `lasso.mask` (legacy format only)

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
| 19 | Same-side knee angle | `stateHistory('knee_extension_{r/l}/value')` |
| 20 | Same-side knee angular velocity | `stateHistory('knee_extension_{r/l}/speed')` |
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
| 7 | `tib_ant` | `tib_ant_r` / `tib_ant_l` |

> **Bilateral symmetry**: the same `beta{p}` and `bias{p}` are used for
> both the right leg and the left leg.  The only difference between right
> and left excitation comes from their own phase and their own same-side
> feature vector.

---

## Quick checklist before saving

- [ ] File saved under `results/` as a `.mat`.
- [ ] Contains either a variable named `lasso`, or exactly one struct variable.
- [ ] **Grouped format (recommended):** `lasso.groups(g).beta`/`.bias`/`.mask`/`.phases`/`.label`.
- [ ] **Legacy format:** `lasso.beta` is `nPhases × 1` cell, each cell `22 × 7`.
- [ ] **Legacy format:** `lasso.bias` (or `lasso.b`) is `nPhases × 1` cell, each cell has 7 values.
- [ ] Every `beta` entry outside its `mask` is exactly `0`.
- [ ] `lasso.phaseIds` matches the phase IDs output by `cal_gait_phase.m` (`0:4`).
- [ ] Feature rows 1–14 follow the canonical muscle order (ham, glu, ili, vas, gas, sol, tib).
