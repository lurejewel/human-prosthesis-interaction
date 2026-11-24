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

> **Also deprecated (renamed with a `DEPRECATED_` prefix, reference only):**
> `DEPRECATED_run_B30_outer.py` — early 30 s B30 outer-loop driver, superseded
> by `run_B30_continue.py`; it has a known shared-segment race because the
> parallel chunk solves share `experiments/moco/outputs/segments/` (see
> [`docs/30s_continuation_runbook.md`](docs/30s_continuation_runbook.md));
> `functions/utils/DEPRECATED_determine_gait_phase_from_motion.m` — uses the
> obsolete 3-argument `cal_gait_phase` API and missing `subject01` assets;
> `functions/utils/DEPRECATED_read_trc_file.m` — legacy `.trc` visualization
> with missing assets.  None of these is part of an active pipeline.
> By contrast,
> [`Adaptive_LASSO_StanceSwing.m`](<Sparse Group LASSO Validation/Adaptive_LASSO_StanceSwing.m>)
> is **not** deprecated: it remains the supported legacy `UN.sto` entry point.

---

## Key Components

| Module | Description |
|--------|-------------|
| **Musculoskeletal Model** | 9 DOF (lower limbs + pelvis), 18 Hill-type muscles (9 per leg). Forward-dynamic simulation via OpenSim. Supported models: `human0918` (18 muscles, active). |
| **Static Optimisation** | Iterative QP solver determines physiologically accurate initial muscle activations from first-frame kinematics and joint moments, replacing the uniform 0.05 default. |
| **Gait Phase Detection** | Real-time per-leg phase identification (5 phases: Early Stance, Late Stance, Liftoff, Swing, Landing). |
| **Muscle Reflex Controller** | LASSO sparse linear per-phase model with optional **grouped format** (v2): phases sharing the same controller reference a single `beta`/`bias`/`mask`, eliminating parameter duplication. Legacy hand-crafted controller is deprecated. |
| **Structure Identification (stability)** | K random seeds × run-level bootstrap → per-entry selection frequency π → consensus 0-1 mask (majority vote) → frozen-mask refit on the full data. The 0-1 map is hard-frozen for CMA-ES; coefficients are only the initial guess. See [`docs/LASSO_STABILITY_IDENTIFICATION.md`](docs/LASSO_STABILITY_IDENTIFICATION.md). |
| **Optimisation (IPOP-CMA-ES)** | Restart strategy with soft sigma-boost before restart; population size doubles on restart. Parallelised over all available cores. Reproducible via `rngSeed` + `threefry`. |
| **Fitness Evaluation** | Multi-objective: gait completeness, velocity tracking, joint hyperextension (knee + ankle), GRF ceiling, metabolic effort. Every term is non-negative by construction: the effort term (cost of transport) keeps the old expression `metabolic/(mass*distance)` for ANY positive forward displacement, and falls back to a 0.5 m distance only when the net displacement is ≤ 0, so a fall/backward drift can never become a negative cost. `measure_simResults` also returns the per-term breakdown and warns if any term goes negative. |
| **Checkpoint / Resume** | Full-state checkpoint saved every 10 generations. Resume from crash or intentional stop by setting `resumeFromCheckpoint = true`. |

---

## Code Structure

```
.
├── demo_predictiveForwardSimulation_humanModel.m   % Main entry script
├── assets/                                         % Marker sets, IK/ID data, muscle CSV
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
(with static-optimisation initial activations), prints the fitness breakdown
(completeness / velocity / knee / ankle / GRF / effort), and displays
diagnostic figures (joint angles, pelvis trajectory, gait phases, muscle
excitations, GRF, joint torques, gait-cycle-normalised plots).
Optionally plays a 3D visualisation of the gait.

Lightweight alternative (no plots, no 3D, just the fitness terms):

```matlab
tests/evaluate_opt_result_fitness
```

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

### 6. MocoInverse → LASSO (measured Vicon/GRF/sEMG data)

The old `UN.sto → LASSO` path is unchanged.  For measured-data MocoInverse
results, run the new two-stage interface:

```matlab
cd('Sparse Group LASSO Validation')
moco_inverse_to_lasso_dataset   % Moco solution + IK + GRF -> LASSO-ready dataset
Adaptive_LASSO_MocoInverse      % fit and export lasso_controller_result_moco.mat
```

`moco_inverse_to_lasso_dataset.m` reconstructs `fiber_length_norm` /
`mtu_force_norm` with `human0918.osim`, rebuilds stance/swing from GRF, and
cleans the raw excitation (Hampel despike + zero-phase low-pass).  The formal
MocoInverse result must contain excitation; if it is missing the interface
errors instead of reconstructing it.  See
[`docs/MOCOINVERSE_LASSO_INTERFACE.md`](docs/MOCOINVERSE_LASSO_INTERFACE.md) for
the short contract and the list of new `.m` files, and
[`docs/MOCOINVERSE_TO_LASSO_WORKFLOW.md`](docs/MOCOINVERSE_TO_LASSO_WORKFLOW.md)
for the detailed end-to-end workflow, equations, parameter tables and
troubleshooting.

### 7. 多随机种子一致性结构辨识（stability selection）

```matlab
cd('Sparse Group LASSO Validation')
Adaptive_LASSO_MocoInverse_Stability
```

同一份实验数据、同一套辨识程序（超参不变），用 `nSeeds`（默认 20）个随机
种子对**步态片段 (run)** 做有放回重采样，各跑一次 adaptive LASSO，得到
K 个 26×9 的 0-1 结构。按逐条目入选频率 π 做**多数投票**得到共识 mask——
"20 个种子里 13 个收敛到同一框架就选它"正是 π ≥ 0.5 的特例（13/20 > 0.5，
每个条目都被多数覆盖）。最后在**全量数据**上冻结 mask、只重新估计系数。

- 0-1 图是**硬约束**：`unpack_lasso_reflex_params` 每帧强制 `beta(~mask)=0`，
  CMA-ES 只能改有值条目的数值；
- 系数只是 CMA-ES 的 **initial guess**（real2sim gap 交给优化吸收）；
- 若最大一致簇占比低于 `minDominantFraction`（默认 0.6），说明这份数据
  没有唯一可辨识的结构，脚本会报警而不是硬选一个。

输出 `results/lasso_controller_result_moco_consensus.mat`（可直接替换
`lassoFile`）与 `results/lasso_stability_diagnostics.mat`（π、Jaccard、
聚类、每个种子的 mask……）。详见
[`docs/LASSO_STABILITY_IDENTIFICATION.md`](docs/LASSO_STABILITY_IDENTIFICATION.md)。

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
| `results/` | `lasso_controller_result_moco_consensus.mat` | Consensus (stability) controller; frozen 0-1 map + `identification`/`stability` metadata |
| `results/` | `lasso_stability_diagnostics.mat` | Per-seed masks, π, Jaccard/clusters, frozen-mask refit + unfrozen reference fit |
| `results/` | `lasso_stability_checkpoint.mat` | Resume point of the multi-seed identification |
| `results/` | `lasso_stability_*.png` | Stability diagnostics figures |
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

## Recent Changes (v2.3)

v2.3 adds an **experimental-data-informed front end** on top of the v2.2
predictive core.  The core itself is unchanged: `human0918` (18 muscles,
9/leg), 26-feature grouped LASSO controller, IPOP-CMA-ES, checkpoint/resume.

- **MocoInverse integration** (`experiments/moco/`): Python/OpenSim Moco
  scripts now solve MocoInverse on measured Vicon/GRF data, optionally with
  sEMG tracking, and export full solutions (activation/excitation,
  residual/reserve actuators), summaries, IK/GRF windows and segment-level
  diagnostics.  Key scripts are `moco_inverse.py`, `moco_env_check.py` and
  `read_moco_results.m`.
- **sEMG tracking with adjustable scale factors**: added the reverse-scaled
  B-form control tracking goal `(e - s·u)²` (with `g = 1/s`) alongside the
  original A-form `(u - g·e)²`, implemented by the
  `MocoReverseScaledControlTrackingGoal` plugin.  Each tracked muscle has
  its own scale `s`, which can be fixed (`--emg-s-values`), initialised from
  a previous solution, or updated by the outer closed-form loop
  (`s_i = Σ e_i u_i / Σ u_i²`, clipped to 0.5–2).  This makes the
  EMG-informed trajectory less sensitive to the EMG/activation amplitude
  mismatch and provides an explicit, interpretable scaling parameter.
- **MocoInverse → LASSO interface**: new `.m` pipeline
  [`moco_inverse_to_lasso_dataset.m`](<Sparse Group LASSO Validation/moco_inverse_to_lasso_dataset.m>)
  converts Moco results + IK + GRF into a SCONE-style LASSO dataset:
  - requires complete excitation and activation for all 18 muscles
    (missing channels raise `MissingExcitation` / `MissingActivation`,
    not NaN filling);
  - reconstructs `fiber_length_norm` / `mtu_force_norm` with
    `human0918.osim` and sub-stepped fibre-length integration;
  - rebuilds the exact `cal_gait_phase.m` 0..4 phases from 0.02·BW GRF plus
    pelvis/calcn COM (no hysteresis or percentage heuristics);
  - cleans raw excitation with Hampel despiking + zero-phase low-pass;
  - writes `moco_lasso_dataset.sto`, `moco_lasso_dataset.mat`, QC figures
    and `sampleWeight` for seam/invalid samples.
- **New Moco-specific LASSO main**:
  [`Adaptive_LASSO_MocoInverse.m`](<Sparse Group LASSO Validation/Adaptive_LASSO_MocoInverse.m>)
  fits the per-phase adaptive LASSO with sample weights and Huber IRLS
  (both on the initial OLS and on the post-selection refit), and writes
  `results/lasso_controller_result_moco.mat`.  The old
  `UN.sto → Adaptive_LASSO_StanceSwing.m` path is untouched.
- **Model update (2026-11)**: all `human0918.osim` muscle
  `max_isometric_force` values were doubled, now matching
  `SQR_simbody.osim`.  Code reads Fmax dynamically; normalized `lCEN` /
  `fATN` features are invariant under a uniform Fmax scaling, but the
  absolute force capacity changes forward dynamics and static optimisation.
  Old checkpoints/optimisation results should not be resumed and CMA-ES
  must be retuned.  See
  [`docs/MOCOINVERSE_TO_LASSO_WORKFLOW.md`](docs/MOCOINVERSE_TO_LASSO_WORKFLOW.md)
  §14.
- **Validation and documentation**: added
  [`test_moco_feature_reconstruction.m`](tests/test_moco_feature_reconstruction.m)
  for the fibre-length reconstruction, plus
  [`docs/MOCOINVERSE_TO_LASSO_WORKFLOW.md`](docs/MOCOINVERSE_TO_LASSO_WORKFLOW.md)
  and [`docs/MOCOINVERSE_LASSO_INTERFACE.md`](docs/MOCOINVERSE_LASSO_INTERFACE.md)
  for the end-to-end workflow, parameter tables and troubleshooting.

## Recent Changes (v2.4) — 多随机种子一致性结构辨识

v2.4 在 LASSO 前端增加稳定性选择（stability selection），解决“0-1 结构到底
可不可信”的问题。核心（`human0918`、26-feature grouped controller、
IPOP-CMA-ES、checkpoint/resume）不变。

- **新功能**：[`Adaptive_LASSO_MocoInverse_Stability.m`](<Sparse Group LASSO Validation/Adaptive_LASSO_MocoInverse_Stability.m>)
  用 `nSeeds` 个随机种子对步态片段做 run-level bootstrap 重采样，各跑一次
  同一套 adaptive LASSO，按逐条目入选频率 π 做多数投票得到共识 0-1 结构，
  再在全量数据上冻结 mask 重拟合系数。输出
  `results/lasso_controller_result_moco_consensus.mat`。
- **稀疏度分析**：[`analyze_lasso_sparsity.m`](<Sparse Group LASSO Validation/analyze_lasso_sparsity.m>)
  量化“多一个非零条目换来多少误差”（lambda path 边际收益 + π 阈值扫描 +
  `tolSelect` 扫描），并用 Kneedle 判据给出膝点：**λ path 膝点 = 2–3 项/肌肉**，
  **π 阈值膝点 = 175 参数（π≥0.70）**。结论：π≥0.7 的 stability 选择给出
  **175** 个参数（比 majority 的 247 少 29%，留出 RMSE +8.5%），见
  `results/lasso_sparsity_tradeoff.png` 与
  `results/lasso_controller_result_moco_consensus_pi70.mat`。
- **拟合效果可视化**：[`plot_lasso_fit_quality.m`](<Sparse Group LASSO Validation/plot_lasso_fit_quality.m>)
  画出右腿 9 块肌肉的实测 vs 预测 excitation（全程 30 s、3 s 放大、逐肌肉 RMSE 柱状），
  175 与 247 参数的预测曲线肉眼几乎重合。另有
  [`plot_controller_fit_quality.m`](<Sparse Group LASSO Validation/plot_controller_fit_quality.m>)：
  直接读取**任意** controller `.mat`（不需要 diagnostics 文件）重建预测并出同样三张图。
- **重采样对照**：[`compare_bootstrap_resampling.m`](<Sparse Group LASSO Validation/compare_bootstrap_resampling.m>)
  证明 i.i.d. 重采样会把稳定性系统性夸大（Jaccard 0.67/0.75 vs run-level 0.53/0.59），
  但不改变 nnz。
- **重构**：`fit_adaptive_lasso_phase_model.m` 拆成
  [`load_lasso_design_matrices.m`](<Sparse Group LASSO Validation/load_lasso_design_matrices.m>)
  + [`fit_adaptive_lasso_core.m`](<Sparse Group LASSO Validation/fit_adaptive_lasso_core.m>)，
  原函数保留为薄 wrapper（签名、默认值、输出字段全部不变）。等价性由
  [`tests/check_lasso_refactor_equivalence.m`](tests/check_lasso_refactor_equivalence.m)
  用重构前基线逐项验证（容差 1e-12）。
- **`sgl_fit.m`**：每次求解的 summary 打印改为受 `opts.verbose` 控制
  （K×60 lambda 的集合运行不再刷屏；数值完全不变）。
- **测试**：[`tests/test_lasso_stability_identification.m`](tests/test_lasso_stability_identification.m)
  覆盖“13/20 就选它”的多数投票性质、π=0.5 tie-break、bootstrap 可复现性、
  合成数据的端到端结构不变量与 loader round-trip。
- **文档**：[`docs/LASSO_STABILITY_IDENTIFICATION.md`](docs/LASSO_STABILITY_IDENTIFICATION.md)
  （参考手册）、
  [`docs/STABILITY_IDENTIFICATION_EXPLAINED.md`](docs/STABILITY_IDENTIFICATION_EXPLAINED.md)
  （从零讲起：词汇、手算例子、每一步的输入输出、参数手册、FAQ），
  [`docs/LASSO_SPARSITY_ANALYSIS_ARCHIVE.md`](docs/LASSO_SPARSITY_ANALYSIS_ARCHIVE.md)
  （论文素材：稀疏度曲线、膝点、与 Meunier 1993 的关系），
  [`docs/REFLEX_CONTROLLER_INTERPRETATION.md`](docs/REFLEX_CONTROLLER_INTERPRETATION.md)
  （`y = Xβ + b` 的生理学解释：反馈通道 + CPG + 相位调制）。

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