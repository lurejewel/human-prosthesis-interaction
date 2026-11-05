# MocoInverse → LASSO → 预测动力学：接口与预处理说明

本文档说明新增的实测数据链路：

```text
Vicon / GRF / sEMG
        │
        ▼
MocoInverse（experiments/moco/moco_inverse.py）
        │  肌肉 activation / excitation
        │  ik_radians.mot（运动学）
        │  grf_window_*.mot（GRF）
        ▼
moco_inverse_to_lasso_dataset.m        ← 新接口（本文档）
        │  moco_lasso_dataset.sto
        │  moco_lasso_dataset.mat
        ▼
Adaptive_LASSO_MocoInverse.m           ← 新 LASSO main
        │  results/lasso_controller_result_moco.mat
        ▼
demo_predictiveForwardSimulation_humanModel.m（CMA-ES 预测仿真）
```

旧的 `UN.sto → Adaptive_LASSO_StanceSwing.m` 路径完全保留，未被修改。

---

## 1. 为什么需要新接口

旧的 LASSO 输入 `Sparse Group LASSO Validation/UN.sto` 来自 SCONE，已经包含
LASSO 需要的全部列：

```text
<muscle>_r.excitation
<muscle>_r.fiber_length_norm
<muscle>_r.mtu_force_norm
hip_flexion_r, knee_extension_r, ankle_dorsiflexion_r, ...
leg1_r.state
```

MocoInverse 的输出不是这个格式：

- 初步 `--sagittal-cost` 解只有 18 个肌肉 **activation** state，没有坐标、没有
  excitation control；
- 即使完整解包含 excitation，也没有 `fiber_length_norm` / `mtu_force_norm`
  （Moco 使用了 `ModOpIgnoreTendonCompliance`，纤维长度不是轨迹 state）；
- 没有 SCONE 的 `leg1_r.state`，步态相位必须从 GRF 重建。

新接口补上这三类缺失，并**按预测仿真模型 human0918 的定义**计算特征。
这一点很关键：LASSO 学到的控制器最终在 `human0918.osim` 里运行，所以训练特征
必须和 `extract_lasso_reflex_features.m` 在运行时看到的量一致。

---

## 2. 新文件

| 文件 | 作用 |
|------|------|
| `Sparse Group LASSO Validation/moco_inverse_to_lasso_dataset.m` | 接口主脚本：读 Moco 解、IK、GRF，清洗 excitation，算特征，写 LASSO-ready `.sto` / `.mat` |
| `Sparse Group LASSO Validation/Adaptive_LASSO_MocoInverse.m` | 新 LASSO main；输出 `results/lasso_controller_result_moco.mat` |
| `read_sto_mot_table.m` | 纯 MATLAB 的 `.sto/.mot` 读取器（不依赖 OpenSim MATLAB 接口） |
| `read_moco_inverse_solution.m` | 从 Moco 解中解析 `/forceset/<muscle>`、`/activation`、`/jointset/.../value|speed` |
| `clean_excitation_signals.m` | despike + 零相位低通；target 始终是 cleaned raw excitation |
| `zero_phase_lowpass.m` | 零相位 Butterworth（无 Signal Processing Toolbox 时退化为移动平均） |
| `compute_com_x_from_kinematics.m` | 计算 pelvis/calcn 质心全局 X，供 `cal_gait_phase` 规则使用 |
| `detect_stance_from_grf.m` | 严格按 `cal_gait_phase.m` 重建 0..4 相位（0.02 BW 单阈值） |
| `build_lasso_features_from_kinematics.m` | 沿给定运动积分纤维长度状态，计算 lCEN / fATN |
| `write_lasso_dataset_sto.m` | 写 SCONE 风格的 `.sto` |
| `fit_adaptive_lasso_phase_model.m` | 可复用的 per-phase adaptive LASSO 核心（含样本权重与 Huber IRLS） |
| `plot_moco_lasso_dataset_qc.m` | 数据 QC 图 |

---

## 3. 运行方法

### 3.1 生成 LASSO 数据集

编辑 `moco_inverse_to_lasso_dataset.m` 顶部的 `CONFIG`：

```matlab
solutionFile = fullfile(repoRoot,'experiments','moco','outputs', ...
    'moco_inverse_solution_B30_right_final.sto');
ikFile       = fullfile(repoRoot,'experiments','moco','outputs','ik_radians.mot');
grfFile      = fullfile(repoRoot,'experiments','moco','outputs', ...
    'grf_window_209p5_240.mot');
modelFile    = fullfile(repoRoot,'model','human0918.osim');
targetDt     = 0.01;
```

> `solutionFile` 必须包含 excitation（`/forceset/<muscle>` 控制列）。
> 如果检测不到 excitation，接口会直接报 `NoExcitation`，不会从 activation
> 反演。正式的 MocoInverse 结果应满足这个条件。

然后：

```matlab
cd('Sparse Group LASSO Validation')
moco_inverse_to_lasso_dataset
```

输出在 `Sparse Group LASSO Validation/moco_dataset/`：

```text
moco_lasso_dataset.sto
moco_lasso_dataset.mat
moco_lasso_dataset_qc.csv
moco_lasso_dataset_qc_drives.png
moco_lasso_dataset_qc_phase.png
moco_lasso_dataset_qc_features.png
moco_lasso_dataset_report.txt
```

### 3.2 拟合 LASSO 控制器

```matlab
cd('Sparse Group LASSO Validation')
Adaptive_LASSO_MocoInverse
```

输出：

```text
results/lasso_controller_result_moco.mat
results/lasso_moco_fit_diagnostics.mat
Sparse Group LASSO Validation/Adaptive_LASSO_MocoInverse.log
```

`lasso_controller_result_moco.mat` 的 `lasso` 变量与旧格式相同（grouped v2），
`load_lasso_reflex_controller.m` 可直接读取。把
`demo_predictiveForwardSimulation_humanModel.m` /
`prepare_fresh_start.m` 里的 `lassoFile` 指向这个文件即可。

> 默认不覆盖旧的 `results/lasso_controller_result.mat`，避免破坏 UN.sto 基准。

---

## 4. 特征是怎么算出来的

`build_lasso_features_from_kinematics.m` 复现 `ModelInfo.update` 的定义：

```matlab
lCEN = muscle.getNormalizedFiberLength(state);
fATN = muscle.getFiberForceAlongTendon(state) / muscle.getMaxIsometricForce();
```

给定 Moco 的 activation 轨迹和 IK 的 q/qdot：

1. 把 q、qdot、activation 插值到统一的 100 Hz 网格；
2. 用 `model.setStateVariableValue` 写入 activation；
3. 以 `lopt` 为初值，`equilibrateMuscles` 求初始纤维长度；
4. 对每个输出时间步，用 ≤2 ms 的子步长显式积分纤维长度状态：

   ```text
   d(fiber_length)/dt = model.computeStateVariableDerivatives(...)
   fiber_length <- fiber_length + dt_sub * d(fiber_length)/dt
   ```

5. 在真实输出时刻计算 lCEN / fATN。

验证结果（对 `results/sim_result_2026-07-27_15-39-45.sto` 的已知 fiber-length
state）：子步积分与真实状态的 RMSE < 1e-3（归一化长度），而准静态平衡近似在
rect_fem / vasti 上误差可达 0.13。因此这里没有用“每帧重新平衡”的简化。

> 默认用 `human0918.osim`，因为它是预测仿真实际使用的模型。若想用
> `SQR_simbody.osim` 的特征，把 `modelFile` 改过去即可；只要 qNames 能匹配，
> 其余流程不变。
>
> 注意：`human0918.osim` 的 `max_isometric_force` 现已调整为与
> `SQR_simbody.osim` 相同的值（例如 hamstrings 5188 N），但两者的
> `optimal_fiber_length` / tendon slack / 几何仍略有差异。`excitation` 是
> 无量纲量；`fATN = 肌肉力 / 该模型 Fmax`，在固定 activation/kinematics 下
> **均匀缩放 Fmax 不改变 fATN**（肌肉力也按同比例缩放），真正改变 fATN 的
> 是 `lopt`/tendon slack/几何等非均匀差异。控制器最终在 human0918 里运行，
> 因此默认必须用 human0918 计算特征；只有当你想分析 SQR 模型自身的控制律
> 时才改 `modelFile`。
>
> **Fmax ×2 后的重要提醒**：`UN.sto` 的 `mtu_force_norm` 和旧 LASSO 控制器
> 的 fATN 行**不需要**缩放（归一化特征不变）。但绝对力容量翻倍会改变前向
> 动力学和静态优化，旧的 `checkpoint.mat` / `opt_result_*.mat` 中的 `a_opt`、
> fitness、CMA-ES 状态都对应旧模型，应删掉/备份后从 fresh start 重新优化；
> 旧 LASSO 控制器可作为 CMA-ES 初值，但需要重新调参。详见
> [`MOCOINVERSE_TO_LASSO_WORKFLOW.md`](MOCOINVERSE_TO_LASSO_WORKFLOW.md)
> 第 14 节。

---

## 5. 噪声 / burst / 不连贯的处理

实测数据进入 MocoInverse 后，excitation 常见三类问题：

1. **高频 chatter**：collocation 网格上的交替分量；
2. **短 burst / spike**：优化器利用 activation 快速通道产生的非生理尖峰；
3. **x 与 a 不一致**：`x` 是 control，`a` 是 state，二者只通过
   `da/dt = (x-a)/tau` 联系；直接拟合 raw x 会把数值伪影学进控制器。

### 5.1 输入侧（接口里做）

`clean_excitation_signals.m` 的顺序：

1. 非有限值插值；
2. **Hampel/中值 despike**：只替换长度 ≤ `maxBurstSamples`（默认 5 个样本，
   100 Hz 下 50 ms）的短 run，持续的真实峰值保留。因为很多 activation 样本
   停在 0.01 下限，鲁棒尺度会退化成 0，所以用 `despikeMinSigma`（默认 0.005）
   设一个绝对噪声地板，避免把正常的小幅波动全判成 spike；
3. **零相位低通**：默认 6 Hz、4 阶 Butterworth（`filtfilt`，无相移）；
4. 目标构造：**只使用清洗后的 raw excitation**；
   - 如果 `xRaw` 为空，或某块肌肉的 excitation 整列缺失，直接报错；
   - 不再提供 `activation_consistent` / `activation` / `auto` 等回退方案；
5. 裁剪到 `[0.01, 1.0]`，其中 0.01 与 LASSO 的 censoring floor `c` 对齐。

对 activation / lCEN / fATN 还额外做 6 Hz 低通（`featureFilterHz`），避免
`fATN` 继承 activation 的高频噪声。

### 5.2 LASSO 侧（新 main 里做）

1. **样本权重**：`moco_lasso_dataset.mat` 的 `ds.sampleWeight` 把
   Moco 分段接缝、首尾边缘、非有限样本标为 0；`sgl_fit` 新增
   `opts.sampleWeight`，加权同时作用于 uncensored 平方项和 censored hinge 项。
   这样不需要删除行，lambda path 仍然连续。
2. **Huber IRLS**：在选定 support 后的 post-selection refit 阶段，用
   `delta = 1.345 * 1.4826 * MAD(residual)` 计算 Huber 权重，迭代
   `options.robustIter`（默认 2）次。短 burst 的残差会被降权，不会把系数拉偏。
3. **自适应权重 + 后选择 OLS**：保留旧流程，保证与 UN.sto 结果可比。
4. **lambda 选择**：仍用“最小 refit RMSE 的 (1+tolSelect) 邻域内取最稀疏”的
   规则，但 RMSE 是**加权** RMSE。

### 5.3 还可以怎么做（按需）

- 把 `cleanOpts.lowpassHz` 从 6 降到 4：更平滑，但可能削掉真实的快速峰值；
- 在 `excludeTimeWindows` 里显式排除分段接缝（如 B30 的
  217.5 / 225 / 232.5 s）或已知坏窗口；
- 把 `featureFilterHz` 设为 0，检查滤波对系数的影响；
- 用 `robustIter = 0` 与 `robustIter = 2` 对比，做消融实验。

---

## 6. 相位标签

`detect_stance_from_grf.m` 严格按 `functions/core/cal_gait_phase.m` 重建
0..4，不再使用滞回/短 run 清理/百分比启发式。

- 阈值只有一个：`th = 0.02 * bodyWeight`，与
  `forward_simulation.m` 的 `simCache.stanceTh = 0.02` 一致；
- 同时需要 `pelvis / calcn_r / calcn_l` 质心的全局 X 坐标，由
  `compute_com_x_from_kinematics.m` 用同一个模型和 q/qdot 计算；
- 0–4 的定义和状态转移完全照抄 `cal_gait_phase.m`；
- `stateR/stateL` 是严格 0..4；
- `stanceR/stanceL` 是 LASSO 的 stance 组（state 0/1）；
- `contactR/contactL` 是真实接触（GRF >= th，包含 phase 2 Liftoff）。

`extractStoMuscleFeatures.m` 用 state 0/1 作为 stance、2/3/4 作为 swing，
所以 phase 2 Liftoff 在 LASSO 里归入 swing 组，尽管它仍有 GRF。

如果 GRF 文件旁边有 ExternalLoads XML，会优先按
`<ExternalForce name="right">` 的 `force_identifier` 找列；否则按 SQR 约定
`ground_force2_vy = right`、`ground_force1_vy = left`。

阈值是体重的百分比，所以换受试者时要改 `phaseOpts.bodyWeight`（默认 78 kg
× 9.81）。如果出现 `invalidPhase`，应检查 GRF 列/符号、COM 计算模型和
坐标单位，而不是改阈值。

---

## 7. 快速检查清单

- [ ] `moco_lasso_dataset.sto` 里 `hamstrings_r.excitation` 等列存在且为 0.01–1；
- [ ] `lCEN` 在 ~0.5–1.4，`fATN` 在 ~0–1（均匀缩放 Fmax 不影响 fATN；负值可能来自 passive force，属正常）；
- [ ] `leg1_r.state` 只有 0–4；
- [ ] `moco_lasso_dataset_qc_drives.png` 里清洗后的 target 不再有尖刺；
- [ ] LASSO log 里每个肌肉/相位的 `nnz` 不为 0；
- [ ] `results/lasso_controller_result_moco.mat` 能被
      `load_lasso_reflex_controller.m` 读取（可用 `test_lasso_reflex_controller.m`
      的同款检查）。
