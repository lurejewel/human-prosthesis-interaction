# MocoInverse Results → LASSO Optimization：完整对接流程与实现细节

本文档描述从 **MocoInverse 实测数据结果** 到 **LASSO 肌肉反射控制器初值** 的
完整链路。它比 [`MOCOINVERSE_LASSO_INTERFACE.md`](MOCOINVERSE_LASSO_INTERFACE.md)
更详细，重点说明：

- 每一类输入文件里到底有什么、缺什么；
- 每个中间量（excitation、lCEN、fATN、相位、样本权重）是怎么算出来的；
- LASSO 的数学形式、lambda path、post-selection refit 和 robust IRLS；
- 代码中每个参数的默认值、物理含义和调节建议；
- 常见失败模式与排查方法。

## 目录

- [0. 一句话总览](#0-一句话总览)
- [1. 输入文件：MocoInverse 结果有哪几种形态](#1-输入文件mocoinverse-结果有哪几种形态)
- [2. 接口主流程：`moco_inverse_to_lasso_dataset.m`](#2-接口主流程moco_inverse_to_lasso_datasetm)
- [3. 特征定义与纤维长度重建](#3-特征定义与纤维长度重建)
- [4. 噪声、burst 与 excitation 目标构造](#4-噪声burst-与-excitation-目标构造)
- [5. 相位重建：`detect_stance_from_grf.m`](#5-相位重建detect_stance_from_grfm)
- [6. LASSO 优化：具体数学与实现](#6-lasso-优化具体数学与实现)
- [7. 输出文件与字段](#7-输出文件与字段)
- [8. 运行步骤](#8-运行步骤)
- [9. QC 与验收标准](#9-qc-与验收标准)
- [10. 常见问题排查](#10-常见问题排查)
- [11. 与旧 `UN.sto` 路径的关系](#11-与旧-unsto-路径的关系)
- [12. 设计取舍与局限](#12-设计取舍与局限)
- [13. 建议的消融/敏感性实验](#13-建议的消融敏感性实验)
- [14. 模型 Fmax 变更后的兼容性（2026-11）](#14-模型-fmax-变更后的兼容性2026-11)

**建议阅读顺序**：第一次接触这条链路时，先读第 0–2 节建立整体流程；
需要改代码/调参数时重点看第 3–6 节；实际运行时按第 8 节操作；
出问题先查第 9–10 节。

---

## 0. 一句话总览

```text
Vicon / GRF / sEMG
        │
        ▼
MocoInverse（experiments/moco/moco_inverse.py）
        │  activation / excitation
        │  ik_radians.mot（运动学）
        │  grf_window_*.mot（GRF）
        ▼
moco_inverse_to_lasso_dataset.m
        │  moco_lasso_dataset.sto   （SCONE 风格，供旧 reader 使用）
        │  moco_lasso_dataset.mat   （原始/清洗信号、特征、相位、权重、QC）
        ▼
Adaptive_LASSO_MocoInverse.m
        │  results/lasso_controller_result_moco.mat
        ▼
demo_predictiveForwardSimulation_humanModel.m
        │  CMA-ES 预测仿真
        ▼
kinematics / kinetics / biomechanics / 最优反射参数
```

旧的 `UN.sto → Adaptive_LASSO_StanceSwing.m` 路径完全保留；新链路只增加
MocoInverse 分支，输出文件名带 `_moco`，不会覆盖旧控制器。

---

## 1. 输入文件：MocoInverse 结果有哪几种形态

MocoInverse 输出并不是只有一种格式。当前项目里至少会出现以下三种：

| 形态 | 包含的列 | 典型来源 | 接口如何补 |
|---|---|---|---|
| **完整解** | `/forceset/<muscle>/activation`、`/forceset/<muscle>`、`/jointset/.../value`、`/jointset/.../speed` | 未加 `--sagittal-cost` 的 Moco 解 | 直接用解里的坐标/速度；有 excitation 时可用 raw control |
| **只有 activation** | 只有 18 个 `/forceset/<muscle>/activation` | 初步 `--sagittal-cost` 解 | **不支持**：接口需要 excitation，缺省会报 `MissingExcitation`；请提供含 excitation 的最终解 |
| **只有 excitation** | 只有 `/forceset/<muscle>` | `moco_inverse_excitations_*.sto` | **不完整**：必须另给 activation 文件；不会由 x 积分补 a |

`read_moco_inverse_solution.m` 的列名匹配规则：

```text
activation : ^/?forceset/([^/]+)/activation$
excitation : ^/?forceset/([^/]+)$
             但排除 residual_* / reserve_*
坐标       : ^/?jointset/[^/]+/([^/]+)/(value|speed)$
```

注意：

- 裸列名 `/forceset/soleus_r` 是 **excitation control**，不是 activation；
- `/forceset/soleus_r/activation` 是 **activation state**；
- 骨盆残差/储备作动器也长成 `/forceset/residual_...`，必须排除，否则会被
  误当成肌肉兴奋；
- Moco 的 `/jointset/.../value` 已经是 SimTK 内部单位（rad / m），**不要**
  再根据 `inDegrees` 转换；只有 IK 文件才需要 deg → rad。

### 1.1 除 Moco 解之外还需要什么

| 文件 | 何时需要 | 作用 |
|---|---|---|
| `ik_radians.mot` | Moco 解没有 `/jointset/.../value` 时 | 提供 q/qdot；已经 15 Hz 滤波、单位 rad |
| `grf_window_*.mot` | 总是需要 | 重建 stance/swing；不参与肌肉特征 |
| `grf_window_*.xml` | 推荐 | 告诉接口 `ground_force2_*`/`ground_force1_*` 哪个是右/左 |
| `model/human0918.osim` | 总是需要 | 计算 lCEN/fATN；必须是预测仿真实际使用的模型 |

注意：`--sagittal-cost` 的初步结果通常只有 `activation + ik + grf`，**没有
excitation**。这种结果不再被接口接受：运行时会直接报 `MissingExcitation`。
正式 MocoInverse 结果必须包含 `/forceset/<muscle>` 的 excitation control。

---

## 2. 接口主流程：`moco_inverse_to_lasso_dataset.m`

主脚本按 8 个阶段执行。下面按代码顺序说明。

### 2.1 路径与日志

```matlab
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
outDir    = fullfile(scriptDir, 'moco_dataset');
diary(fullfile(outDir, 'moco_lasso_dataset_report.txt'));
```

所有输出写到：

```text
Sparse Group LASSO Validation/moco_dataset/
```

### 2.2 读 Moco 解与可选额外文件

```matlab
S = read_moco_inverse_solution(solutionFile, muscleNames);
```

`muscleNames` 固定为 canonical 顺序：

```text
hamstrings, bifemsh, glut_max, iliopsoas, rect_fem,
vasti, gastroc, soleus, tib_ant
```

先右腿 9 块，再左腿 9 块。LASSO 只拟合右腿 9 块，并假设双侧对称；左腿列
仍然写入 `.sto`，方便以后做双侧分析。

如果 `solutionFile` 完全没有 excitation/activation，可以用
`excitationFile` / `activationFile` 作为完整替代；接口不会按列合并，
也不会用 NaN 填充。主解/替代文件都必须覆盖全部 18 块肌肉。

### 2.3 公共时间网格

所有输入必须取交集：

```matlab
starts = max(各文件起始时间);
ends   = min(各文件结束时间);
tOut   = (starts:targetDt:ends).';
```

默认 `targetDt = 0.01 s`（100 Hz）。原因：

- 线性反射模型的主要带宽在 1–10 Hz；
- 100 Hz 对 10 ms 反射延迟仍有 1 个样本的分辨率；
- 纤维长度积分的 Java/OpenSim 调用次数与采样率成正比，100 Hz 在
  30–60 s 数据上仍然可控；
- LASSO 的样本数不会因为 1 kHz 而膨胀 10 倍。

如果 Moco 解自带 `/jointset` 坐标，则不需要把 IK 纳入交集；否则 IK 必须
覆盖整个时间窗。

### 2.4 运动学：Moco 坐标或 IK

两条分支：

1. **Moco 解有 `/jointset/.../value|speed`**
   - 直接插值到 `tOut`；
   - 单位已经是 rad/m；
   - 如果某些 `/speed` 缺失，用 q 的中心差分补。

2. **Moco 解没有 `/jointset` 坐标列**
   - 读 `ik_radians.mot`；
   - 检查 `inDegrees`：只对非平移坐标做 `deg2rad`；
   - qdot 用 `gradient(q, t)` 计算；
   - 再插值到 `tOut`。

接口会显式检查 `qOut/qdOut` 是否含 NaN/Inf；有则报 `BadKinematics`，不会
把 NaN 传进 OpenSim。

### 2.5 清洗 muscle drives

调用：

```matlab
[xClean, aClean, qcClean] = clean_excitation_signals(tOut, xRaw, aRaw, cleanOpts);
```

目标只有一个来源：**清洗后的 raw excitation**。

- 接口要求 `xRaw` 和 `aRaw` 都非空，且每块肌肉都有 finite excitation/activation；
- 缺失 excitation → `MissingExcitation`；缺失 activation → `MissingActivation`；
- 不用 NaN 填充，也不从另一个信号反演/积分补全；
- `aClean` 只用于下一步计算 `fATN/lCEN`，不参与 target 构造。

详见第 4 节。

### 2.6 肌肉特征 lCEN / fATN

调用：

```matlab
F = build_lasso_features_from_kinematics(modelFile, tOut, qOut, qdOut, ...
        aClean, muscleNames, qNamesOut, featureOpts);
```

默认 `modelFile = model/human0918.osim`。详见第 3 节。

### 2.7 相位

调用：

```matlab
comX = compute_com_x_from_kinematics(modelFile, tOut, qOut, qdOut, qNamesOut);
P = detect_stance_from_grf(grfFile, tOut, comX, phaseOpts);
```

得到严格按 `cal_gait_phase.m` 规则计算的 `stateR/stateL`（0..4），
以及 LASSO 用的 `stanceR/stanceL`（state 0/1）和真实接触
`contactR/contactL`（GRF >= 0.02 BW）。详见第 5 节。

### 2.8 有效样本与 `sampleWeight`

```matlab
validMask = true(nT, 1);
validMask(1:edge) = false;          % 首尾 0.05 s
validMask(end-edge+1:end) = false;
validMask = validMask & all(isfinite([xClean, aClean, F.lCEN, F.fATN]), 2);
% 再排除用户给的 excludeTimeWindows
sampleWeight = double(validMask);
```

`sampleWeight` 传到 LASSO 后，`sgl_fit` 会把它同时乘到 uncensored 平方项和
censored hinge 项；**不删除行**，所以 lambda path 的样本顺序和相位索引保持
连续。

### 2.9 输出

```text
moco_lasso_dataset.sto
moco_lasso_dataset.mat
moco_lasso_dataset_qc.csv
moco_lasso_dataset_qc_drives.png
moco_lasso_dataset_qc_phase.png
moco_lasso_dataset_qc_features.png
moco_lasso_dataset_report.txt
```

---

## 3. 特征定义与纤维长度重建

### 3.1 预测仿真运行时的 26 个特征

`extract_lasso_reflex_features.m` 在运行时使用的顺序：

| 行 | 特征 | 来源 |
|---|---|---|
| 1–9 | 归一化纤维长度 `lCEN` | 同侧 9 块肌肉 |
| 10–18 | 归一化 MTU 力 `fATN` | 同侧 9 块肌肉 |
| 19 | `pelvis_tilt` | 骨盆倾角 |
| 20 | `pelvis_tilt_u` | 骨盆倾角速度 |
| 21 | `hip_flexion_r/l` | 同侧髋角 |
| 22 | `hip_flexion_r/l_u` | 同侧髋角速度 |
| 23 | `knee_extension_r/l` | 同侧膝角 |
| 24 | `knee_extension_r/l_u` | 同侧膝角速度 |
| 25 | `ankle_dorsiflexion_r/l` | 同侧踝角 |
| 26 | `ankle_dorsiflexion_r/l_u` | 同侧踝角速度 |

肌肉列顺序固定为：

```text
hamstrings, bifemsh, glut_max, iliopsoas, rect_fem,
vasti, gastroc, soleus, tib_ant
```

`lasso.groups(g).beta` 是 26×9，列对应这 9 块肌肉。

### 3.2 为什么必须用 human0918 算特征

当前 `human0918.osim` 的 `max_isometric_force` 已经调整为与
`SQR_simbody.osim` 相同的值（例如 hamstrings 5188 N），但两者的
`optimal_fiber_length`、tendon slack 和部分几何仍略有差异。

- `excitation` 是无量纲量，可以跨模型比较；
- `fATN = F_fiber_along_tendon / Fmax` 是“该模型的归一化力”。在固定
  activation / kinematics 下，**均匀缩放 Fmax 不会改变 fATN**（肌肉力也按
  相同比例缩放）；真正会改变 fATN 的是 `lopt`、tendon slack、pennation 和
  肌肉路径等非均匀差异；
- 控制器最终在 `human0918.osim` 里运行，所以训练特征必须用 human0918 的
  肌肉参数计算，否则 LASSO 学到的系数在运行时特征空间里不成立。

> **注意**：如果以后再次修改 `Fmax` / `lopt` / tendon slack，`fATN` 的
> 数值尺度会改变，已有 LASSO 控制器和旧数据都必须重新生成/重新拟合。
> 详见第 14 节。

因此默认：

```matlab
modelFile = fullfile(repoRoot, 'model', 'human0918.osim');
```

只有当你明确要分析 SQR 模型自身的控制律时，才改 `modelFile`。

### 3.3 为什么不能每帧静力平衡

先澄清一个容易混淆的点：`Millard2012EquilibriumMuscle` 虽然叫
“equilibrium muscle”，但它的 `fiber_length` 仍然是**带记忆的状态**，不是
每帧由 q/activation 代数决定的量。

它的平衡方程是：

$$
F_{iso}\Big(a\,f_L(\hat l_{CE})\,f_V(\hat v_{CE})
+ f_{PE}(\hat l_{CE}) + \beta \hat v_{CE}\Big)\cos\phi
= F_{iso}\,f_{SE}(\hat l_T)
$$

OpenSim 并不是每帧去解这个代数方程，而是把它整理成纤维速度的 ODE。下面
写的是**无纤维阻尼**（`beta = 0`）的形式；有阻尼时方程里多一个
`beta * v_CE` 项，仍然要解出 `v_CE` 再积分：

$$
\hat v_{CE}
= f_V^{-1}\left(
\frac{f_{SE}(\hat l_T)/\cos\phi - f_{PE}(\hat l_{CE})}
     {a\,f_L(\hat l_{CE})}
\right)
$$

然后对 `fiber_length` 做数值积分：

$$
l_{CE}(t) = l_{CE}(0) + \int_0^t v_{CE}(\tau)\,d\tau
$$

所以 forward simulation 里的 `l_CE(t)` 取决于**从初始时刻到现在的整段
activation / q / qdot 历史**，而不是只取决于当前帧。

#### `equilibrateMuscles` 到底做了什么

`equilibrateMuscles(state)` 做的是：

> 在当前 q、qdot、activation 下，解上面的**代数平衡方程**，找到“如果
> q 和 activation 从此不再变化，纤维最终会稳定到的长度”。

它相当于把纤维长度当成“瞬时平衡”的结果，丢掉了之前的速度/长度历史。
在慢速、准静态动作里这没什么问题；但步行中 activation 和关节角一直在变，
纤维长度会**滞后**于瞬时平衡点。

#### 一个直观类比

把肌肉-肌腱想成“弹簧 + 阻尼器”：

- 你突然把关节转到一个新角度，肌腱张力不会瞬间跳到新平衡；
- 纤维长度需要一段时间“松弛”到新的平衡点；
- `equilibrateMuscles` 相当于假设“松弛速度无限快”，直接跳到最终平衡；
- 真实 forward simulation 有有限松弛速度，所以中途的长度和平衡长度不同。

#### 为什么 rect_fem / vasti 特别明显

这两个肌肉在测试的步态里力量大、activation/长度变化快，而且工作点落在
力-长度曲线较陡的区域：同样的力不平衡会对应更大的长度偏差。所以
“每帧平衡”与真实动态状态的差异在这两块肌肉上最大。

在 `results/sim_result_*.sto` 上实测（把每帧平衡得到的 lCEN 与文件里真实
`fiber_length/lopt` 比较）：

| 方法 | lCEN 误差（rect_fem / vasti） |
|---|---|
| 每帧准静态平衡 | RMSE 可达 0.13，max 0.34 |
| ≤2 ms 子步积分 | RMSE < 1e-3 |

#### 这对 LASSO 为什么重要

LASSO 的特征 `lCEN`、`fATN` 在运行时来自 forward simulation 的动态状态。
如果训练时用“每帧平衡”的 lCEN，就等于给 LASSO 喂了一个**控制器运行时
根本看不到的纤维长度**。对 rect_fem / vasti，这个误差不是随机噪声，而是
与 activation 相位相关的系统偏差；学到的 beta 在运行时会被应用到另一个
lCEN 上，导致 excitation 出现系统性偏差，尤其是在 stance/swing 过渡。

子步积分则是在每个子步里：

1. 设当前 q/qdot/activation；
2. 用模型自己的 ODE 算 `d(l_CE)/dt`；
3. 用小步长积分 `l_CE`。

它复现的正是 forward simulation 的“记忆”，只是 q/qdot/activation 来自
实测/Moco 数据而不是闭环仿真。因此特征和运行时一致。

> `equilibrateMuscles` 仍然用于第一帧：此时没有历史，解平衡方程是一个合理
> 初值；随后的 transient 会被接口的 `edgeMargin` 排除在 LASSO 之外。

### 3.4 子步积分算法

对每个相邻输出时刻 `[t_{k-1}, t_k]`：

```text
dt = t_k - t_{k-1}
nSub = ceil(dt / substepMax)          % 默认 substepMax = 0.002 s
for s = 1:nSub
    alpha = s / nSub
    q_sub  = (1-alpha) q_{k-1}  + alpha q_k
    qd_sub = (1-alpha) qd_{k-1} + alpha qd_k
    a_sub  = (1-alpha) a_{k-1}  + alpha a_k

    set coordinates(q_sub, qd_sub)
    set activation(a_sub)
    set fiber_length(fl)

    realizeDynamics
    computeStateVariableDerivatives
    dfl = getStateVariableDerivativeValue('/forceset/<muscle>/fiber_length')
    fl  = fl + dt/nSub * dfl
    fl  = clip(fl, 0.20*lopt, 2.50*lopt)
end
```

在真实输出时刻 `t_k`：

```matlab
lCEN(k, m) = muscle.getNormalizedFiberLength(state);
fATN(k, m) = muscle.getFiberForceAlongTendon(state) / muscle.getMaxIsometricForce();
```

初始纤维长度：

1. 第一帧 `fiber_length = lopt`；
2. `equilibrateMuscles(state)`；
3. 失败则逐肌肉 `computeEquilibrium(state)`；
4. 仍失败则保留 `lopt`。

首帧 transient 由接口的 `edgeMargin=0.05 s` 排除在 LASSO 之外。

### 3.5 特征再低通

`featureFilterHz = 6`：

- `lCEN` 和 `fATN` 都会再做一次零相位 6 Hz 低通；
- `fATN` 对 activation 的高频很敏感，低通能显著降低设计矩阵噪声；
- 步态主频约 1 Hz，6 Hz 不会削掉真实步态结构。

如果研究快速反射/跑步，可把 `featureFilterHz` 提高到 10–15 Hz，或设为 0
做敏感性分析。

### 3.6 特征验证

[`tests/test_moco_feature_reconstruction.m`](../tests/test_moco_feature_reconstruction.m)
用 `results/sim_result_*.sto` 的已知 `fiber_length` 验证：

- 只输入 q/qdot/activation；
- 重建 lCEN；
- 与真实 `fiber_length/lopt` 比较；
- 判据：max RMSE ≤ 0.01（正常应 ≲1e-3）。

---

## 4. 噪声、burst 与 excitation 目标构造

### 4.1 三类问题

MocoInverse 的 excitation 可能包含：

1. **高频 chatter**：collocation 网格尺度的交替分量；
2. **短 burst/spike**：优化器利用 10 ms activation 快速通道产生的非生理尖峰；
3. **x/a 不一致**：x 是 control，a 是 state，二者只近似满足
   `da/dt = (x-a)/tau`。

LASSO 用平方误差 + censoring hinge，少数大残差就能主导系数，因此必须先处理。

### 4.2 清洗顺序：每一步到底在做什么

`clean_excitation_signals.m` 的顺序是：

```text
① 非有限值按列插值
        ↓
② Hampel / 中值 despike（只替换短 run）
        ↓
③ 零相位 Butterworth 低通（默认 6 Hz）
        ↓
④ clip 到 [0.01, 1]
```

| 步骤 | 做什么 | 目的 | 关键参数 |
|---|---|---|---|
| ① 插值 | NaN/Inf 按列线性插值 | 让后续滤波有定义 | — |
| ② despike | 用滑动中值判断“局部异常点”，只替换短 run | 去 spike/burst，保留真实持续峰值 | `medianWindow=5`，`hampelSigma=4` |
| ③ 低通 | 4 阶 Butterworth，前后向滤波 | 去 chatter，且不引入相位延迟 | `lowpassHz=6` |
| ④ clip | 限制到 `[0.01, 1]` | 生理范围 + 与 LASSO censoring floor 对齐 | `floor=0.01` |

顺序不能随便换：

- **先 despike 再低通**：如果先低通，spike 会被抹成一个较宽的“鼓包”，
  Hampel 反而更难判断它是不是异常；
- **target 就是 raw excitation**：不再从 activation 反演，所以没有
  “先低通 a 再求导”这一步；`aClean` 只用于特征，不用于 target。

### 4.3 什么是 Hampel / 中值 despike

#### 4.3.1 一句话版本

Hampel 滤波 = **用一个滑动窗口的中值作为“局部正常水平”，用 MAD 作为
“局部正常波动尺度”，只把明显偏离这个正常水平的点判为 spike。**

“中值 despike” = 把判出来的 spike **替换成窗口中值**，而不是把整条信号
都做中值平滑。

#### 4.3.2 数学定义

对每个样本 `x_i`：

```text
1. 取窗口 W_i = {x_{i-k}, ..., x_i, ..., x_{i+k}}     % 默认 k=2，窗口 5 点
2. m_i = median(W_i)                                  % 局部正常水平
3. r_i = x_i - m_i                                    % 相对中值的残差
4. MAD_i = median(|r_j - median(r)|)                  % 残差的鲁棒尺度
5. sigma_i = 1.4826 * MAD_i                           % 换算成近似标准差
6. 若 |r_i| > nSigma * max(sigma_i, minSigma)          % 默认 nSigma=4
      则 x_i 判为 spike
```

几个关键点：

- **为什么用中值而不是均值？**  
  一个 0.95 的 spike 会把窗口均值拉高，但不会明显改变中值。中值对离群点
  不敏感，所以适合做“局部正常水平”。

- **为什么是 `1.4826 * MAD`？**  
  对高斯分布，`MAD ≈ 0.6745 * sigma`，所以乘 `1.4826` 后 `1.4826*MAD`
  近似等于标准差。这样 `nSigma=4` 就近似是“4 倍标准差”。

- **为什么还要 `minSigma=0.005`？**  
  activation 大量停在 0.01 下限，窗口里很多值完全相同，MAD 会变成 0。
  此时任何微小波动都满足 `|r_i| > 4*0`，flag 率可达 20–30%，把正常信号
  也当成 spike。`minSigma=0.005` 是一个绝对噪声地板：小于 0.5% activation
  的偏差不当作 spike。

#### 4.3.3 一个数值例子

单个 spike：

```text
原始:   0.10  0.10  0.10  0.95  0.10  0.10  0.10
窗口:         [0.10  0.10  0.95  0.10  0.10]
中值:                        0.10
残差:                        0.85   -> 远大于 4*sigma -> 判为 spike
替换后: 0.10  0.10  0.10  0.10  0.10  0.10  0.10
```

持续峰值：

```text
原始:   0.10  0.10  0.60  0.62  0.64  0.66  0.64  0.62  0.10
```

这里 0.60–0.66 连续 6 个样本都偏高，窗口内中值也被抬高，残差并不大；
即使被 flag，连续长度 6 > `maxBurstSamples=5`，也会被判定为“真实持续
峰值”而保留。

这就是“只替换短 run”的意义：

- 长度 ≤5 的连续 spike：替换成局部中值；
- 长度 >5 的持续高值：保留，不把真实的肌肉 burst 削掉。

#### 4.3.4 Hampel、普通中值滤波、低通的区别

| 方法 | 行为 | 副作用 |
|---|---|---|
| 普通中值滤波 | **每个点**都替换成窗口中值 | 真实峰值也被削平 |
| Hampel despike | 只替换“相对局部中值异常”的点 | 正常波动和持续峰值保留 |
| 低通滤波 | 对所有频率高于截止频率的成分衰减 | 会削弱快速但真实的峰值 |

所以流程是“Hampel 去短 spike → 低通去剩余 chatter”，而不是只用其中一种。

### 4.4 什么是零相位低通

`zero_phase_lowpass.m` 用的是 4 阶 Butterworth + `filtfilt`：

```text
正向滤波一次 -> 反向再滤波一次
```

- 普通 `filter` 是因果滤波，会引入相位延迟；对反射控制来说，相位延迟等价于
  人为改变反射延迟；
- `filtfilt` 正向+反向，净相位为 0，只改变幅频特性，不移动峰值时刻；
- 默认 6 Hz：步态主频约 1 Hz，6 Hz 保留步态谐波，去掉 Moco 的高频 chatter；
- 如果研究跑步/快速反射，可提高到 10–15 Hz；设为 0 可做“不滤波”的消融。

### 4.5 target 只用 raw excitation；没有 excitation 直接报错

#### 4.5.1 为什么不再用 activation-consistent

正式 MocoInverse 结果一定包含 excitation control，所以 LASSO 的 target
直接使用清洗后的 `x` 即可。`activation-consistent` 是“当解里只有 a、没有 x
时，用 `x = a + tau da/dt` 反推 x”的补丁；既然正式结果有 x，这个补丁不需要，
而且会引入额外假设。

因此当前逻辑是：

```text
excitation + activation 都完整 -> 清洗后拟合
缺 excitation              -> MissingExcitation，停止
缺 activation              -> MissingActivation，停止
不存在从 activation 反演 excitation 的路径
```

如果之前的 `results/lasso_controller_result_moco.mat` 是在
`activation_consistent` target 下拟合的，应重新跑
`Adaptive_LASSO_MocoInverse`，不要直接把它当作最终控制器。

#### 4.5.2 x 和 a 的关系仍然需要理解（用于 QC）

MocoInverse 同时有：

- `x(t)`：excitation，**control**；
- `a(t)`：activation，**state**。

二者由一阶激活动力学联系：

$$
\frac{da}{dt} = \frac{x - a}{\tau}
$$

其中上升用 `tau_act = 0.01 s`，下降用 `tau_deact = 0.04 s`。

虽然 target 不再由 a 反演，但 `aClean` 仍然用于计算 `fATN/lCEN`；
接口还会把 `xClean` 正向积分回 activation，做一致性检查：

```text
aPred = integrate_activation(xClean, tau_act, tau_deact)
activationReconRMSE = RMSE(aPred, aClean)
```

- 如果很小，说明清洗后的 x 和 Moco 的 a 在激活动力学下仍然一致；
- 如果很大，说明滤波/裁剪把 x-a 关系破坏得比较厉害，应该检查
  `lowpassHz` 或确认 Moco 解本身是否可信。

这只是 QC，不用于构造 target。

#### 4.5.3 程序行为

| 情况 | 行为 |
|---|---|
| `xRaw` 为空 | 报 `MissingExcitation` |
| 某块肌肉 excitation 整列缺失 | 报 `MissingExcitation`，并列出肌肉名 |
| `aRaw` 为空 | 报 `MissingActivation` |
| 某块肌肉 activation 整列缺失 | 报 `MissingActivation`，并列出肌肉名 |
| 任一信号含 NaN/Inf | 报 `NonFiniteExcitation` / `NonFiniteActivation` |
| `xRaw` 和 `aRaw` 都完整 | target 用 x；a 只用于特征 |

这里不会用 NaN 填充，也不会从 excitation 积分补 activation；正式
MocoInverse 解必须同时包含完整的 excitation 和 activation。

### 4.6 裁剪与 censoring floor

```matlab
xClean = min(max(xClean, 0.01), 1.0);
```

- 上限 1：excitation 不能超过 1；
- 下限 0.01：与 LASSO 的 `options.c = 0.01` 对齐。

LASSO 的 censored 模型把 `y <= c` 的样本视为“被 floor 在 c”，而不是普通
平方误差。如果 target 的下限和 `c` 不一致（例如 target 到 0，但 c=0.01），
censored 样本的统计假设就会错位。

### 4.7 QC 统计

`qc` 里保存：

- `targetSource` 固定为 `'excitation'`（仅作记录）；
- 每块肌肉 raw/clean 的均值、峰值；
- 被 flag 的样本数与百分比；
- `activationReconRMSE`：由 xClean 积分回 activation 与 aClean 的 RMSE；
- `flaggedExcitation/flaggedActivation` 逻辑矩阵。

这些会写进 `moco_lasso_dataset_qc.csv` 和 `.mat`。

---

## 5. 相位重建：`detect_stance_from_grf.m`

### 5.1 为什么用 GRF

旧 SCONE `.sto` 有 `leg1_r.state`，Moco 解没有。可选信号：

- 竖直 GRF：直接反映脚是否着地，最可靠；
- 脚高/踝角：受模型/IK 误差影响；
- 加速度微分：噪声大。

因此用竖直 GRF。

### 5.2 找到右/左竖直力列

优先级：

1. 用户显式 `rightColumn/leftColumn`；
2. 同名 ExternalLoads XML 的 `<ExternalForce name="right/left">`，
   `force_identifier` + `y`；
3. SQR 约定：`ground_force2_vy = right`、`ground_force1_vy = left`；
4. 列名含 right/left 或 `_r/_l`；
5. 最后兜底：第一个/第二个 `_vy` 列，并 warning。

### 5.3 单一阈值：0.02 BW，不用滞回

相位阈值必须与 `forward_simulation.m` 中的

```matlab
simCache.stanceTh = 0.02;   % 归一化到 body weight
```

一致。因此本接口使用：

```text
th = 0.02 * bodyWeight   [N]
GRF >= th  -> 与 stance/liftoff 相关的判断
GRF <  th  -> 与 swing/landing 相关的判断
```

不再使用 `stanceEnterBW / stanceExitBW` 双阈值，也不做
`minStance/minSwingGap` 的短 run 清理。相位完全由
`cal_gait_phase.m` 的状态转移规则决定。

`bodyWeight` 默认 `78 * 9.81`；换受试者必须改。

### 5.4 0..4 标签：严格按 `cal_gait_phase.m`

0..4 不是启发式百分比划分，而是 `functions/core/cal_gait_phase.m` 里定义的
**五相位状态机**。接口用 `compute_com_x_from_kinematics.m` 计算
`pelvis / calcn_r / calcn_l` 的质心全局 X 坐标，然后逐帧复现该函数。

#### 5.4.1 每帧定义

令 `th = 0.02 * bodyWeight`，`pelvisCOM`、`calcnCOM` 是全局 X 坐标：

| state | 名称 | 条件 |
|---|---|---|
| 0 | Early Stance | `GRF >= th` 且 `calcnCOM >= pelvisCOM` |
| 1 | Late Stance | `GRF >= th` 且 `calcnCOM < pelvisCOM` 且 opposite `GRF < th` |
| 2 | Liftoff | `GRF >= th` 且 `calcnCOM < pelvisCOM` 且（opposite `GRF > th` 或 `calcnCOM + 1 < pelvisCOM`） |
| 3 | Swing | `GRF < th` 且 `calcnCOM < pelvisCOM` |
| 4 | Landing | `GRF < th` 且 `calcnCOM >= pelvisCOM` |

注意 `calcnCOM + 1` 里的 `1` 是 **1 米**，与 `cal_gait_phase.m` 源码一致。

#### 5.4.2 后续帧转移

| 当前 state | 转移条件 | 下一 state |
|---|---|---|
| 0 Early Stance | `calcnCOM < pelvisCOM` | 1 Late Stance |
| 1 Late Stance | opposite `GRF >= th` 或 `calcnCOM + 1 < pelvisCOM` | 2 Liftoff |
| 2 Liftoff | own `GRF < th` | 3 Swing |
| 3 Swing | `calcnCOM > pelvisCOM` | 4 Landing |
| 4 Landing | own `GRF >= th` | 0 Early Stance |

第一帧则用 `cal_gait_phase.m` 里的 `if/elseif` 初始化规则（不是从 0 开始
硬推）。如果某一帧无法落入 0..4，接口会报 `invalidPhase`，而不是静默猜测。

#### 5.4.3 `stanceR` 与 `contactR` 的区别

- `P.stateR/stateL`：严格 0..4；
- `P.stanceR/stanceL`：**LASSO 的 stance 组**，即 `state 0/1`；
- `P.contactR/contactL`：`GRF >= th`，表示脚是否仍接触地面（包含 phase 2）。

`extractStoMuscleFeatures.m` 使用 `state 0/1` 作为 stance、`state 2/3/4`
作为 swing，所以 phase 2 Liftoff 在 LASSO 里归入 swing 组，尽管它仍然
满足 `GRF >= th`。这是旧 reader 的约定，不是新定义的启发式。

---

## 6. LASSO 优化：具体数学与实现

### 6.1 数据组织

`extractStoMuscleFeatures.m` 从 `.sto` 读出：

- `Y`：nSamples × 9 的右腿 excitation；
- `X`：nSamples × 26 的特征；
- `M`：mask（当前 identity）；
- `phase.stanceIdx/swingIdx`：来自 `leg1_r.state`。

然后按 `delaySteps` 做延迟对齐：

```text
target : y(t_k)        -> 去掉第一个样本
feature: A(t_{k-1})    -> 去掉最后一个样本
```

默认 `delaySteps = 1`，100 Hz 下就是 10 ms 反射延迟。

### 6.2 设计矩阵

```text
A = X * M
A_orig = A
muA = mean(A_orig, 1)
sA  = std(A_orig, 1, 1)
A   = (A_orig - muA) ./ sA
```

标准化只改变优化条件数；最终系数会反变换回原始特征尺度：

```text
beta_orig = beta_std ./ sA'
b_orig    = b_std - muA * beta_orig
```

### 6.3 censored 回归目标

`sgl_fit.m` 求解（每个 muscle/phase 独立）：

$$
\min_{\beta,b}\;
\frac12 \sum_{i \in U} w_i \left(y_i - (A_i\beta+b)\right)^2
+
\frac{\rho}{2} \sum_{i \in C} w_i \left((A_i\beta+b-c)_+\right)^2
+
\lambda_1 \sum_g \|\beta_g\|_2
+
\lambda_2 \|\beta\|_1
$$

其中：

- `U = {i: y_i > c + tolCensor}`，uncensored；
- `C = {i: y_i <= c + tolCensor}`，被 floor 在 `c` 的样本；
- `c = 0.01`；
- `(z)_+ = max(z,0)`，所以 censored 项只惩罚
  `A_i beta + b > c`，**不惩罚** `A_i beta + b <= c`；
- `rho = 50` 越大，越强制 censored 样本的 latent prediction
  `A_C beta + b <= c`，即把这一项压到 0；它不是强制 `>= c`；
- `w_i` 是 `sampleWeight`（0 = 排除，1 = 正常，Huber 时 0–1）；
- group 结构固定为：1–9 纤维长度、10–18 MTU 力、19–26 运动学。

> 注意：代码里计算 RMSE 时常用 `max(A beta + b, c)` 作为“报告值”，
> 因此打印/画图看起来总是 `>= c`；那是评估时的 clip，不是优化目标对
> censored 样本的要求。优化目标本身是让 censored 样本的 latent
> prediction 不超过 `c`。

在 adaptive LASSO 阶段：

```text
lambda1 = 0      % 不用 sparse-group，只用 L1
lambda2 = lambda % 沿 path 变化
```

### 6.4 adaptive 权重

先用 uncensored 样本做（加权）OLS 初估：

```text
beta_OLS = argmin sum_{U} w_i (y_i - A_i beta - b)^2
```

这个初始 OLS **默认也做 Huber IRLS**（`robustInitIter=2`）：

```text
1. 普通加权 OLS（一次闭式解）
2. 对 uncensored 残差做 Huber 权重
3. 再解一次加权 OLS
4. 重复 robustInitIter 次
```

为什么：

- 这一步只跑一次/每肌肉相位，计算量远小于 lambda path；
- 如果初始 OLS 被 burst 拉偏，adaptive 权重
  `w_j = 1/max(|beta_OLS_j|,eps)^gamma` 会跟着错；
- 先 robust 一下，adaptive 权重更稳；
- `robustInitIter=0` 可退回普通加权 OLS（对照实验）。

注意：这里只使用 uncensored 样本，且初始 OLS 的 Huber 权重只影响
`beta_OLS` → adaptive 权重；最终系数仍由后面的 lambda path + robust
post-selection refit 决定。

然后：

$$
w_j = \frac{1}{\max(|\beta_{OLS,j}|,\epsilon)^\gamma}
$$

默认 `gamma = 1`，`epsilon = 1e-8`。

列缩放：

```text
A_star = A ./ w
```

这样对 `beta_star` 的普通 L1 等价于对原始 `beta` 的加权 L1：

```text
lambda * sum_j w_j |beta_j|
```

### 6.5 lambda path

`lambda2Max` 由 `compute_lambda2_max_censored.m` 计算：使 `beta=0` 仍最优的
最小 `lambda`。

```text
alphaGrid = logspace(0, log10(1e-4), 60)
lambdaGrid = lambda2Max * alphaGrid
```

对每个 `lambda_k`：

1. `sgl_fit(A_star, y, 0, lambda_k)` → `beta_star`；
2. 换回 `beta_ada = beta_star ./ w`；
3. support `S_k = {j: |beta_ada_j| > nzTol}`，`nzTol=1e-5`；
4. 在 `S_k` 上做 post-selection refit（`lambda=0`），得到 `beta_refit`；
5. 计算加权 refit RMSE：

$$
\mathrm{RMSE}_k =
\sqrt{
\frac{\sum_i w_i (y_i - \max(A_i\beta_k+b_k,c))^2}
     {\sum_i w_i}
}
$$

### 6.6 lambda 选择

```text
minErr = min(RMSE_refit)
threshold = (1 + tolSelect) * minErr      % tolSelect=0.05
selected = 第一个 RMSE_refit <= threshold 的 lambda
```

也就是“在最小 RMSE 的 5% 容差内，取最稀疏（lambda 最大）的解”。这比单纯
取 min RMSE 更不容易过拟合，也符合旧 main 的行为。

### 6.7 robust Huber IRLS

#### 6.7.1 什么是 Huber loss

普通最小二乘对残差 `r` 的惩罚是 `0.5 r^2`，一个很大的 burst 残差会产生
极大的梯度，把系数拉偏。

Huber loss 是“二次 + 线性”的组合：

$$
\rho_\delta(r)=
\begin{cases}
\frac12 r^2, & |r| \le \delta \\
\delta |r|-\frac12\delta^2, & |r| > \delta
\end{cases}
$$

- 小残差（`|r| <= delta`）：仍是二次，保持 OLS 的效率；
- 大残差（`|r| > delta`）：变成线性，单个 burst 的“拉力”被限制。

对 `rho_delta` 求导，得到 IRLS 权重：

$$
w(r)=\frac{\psi(r)}{r}=
\begin{cases}
1, & |r| \le \delta \\
\delta/|r|, & |r| > \delta
\end{cases}
$$

所以“robust IRLS”就是：**用当前残差算一组权重 `w`，做一次加权最小二乘，
再用新系数重算残差/权重，迭代几次。**

#### 6.7.2 本接口的迭代步骤

选定 support 后，迭代 `robustIter=2` 次：

```text
1. 用当前权重 wCur 做 support 上的 WLS（lambda=0）
2. 计算 latent prediction: latent = A*beta + b
3. 只取 uncensored 样本的残差:
       rU = y(U) - latent(U)
4. scale = 1.4826 * MAD(rU)
5. delta = 1.345 * scale
6. wH = 1 for all samples
   对 uncensored: wH = min(1, delta ./ max(|rU|, eps))
7. wCur = wBase .* wH
8. 重复
```

其中：

- `1.4826 * MAD` 把“中值绝对偏差”换算成近似标准差；
- `1.345` 是正态分布下 95% 渐近效率的经典 Huber 常数；
- `wBase` 是接口的 `sampleWeight`（接缝/无效样本为 0，正常为 1）。

#### 6.7.3 为什么只对 uncensored 样本降权

censored 样本的观测值是 `y_i = c`。如果 latent prediction
`A_i beta+b > c`，这是**需要被 `rho` 惩罚的 hinge 违反**，不是“y 的离群点”。

如果也对 censored 样本按 `|y - max(A beta+b,c)|` 降权，就会把本来应该被
`rho` 压下去的违反项也削弱，等于削弱 censoring 约束。因此：

- uncensored 样本：按 Huber 权重降权 burst；
- censored 样本：保持 `wBase`（通常为 1），让 `rho` 继续约束
  `A beta+b <= c`；
- `wBase=0` 的接缝/无效样本始终保持 0。

#### 6.7.4 Huber 用在两处，但不在 lambda path

当前代码里 Huber 用于：

1. **初始 OLS**（`robustInitIter=2`）：为 adaptive 权重提供更稳的
   `beta_OLS`；
2. **post-selection refit**（`robustIter=2`）：得到最终 robust 系数。

lambda path 阶段仍然只用 `wBase`（不跑 Huber）。原因是：

- 每个 lambda 都跑 IRLS 会让计算量成倍增加；
- 先用普通加权 LASSO 选 support/lambda，再在固定 support 上做 robust
  refit，已经能把 burst 对最终系数的影响降下来；
- 如果数据 burst 极多、连 support 选择都被带偏，可以再考虑把 robust
  放到 path 阶段（计算代价更高）。

> 历史说明：旧的 `Adaptive_LASSO_StanceSwing.m`（UN.sto 路径）没有
> Huber IRLS，只有 adaptive LASSO + post-selection OLS。上面的
> `robustInitIter` / `robustIter` 是这次为 MocoInverse 实测数据新增的。

### 6.8 输出系数

最终保存的是 `Beta_refit`（robust refit 后的原始尺度系数），并做：

```matlab
beta(abs(beta) <= nzTol) = 0;
mask = abs(beta) > nzTol;
```

这是为了满足 `load_lasso_reflex_controller.m` 的检查：

```text
beta(~mask) 必须严格为 0
```

### 6.9 伪代码

```text
for muscle m = 1..9
  for phase p in {stance, swing}
    idx = phase samples (delay-adjusted)
    wBase = sampleWeight(idx)

    OLS on uncensored samples -> beta_OLS
      (Huber IRLS, robustInitIter=2)
    adaptive w = 1/max(|beta_OLS|, eps)^gamma
    A_star = A ./ w

    for lambda in path
      beta_star = sgl_fit(A_star, y, 0, lambda, weights=wBase)
      beta_ada  = beta_star ./ w
      S = support(beta_ada)
      beta_refit = sgl_fit(A_refit, y, 0, 0, weights=wBase)
      rmse = weighted_rmse(y, max(A beta_refit + b, c), wBase)
    end

    select lambda by (1+tolSelect)*min(rmse), first feasible

    for iter = 1..robustIter
      beta_refit = sgl_fit(A_refit, y, 0, 0, weights=wCur)
      latent = A beta_refit + b
      U = y > c + tol                 % 只对 uncensored 样本算 Huber
      rU = y(U) - latent(U)
      wHuber = ones(size(y))
      wHuber(U) = min(1, 1.345*1.4826*MAD(rU) ./ max(|rU|, eps))
      wCur = wBase .* wHuber          % censored 样本保持 wBase
    end
  end
end
```

---

## 7. 输出文件与字段

### 7.1 `moco_lasso_dataset.sto`

SCONE 风格，列名与 `UN.sto` 兼容。主要列：

```text
<muscle>.activation
<muscle>.excitation
<muscle>.fiber_length_norm
<muscle>.mtu_force_norm
pelvis_tilt / pelvis_tilt_u
hip_flexion_r / _u
knee_extension_r / _u
ankle_dorsiflexion_r / _u
leg1_r.state
leg0_l.state
```

另外会写左腿列，以及可选的 `pelvis_tx/ty`（仅当 `qNames` 里存在）。
`extractStoMuscleFeatures.m` 可以直接读取，不需要修改旧 reader。

### 7.2 `moco_lasso_dataset.mat`

`ds` 结构的主要字段：

| 字段 | 含义 |
|---|---|
| `time` | 100 Hz 时间网格 |
| `muscleNames` | 18 块肌肉（right 9 + left 9） |
| `rawActivation/rawExcitation` | 清洗前信号（可能为空） |
| `activation/excitation` | 清洗后信号（excitation 是 LASSO target） |
| `lCEN/fATN` | 18 块肌肉的归一化长度/力 |
| `q/qdot/qNames` | 坐标与速度 |
| `phase` | stance/state/GRF/events |
| `validMask/sampleWeight` | LASSO 有效样本权重 |
| `qc` | 清洗/特征/相位 QC |
| `meta` | 来源文件、模型、参数、时间戳 |

### 7.3 `results/lasso_controller_result_moco.mat`

`lasso` 是 grouped v2：

```matlab
lasso.format    = 'grouped';
lasso.nPhases   = 5;
lasso.phaseIds  = 0:4;
lasso.nGroups   = 2;
lasso.groups(1).label  = 'stance';
lasso.groups(1).phases = [0 1];
lasso.groups(1).beta   = 26×9;
lasso.groups(1).bias   = 1×9;
lasso.groups(1).mask   = 26×9 logical;
lasso.groups(2).label  = 'swing';
lasso.groups(2).phases = [2 3 4];
...
lasso.featureNames = 26×1 cellstr;
lasso.muscleNames  = 9×1 cellstr;
lasso.source       = 'MocoInverse';
lasso.datasetFile  = ...;
lasso.preprocessing = ds.meta;
lasso.qc           = ds.qc;
```

`load_lasso_reflex_controller.m` 会忽略额外字段，只使用
`groups(g).beta/bias/mask/phases/label`、`nPhases`、`phaseIds`。

### 7.4 `results/lasso_moco_fit_diagnostics.mat`

`fit` 结构：

```text
fit.Beta_refit / Intercept_refit     最终系数
fit.YPred_refit / RMSE_refit         预测与误差
fit.Beta_adaptive / YPred_adaptive   未 robust 的 adaptive LASSO
fit.LambdaPath.<phase>{m}            每个肌肉/相位的 lambda path
fit.S                                每个肌肉/相位的 support
fit.sampleWeight / options           复现信息
fit.featureNames                     26 个特征名
```

lambda path 的 `stats` 已经 compact，不会把每次 FISTA 的 20k 步 objHist
全部写入磁盘。

---

## 8. 运行步骤

### 8.1 准备文件

确认以下文件存在（路径在
[`moco_inverse_to_lasso_dataset.m`](<../Sparse Group LASSO Validation/moco_inverse_to_lasso_dataset.m>)
的 CONFIG 中修改）：

```text
experiments/moco/outputs/moco_inverse_solution_*.sto
experiments/moco/outputs/ik_radians.mot
experiments/moco/outputs/grf_window_*.mot
experiments/moco/outputs/grf_window_*.xml
model/human0918.osim
```

### 8.2 生成数据集

```matlab
cd('Sparse Group LASSO Validation')
moco_inverse_to_lasso_dataset
```

检查控制台：

```text
activation present     : 1
excitation present     : 1
coordinates in solution: 0
LASSO time grid        : 210.0000 .. 240.0000 s, dt = 0.0100 s, n = 3001
kinematics source      : .../ik_radians.mot (filtered IK)
target                 : raw excitation
right stance fraction  : ...   % 由 cal_gait_phase 规则决定
...
```

### 8.3 拟合 LASSO

```matlab
Adaptive_LASSO_MocoInverse
```

输出：

```text
results/lasso_controller_result_moco.mat
results/lasso_moco_fit_diagnostics.mat
Sparse Group LASSO Validation/Adaptive_LASSO_MocoInverse.log
```

### 8.4 接入预测仿真

把 `demo_predictiveForwardSimulation_humanModel.m` /
`prepare_fresh_start.m` 里的：

```matlab
lassoFile = fullfile('results', 'lasso_controller_result.mat');
```

改成：

```matlab
lassoFile = fullfile('results', 'lasso_controller_result_moco.mat');
```

其余 CMA-ES 流程不变。

### 8.5 可选：特征重建自检

```matlab
cd('tests')
test_moco_feature_reconstruction
```

预期输出：

```text
max RMSE over muscles = ... (< 0.01)
Feature reconstruction test PASSED.
```

---

## 9. QC 与验收标准

### 9.1 数据 QC

| 检查 | 期望 |
|---|---|
| `lCEN` 范围 | 约 0.5–1.4；不应大量贴 0.2 或 2.5 |
| `fATN` 范围 | 约 0–0.8；均匀缩放 Fmax 不改变 fATN；小负值可接受（passive force） |
| stance fraction | 约 0.45–0.65（取决于 cal_gait_phase 状态机，不是固定验收值） |
| stance 碎片 | 每个 stance run 不应只有 1–2 个样本 |
| `pct_flagged` | 通常 < 10%；过高说明 despike 太激进 |
| `activationReconRMSE` | 通常 < 0.05–0.1；过大说明低通/目标不自洽 |
| `validMask` | 不应少于总样本的 80%（除非显式排除了很长窗口） |

### 9.2 拟合 QC

| 检查 | 期望 |
|---|---|
| `RMSE_refit` | 每块肌肉/相位应有有限值 |
| `nnz` | 不应全为 0；全 0 说明 target 太接近 `c` 或相位样本太少 |
| `lambda2Sel` | 不应总是等于 `lambda2Max`（说明没有进入有效路径） |
| `mask` | `beta(~mask)` 必须严格 0 |
| `lasso.groups` | 每个 group 的 `beta` 为 26×9、`bias` 为 1×9 |

### 9.3 在 B30 初步结果上的参考数值

这些只是当前数据的量级参考，不是判据。注意：`B30_right_final.sto` 是
activation-only 的初步解，新逻辑会直接报 `MissingExcitation`；这里的数值只用于
说明特征/清洗的量级，不代表该文件还能跑通完整接口。

- `lCEN ≈ 0.58–1.33`（不受 Fmax 影响）；
- `fATN ≈ 0–0.79`（均匀缩放 Fmax 不改变该范围）；
- stance fraction 由严格的 `cal_gait_phase.m` 规则决定，不是固定验收值；
- 6 Hz 清洗后最大逐样本跳变：
  - bifemsh `0.341 → 0.085`
  - soleus `0.273 → 0.112`
  - hamstrings `0.150 → 0.034`

---

## 10. 常见问题排查

| 现象 | 可能原因 | 处理 |
|---|---|---|
| `MissingExcitation` | 解里没有 excitation，或某些肌肉的 excitation 整列缺失 | 重新导出完整 18 块 excitation；不要用部分通道文件 |
| `MissingActivation` | 解里没有 activation，或某些肌肉的 activation 整列缺失 | 重新导出完整 18 块 activation；不会从 excitation 积分补全 |
| `NonFiniteExcitation` / `NonFiniteActivation` | 某些肌肉的信号含 NaN/Inf | 检查 Moco 求解是否失败；不要在接口里填 NaN |
| `noMuscleColumns` | 读错文件（例如 IK/GRF） | 确认 `solutionFile` 是 Moco solution |
| `columnNotFound` / GRF 列找不到 | XML 缺失或列名不是 `_vy` | 显式设 `rightColumn/leftColumn` 或提供 XML |
| `invalidPhase` | GRF/COM 数据不满足 `cal_gait_phase.m` 规则 | 检查 GRF 列/符号、`bodyWeight`、COM 计算模型和坐标单位；不要改 0.02 阈值 |
| 相位看起来不对 | 严格规则对 GRF/COM 很敏感 | 检查 pelvis/calcn COM、GRF 竖直分量、模型与 IK 是否匹配 |
| `lCEN` 大量贴 2.5 | 子步太大 / 单位错 / 积分发散 | 减小 `substepMax`，检查 q/qdot 单位 |
| `lCEN` 大量贴 0.2 | 同上，或 activation 全 0 | 检查 Moco 解是否有效 |
| LASSO 全 0 support | 某相位 target 几乎全在 c / 样本太少 | 检查相位比例、`c`、`sampleWeight` |
| loader 报 `beta outside mask` | 手工修改了 lasso | 调用 `fit_adaptive_lasso_phase_model` 的导出逻辑，或手动清零 |
| `fit` 诊断文件过大 | 保存了完整 FISTA objHist | 已用 `compact_stats`；不要改回完整 stats |
| demo 没读到新控制器 | `lassoFile` 仍指向旧文件 | 改成 `lasso_controller_result_moco.mat` |

---

## 11. 与旧 `UN.sto` 路径的关系

| 项目 | 旧路径 | 新路径 |
|---|---|---|
| 入口 | `Adaptive_LASSO_StanceSwing.m` | `Adaptive_LASSO_MocoInverse.m` |
| 数据 | `Sparse Group LASSO Validation/UN.sto` | `moco_dataset/moco_lasso_dataset.sto` |
| 特征来源 | SCONE 已算好的 `fiber_length_norm/mtu_force_norm` | 接口用 human0918 重建 |
| 相位 | SCONE `leg1_r.state` | GRF 重建 |
| 目标 | SCONE `.excitation` | 清洗后的 raw excitation |
| 权重 | 全 1 | `sampleWeight` + Huber IRLS |
| 输出 | `results/lasso_controller_result.mat` | `results/lasso_controller_result_moco.mat` |
| 旧文件是否修改 | 否 | 否（只新增） |

---

## 12. 设计取舍与局限

1. **只拟合右腿，假设双侧对称**：与旧流程一致；左腿列保留，但当前不参与
   拟合。以后可扩展为左右独立参数。
2. **特征用 human0918 而不是 SQR**：保证与运行时控制器特征空间一致。两者的
   `Fmax` 现已对齐，但 `lopt` / tendon slack / 几何仍有差异，所以 `fATN`
   仍不是完全同源；这是为了保证“训练特征 = 运行时特征”而接受的近似。
3. **相位 0..4 严格遵循 `cal_gait_phase.m`**：不是 60/40、20/60/20 的
   启发式，而是“GRF + pelvis/calcn COM”的五相位状态机。当前 LASSO 的
   stance 组是 state 0/1；phase 2 Liftoff 仍接触地面，但旧 reader 把它归入
   swing 组（state>=2）。
4. **6 Hz 低通是经验值**：对步行合适；跑步/快速反射应做 4/6/10 Hz 敏感性分析。
5. **Huber IRLS 只作用于 post-selection refit**：lambda path 仍用 base weight，
   计算量可控；如果数据 burst 极多，可考虑在 path 阶段也做 robust。
6. **纤维长度积分依赖 OpenSim MATLAB API**：接口的纯 MATLAB reader 只能读表，
   真正算 lCEN/fATN 需要 OpenSim 类可用。
7. **不自动判断 Moco 解质量**：残差/储备作动器是否饱和、activation 是否饱和，
   仍应先用 `experiments/moco/check_solution_quality.py` / `read_moco_results.m`
   检查；接口只负责“把解转换成 LASSO 可用数据”，不负责判断物理解是否可信。

---

## 13. 建议的消融/敏感性实验

1. `doDespike = true` vs `false`（Hampel 的影响）；
2. `robustIter = 0` vs `2`；
3. `lowpassHz = 0 / 4 / 6 / 10`；
4. `featureFilterHz = 0 / 6 / 10`；
5. `delaySteps = 1` vs `2`（10 ms vs 20 ms）；
6. `c = 0.01` 与 `cleanOpts.floor` 的一致性检查；
7. 不同 `lambdaPathConfig.tolSelect`（0.02 / 0.05 / 0.10）；
8. 用 `excludeTimeWindows` 排除/不排除 Moco 分段接缝。

每组实验都建议保存独立的 `lasso_controller_result_moco_<tag>.mat`，避免覆盖。

---

## 14. 模型 Fmax 变更后的兼容性（2026-11）

### 14.1 这次改了什么

`model/human0918.osim` 中 18 块肌肉的 `max_isometric_force` 都被改为原来的
**2 倍**。改完后，human0918 的 Fmax 与 `SQR_simbody.osim` 的对应肌肉一致，
例如：

| 肌肉 | 旧 human0918 Fmax | 新 human0918 / SQR Fmax |
|---|---:|---:|
| hamstrings | 2594 N | 5188 N |
| bifemsh | 804 N | 1608 N |
| glut_max | 1944 N | 3888 N |
| iliopsoas | 2186 N | 4372 N |
| rect_fem | 1169 N | 2338 N |
| vasti | 4530 N | 9060 N |
| gastroc | 2241 N | 4482 N |
| soleus | 3549 N | 7098 N |
| tib_ant | 1579 N | 3158 N |

注意：`optimal_fiber_length`、tendon slack 和肌肉路径**没有**跟着改，因此
human0918 与 SQR 仍然不是完全相同的模型。

### 14.2 关键结论：均匀缩放 Fmax 不改变归一化特征

`fATN` 的定义是：

$$
\mathrm{fATN} = \frac{F_{\text{fiber along tendon}}}{F_{\max}}
$$

Hill-type 肌肉在固定 activation、固定 lCEN/lCEV 时，肌肉力与 `Fmax` 成比例；
分母也按同一比例变化，因此 **fATN 对“所有肌肉统一乘同一个系数”的 Fmax
缩放是不变的**。归一化纤维长度 `lCEN` 由肌肉路径/肌腱平衡决定，也不因均匀
Fmax 缩放而改变。

在 human0918 上做了数值验证：把当前模型在内存中缩放到旧 Fmax（×0.5），
在相同 q / activation 下比较：

```text
a=0.05: lCEN maxdiff=0, fATN maxdiff=0
a=0.20: lCEN maxdiff=0, fATN maxdiff=0
a=0.50: lCEN maxdiff=0, fATN maxdiff=0
a=0.90: lCEN maxdiff=0, fATN maxdiff=0
```

因此：

- MocoInverse → LASSO 接口输出的 `lCEN` / `fATN` **不会**因为这次 Fmax
  加倍而整体减半；
- `UN.sto` 里的 `*.mtu_force_norm` 仍然是有效的归一化力，**不要**把它乘以
  0.5；
- 旧 LASSO 控制器中 fATN 行的系数仍然对应同一个特征尺度，**不要**把
  beta 第 10–18 行乘以 2。

> 只有在“非均匀缩放”（不同肌肉乘不同比例）或同时修改了 `lopt`/tendon
> slack/肌肉路径时，fATN 才会真正改变；那时才需要按肌肉重新拟合/缩放。

### 14.3 真正变化的是什么

Fmax 翻倍改变的是**绝对力容量**，不是归一化特征：

1. **前向动力学**：同一 excitation 下肌肉力约翻倍，模型加速度/运动学会变；
   原来能走好的控制器可能需要重新调参。
2. **静态优化**：同样关节力矩下，所需 activation 降低；`a_opt` 会变。
3. **MocoInverse**：如果 Moco 本身用 human0918 跑，excitation/力轨迹会变；
   当前示例的 Moco 解用 SQR_simbody，**不受 human0918 Fmax 影响**。
4. **CMA-ES 最优解**：搜索空间的特征没变，但适应度地形变了，旧最优参数
   不再是最优。

### 14.4 哪些程序会自动适配

以下代码都是**从 .osim 动态读取** Fmax，不存硬编码数值，因此不需要改代码：

| 程序 | 行为 |
|---|---|
| `read_muscle_static_prop.m` | `fopt = getMaxIsometricForce`，自动读新值 |
| `iterative_static_optimization.m` | 用 `fopt` 构造 QP，自动用新 Fmax |
| `forward_simulation.m` | OpenSim 模型自身的力上界自动用新 Fmax |
| `build_lasso_features_from_kinematics.m` | `fATN = F / getMaxIsometricForce`，自动用新 Fmax |
| `moco_inverse_to_lasso_dataset.m` | 通过上面的特征函数自动得到 fATN |

所以 **MocoInverse → LASSO 接口本身不需要改代码**。

### 14.5 哪些数据和结果需要重做

| 文件/结果 | 是否受 Fmax 特征尺度影响 | 建议 |
|---|---|---|
| `Sparse Group LASSO Validation/UN.sto` | 否（归一化力不变） | 不需要缩放；若想得到新模型下的真实力轨迹，重新跑 SCONE |
| `results/lasso_controller_result.mat` | 否（特征尺度不变） | 可作为 CMA-ES 初值；但闭环动力学变了，应重新优化 |
| `results/lasso_controller_result_moco.mat` | 否（特征尺度不变） | 同上；可重新拟合 LASSO 做一致性检查 |
| `results/checkpoint.mat` / `opt_result_*.mat` | 是（`a_opt`、fitness、CMA 状态对应旧力容量） | 不要直接 resume；从 fresh start 跑 |
| `results/sim_result_*.sto` | 否（lCEN/fATN 归一化特征仍可复现） | 可作为特征管线的回归测试；不要作为新模型动力学结论 |

### 14.6 推荐操作顺序

1. 确认 `model/human0918.osim` 的 Fmax 已经是目标值；
2. 重新跑 `moco_inverse_to_lasso_dataset`（特征数值应基本不变，用于确认
   fATN/lCEN 没有被非预期地改动）；
3. 重新跑 `Adaptive_LASSO_MocoInverse`（可选，用于确认系数稳定）；
4. 旧 UN.sto 路径不需要缩放 `mtu_force_norm`；如果想完全重做，可用新模型
   重新跑 SCONE 后再跑 `Adaptive_LASSO_StanceSwing`；
5. 备份/删除 `results/checkpoint.mat`，从 fresh start 跑 CMA-ES；
6. 用 `tests/test_moco_feature_reconstruction.m` 和 QC 图确认 lCEN/fATN
   量级正常。

### 14.7 如果你观察到 fATN 真的变了

说明不是“所有肌肉统一 ×2”这么简单，或者 `lopt`/tendon slack/几何也被改了。
此时按肌肉换算：

```text
fATN_new = (Fmax_old / Fmax_new) * fATN_old
```

并把旧 LASSO 控制器中对应肌肉的 fATN 行乘以 `(Fmax_new/Fmax_old)` 以保持
同一特征输入下的输出；但更推荐直接用新模型重新生成数据并重新拟合，因为
非均匀 Fmax 改变也会改变最优控制。
