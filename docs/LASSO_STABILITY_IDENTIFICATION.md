# 多随机种子一致性结构辨识（LASSO stability identification）

> 关键词：stability selection / consensus structure identification /
> run-level bootstrap / 0-1 mask / 冻结结构
>
> **本文是“参考手册”版；如果你觉得这里的行话太多，请先读从零讲起的
> [`STABILITY_IDENTIFICATION_EXPLAINED.md`](STABILITY_IDENTIFICATION_EXPLAINED.md)**
> （含 π/Jaccard/RMSE 的通俗解释、手算例子、每一步的输入输出、参数手册、FAQ）。
>
> 相关代码：
> [`bootstrap_phase_runs.m`](<../Sparse Group LASSO Validation/bootstrap_phase_runs.m>)、
> [`consensus_reflex_structure.m`](<../Sparse Group LASSO Validation/consensus_reflex_structure.m>)、
> [`lasso_stability_identification.m`](<../Sparse Group LASSO Validation/lasso_stability_identification.m>)、
> [`Adaptive_LASSO_MocoInverse_Stability.m`](<../Sparse Group LASSO Validation/Adaptive_LASSO_MocoInverse_Stability.m>)

---

## 0. 一句话总览

同一份实验数据、同一套辨识程序，用 **K 个随机种子**各跑一次 adaptive LASSO，
得到 K 个 26×9 的 **0-1 结构（mask）**；统计每个
`(feature, muscle)` 条目被选中的频率 π，用**逐条目多数投票**得到共识结构，
最后在**全量数据**上把 mask 冻结、只重新估计系数。

输出的 controller 里：

- **0-1 图**是共识结构，CMA-ES **不能改**（`unpack_lasso_reflex_params.m`
  每帧强制 `beta(~mask) = 0`，向量长度在 `load` 时就固定）；
- **系数**是全量数据重拟合值，只作为 CMA-ES 的 **initial guess**（real2sim
  gap 由优化吸收）。

---

## 1. 为什么需要它：现在的链路是确定性的

`fit_adaptive_lasso_phase_model.m` 的整条链——

```
extractStoMuscleFeatures → 加权 OLS（Huber IRLS）→ adaptive 权重
→ lambda path（零初值 + 确定性 warm start）→ Huber refit
```

——**没有任何随机数**（`Sparse Group LASSO Validation/` 下检索
`rng|rand|randn|randperm|seed` 零命中）。因此：

> 直接“换 20 个种子跑 20 次”，会得到 20 份**逐位相同**的 mask，
> 13/20 会退化成 20/20，这个方法就失去意义。

真正让结构抖动的是 lambda 选择那条规则：

```
minErr    = min(RMSE_refit)
threshold = (1 + tolSelect) * minErr      % tolSelect = 0.05
selected  = path 上第一个 RMSE_refit <= threshold 的 lambda
```

这是一个**离散跳变**：数据/权重动一点点，选中的 lambda 就换一档，support
整块变化。所以本方法要做的事是：

1. **定义随机性从哪来**（对同一份实验数据做重采样）；
2. 在这个随机性上做**一致性统计**。

这在统计学习里就是 **stability selection**（Meinshausen & Bühlmann, 2010）
的思想，本仓库称之为“多随机种子一致性结构辨识”。

---

## 2. 随机性的定义：run-level bootstrap

### 2.1 为什么按“步态片段”重采样

步态数据的相邻样本高度自相关。如果逐样本 i.i.d. 重采样，等于假设样本独立，
会**系统性低估**结构的不稳定性（看起来“很稳”，实际不是）。

`D.phaseIdxMap.stance` / `.swing` 天然是一段段连续区间（一个步态周期里
stance 一段、swing 一段）。`bootstrap_phase_runs.m` 先把索引切成连续 run，
再**整段**有放回抽样，段内的时间相关性完整保留：

```
[1..60] [101..160] [201..260] ...        % 原始 run
   ↓ 按 seed 抽 nRuns 个 run（有放回）
[1..60] [201..260] [101..160] ...        % 一个 replicate 的样本集合
```

每个种子的样本总数**接近但不等于**原始样本数（等于被抽中 run 的长度之和），
不截断任何 run——这是 block bootstrap 的标准做法。

### 2.2 三种模式

| 模式 | 行为 | 用途 |
|---|---|---|
| `run_bootstrap`（默认） | 整段有放回抽 `nRuns` 段 | 主方法 |
| `run_subsample` | 不放回抽 `ceil(0.8*nRuns)` 段 | 消融：会把“样本变少”混进来，不作默认 |
| `iid_bootstrap` | 逐样本有放回 | 反面对照：预期过于“稳定” |

### 2.3 什么随种子变，什么不变

**只有“抽到哪些样本”随种子变。** 这些全部固定不变：

- 模型、特征定义、`.sto` 读取、delay 裁剪、相位划分；
- `gamma`、`epsilonOLS`、`rho`、`c`、`tolSelect`、lambda 网格、Huber 常数；
- 标准化常数 `muX / sX`（**用全量数据算一次**）。

固定 `muX/sX` 是刻意的：它只影响条件数，不改变预测；固定后不同 replicate 的
系数在**同一尺度**上可比，`lambda2Max` 的量级也稳定，避免把“预处理噪声”
混进“样本噪声”。

`lambda2Max` 则**每个 replicate 重算**，因为它就是 lambda 网格的起点，必须
跟当前样本走。

---

## 3. 算法

```
读取 .sto（一次）  →  X, X_std, Y, w, phaseRuns           [load_lasso_design_matrices]

for k = 1..K                                               [lasso_stability_one_replicate]
    seed_k = baseSeed + k
    idx_k  = run_bootstrap(stance) ∪ run_bootstrap(swing)   ← 唯一随种子变的东西
    拟合 idx_k 上的 9 肌肉 × 2 组 adaptive LASSO
    mask_k(g) = support 拼成的 26×9 逻辑矩阵
    （只保留 mask / RMSE / support size，不保留完整 fit）

Pi(g) = mean_k mask_k(g)                                   [consensus_reflex_structure]
C(g)  = (Pi(g) >= 0.5)                    % 多数投票（默认）

β_final = 全量数据 + 冻结 mask C 的 Huber refit            [fit_adaptive_lasso_core
                                                            mode = 'mask_override']
```

> 记号：设计矩阵直接用 `X`（旧版写作 `A = X*M`，`M` 为单位阵，已移除）。

第 3 步是**一次**独立的全量重拟合，不是 bootstrap 平均。原因：bootstrap
平均会收缩系数、可能把条目平均成 0 破坏结构；全量重拟合用的是真实数据、
真实尺度，才适合当 CMA-ES 的初值。

---

## 4. 共识规则

| rule | 定义 | 说明 |
|---|---|---|
| **`majority`**（默认） | `C = (Pi >= 0.5)` | 逐条目多数投票 = mask 的**坐标中位数**，即到所有 replicate 的 Hamming 距离之和最小的 mask |
| `mode` | 出现次数最多的**那个真实 mask** | 你最初描述的做法；并列时取平均 Jaccard 最高者 |
| `medoid` | 与其余 replicate 平均 Jaccard 最高的那次运行 | “最典型的一次运行” |
| `stability` | `C = (Pi >= piThreshold)`，默认 0.7 | Meinshausen–Bühlmann 风格的严格版 |

### 为什么多数投票正好等于“13/20 就选它”

设 K 个 replicate 中有 m 个产出**完全相同**的 mask M。对任意条目，
M 在该条目上的取值至少被 m 个 replicate 支持，因此

```
Pi(entry) >= m/K   或   Pi(entry) <= 1 - m/K
```

只要 `m/K > 0.5`，多数投票就逐条目还原 M。

**所以“20 个种子里 13 个收敛到同一框架 → 选这个框架”不需要额外逻辑，
它就是多数投票的自然结果。** 这条性质写进了单元测试
（`tests/test_lasso_stability_identification.m` Test 1），是不可回退的契约。

> 只有 `Pi == 0.5` 恰好成立时（K 为偶数）需要表态：
> `tieBreak = 'include'`（默认，保守地保留参数，让 CMA-ES 自己把系数压到 0）
> 或 `'exclude'`（更稀疏）。

---

## 5. 统计量与质量门限

对每个 phase group 单独统计（stance 与 swing 用不同的样本集，必须分开）：

| 量 | 含义 |
|---|---|
| `Pi` (26×9) | 逐条目入选频率，即 0-1 图的“置信度图” |
| `exactMaskCounts` | 不同 mask 的出现次数（降序），即“多少种子共享同一个框架” |
| `modeCount / modeFraction` | 众数 mask 的票数 / 占比 |
| `jaccardMean` | 两两 Jaccard 均值（整体稳定性标量） |
| `jaccardToConsensus` | 每个 seed 与**最终共识 mask** 的 Jaccard |
| `conformingFraction` | Jaccard ≥ `conformJaccard` 的种子占比，即“有多少种子支持选中的框架” |
| `dominantClusterFraction` | Jaccard ≥ `clusterJaccard`(0.9) 建图后最大连通簇占比（严格版诊断） |
| `nStableOn / nStableOff / nAmbiguous` | π≥0.9 / π≤0.1 / 其余 的条目数 |
| `supportSizes` | 各 replicate 的 nnz |

### ⚠️ 三个“门限”不是同一个东西

很容易混淆，这里明确区分：

| 参数 | 是什么 | 默认 | 触发时 |
|---|---|---|---|
| `piThreshold` | **π 阈值**，只被 `stability` 规则使用 | 0.7 | — |
| `conformJaccard` | 判定“某个 seed 与共识 mask 同框架”的 Jaccard 阈值 | 0.8 | — |
| `minStabilityJaccard` | **种子之间两两 Jaccard 的均值下限**（真正的稳定性门限） | 0.35 | 警告“结构可重复性弱” |
| `minDominantFraction` | 仅报告：严格“同一框架”的占比 | 0.6 | 不触发警告 |

多数投票用的 π 阈值**固定为 0.5**，不暴露成参数——否则“13/20”的性质就没了。

**为什么门限不是“完全相同的 mask 占比”？** 234 个 0/1 条目几乎不可能逐位重复：
本数据 K=20 时众数 mask 只出现 1 次（`modeFraction = 0.05`），但共识仍有
123/106 个条目、且冻结后 RMSE 与未冻结全量拟合相差 <1%。因此：

- **稳定性**用 `jaccardMean`（种子两两相似度）衡量；
- **共识的代表性**用 `meanCoverage`（共识保留了多少每个种子的条目）衡量。
  注意共识是**超集**（多数投票会纳入半数以上种子都选的条目），所以它与单个
  种子的 Jaccard 天然偏低，`conformingFraction` 只作诊断，不要当门限；
- **保真度**用 `rmseDelta`（相对未冻结全量拟合）衡量。

`minStabilityJaccard = 0.35` 的含义：如果种子之间的平均 Jaccard 低于 0.35，
说明不同重采样选出的结构几乎不重叠，这份数据的 0-1 图**不可重复**，应当
先检查数据/清洗/`tolSelect`，而不是直接拿去优化。

---

## 6. 输出

| 文件 | 内容 |
|---|---|
| `results/lasso_controller_result_moco_consensus.mat` | `lasso`（grouped v2，可直接给 `load_lasso_reflex_controller.m`）+ `lasso.identification`（K、baseSeed、规则、modeFraction、dominantClusterFraction、maskSignature、totalParameters…）+ `lasso.stability`（`Pi`、各规则 mask、clusterSizes） |
| `results/lasso_stability_diagnostics.mat` | `R`：`R.replicates`（每个种子的 mask/RMSE/supportSize/seed）、`R.consensus`（全部一致性统计）、`R.fitFinal`（全量重拟合）、`R.rmseDelta`、`R.warnings`、`R.timing` |
| `results/lasso_stability_checkpoint.mat` | 断点续跑（`options.resumeFromCheckpoint = true` 时复用已完成的种子） |
| `results/lasso_stability_*.png` | π 热图、mask 频次谱、Jaccard 矩阵、support size 分布、RMSE 对比、π 直方图 |
| `Sparse Group LASSO Validation/Adaptive_LASSO_MocoInverse_Stability.log` | 运行日志 |

`load_lasso_reflex_controller.m` 只读取 `groups(g).beta/bias/mask/phases/label`
与 `nPhases/phaseIds`，**忽略** `identification` / `stability` 等额外字段，
因此完全向后兼容。

---

## 7. 运行

```matlab
cd('Sparse Group LASSO Validation')
Adaptive_LASSO_MocoInverse_Stability
```

前提：先跑过 `moco_inverse_to_lasso_dataset.m`（生成
`moco_dataset/moco_lasso_dataset.sto/.mat`）。

关键参数（脚本顶部 CONFIG）：

| 参数 | 默认 | 说明 |
|---|---|---|
| `nSeeds` | 20 | 种子数 K |
| `baseSeed` | 2026 | 第 k 个种子 = `baseSeed + k` |
| `resampleMode` | `run_bootstrap` | 重采样方式 |
| `consensusRule` | `majority` | 共识规则 |
| `tieBreak` | `include` | π=0.5 时取“有值” |
| `clusterJaccard` | 0.9 | 聚类连边阈值（严格诊断） |
| `conformJaccard` | 0.8 | 判定“某 seed 与共识同框架”的 Jaccard 阈值 |
| `minStabilityJaccard` | 0.35 | 种子两两 Jaccard 均值下限（稳定性警告线） |
| `maxRMSEIncrease` | 0.15 | 冻结 mask 后全量 RMSE 涨幅警告线 |
| `computeReferenceFit` | true | 额外跑一次**未冻结**的全量单次拟合，作为 RMSE 的公平参照（种子 RMSE 是 in-bag 的，偏乐观，不能当基准） |
| `useParallel` | true | 种子之间 `parfor` |
| `workerThreads` | 1 | 每个 worker 的 BLAS 线程数。26×26 的最小二乘不需要多线程；16 个 worker 各自多线程会严重超订，反而更慢 |
| `nRetry` | 3 | 单个种子失败重试（seed + 1000·attempt） |

**成本**：单次完整拟合实测约 **126 s**（16 核机器，`moco_lasso_dataset`，
nLambda=60）。20 个种子在 16 核上是 2 波 `parfor`，约 **4–6 min**；
串行约 42 min。

中断后把 `options.resumeFromCheckpoint = true` 再跑，已完成的种子不会重算。

---

## 8. 与 CMA-ES 的接口

```
lasso.identification.maskSignature      % 0-1 图的指纹（两个 group 各一个整数）
lasso.identification.totalParameters    % = nnz(mask_stance) + 9 + nnz(mask_swing) + 9
```

`prepare_fresh_start.m` 里把 `lassoFile` 指向新的 `.mat` 即可：

```matlab
lassoFile = fullfile('results', 'lasso_controller_result_moco_consensus.mat');
```

运行时：

- `load_lasso_reflex_controller` 把 mask 展开成每相位 cell，并构造扁平向量布局；
- `unpack_lasso_reflex_params` 每次调用都执行 `beta(~mask) = 0`；
- `CMAES_optimization` 只看到固定长度的向量，**没有任何路径能改变 0-1 图**。

因此“0-1 图严格冻结、系数只作初值”在架构上成立；`maskSignature` 的存在
只是为了让“这个 controller 的结构到底是哪一版”可追溯、可比对。

---

## 9. 验收标准

| 检查 | 期望 |
|---|---|
| 有效种子数 | `nValid == K`（若失败会自动重试；仍失败会在日志里列出） |
| 两个 group 的 `nnz(consensus)` | 均 > 0，且与 `median(supportSizes)` 同量级 |
| `jaccardMean` | ≥ 0.35；低于它说明 0-1 图不可重复，**不要直接拿去优化** |
| `meanCoverage` | 共识保留了多少每个种子的条目；超集共识通常 0.85–1.0 |
| `modeFraction` | 仅作诊断：234 个 0/1 条目很难逐位重复，这项通常很低，不代表失败 |
| `nAmbiguous` | 越小越好；它是“结构在抖”的条目数 |
| `rmseDelta` | 冻结 mask 的全量重拟合 vs **未冻结**的全量单次拟合的相对 RMSE 变化；≤ 0.15 |
| `beta(~mask)` | 严格 0（loader 会检查） |
| `totalParameters` | 与新 controller 的 `reflexParamMap.totalLen` 一致（脚本会 assert） |

---

## 10. 局限与注意事项

1. **bootstrap 仍偏乐观**：即使按 run 分块，约 20 段的重采样粒度有限。π 只用于
   **结构决策**，**不要当 p 值**报告。
2. **run 数少**：若某条腿只有个位数个步态周期，簇占比的估计会很粗。可把
   `maxBlockLen` 设成更小的值，把 run 再切细。
3. **共识 mask 可能更稀**：导致全量 RMSE 上升，所以有 `rmseDelta` 门限。
4. **`Pi` 是“入选频率”，不是“正确概率”**：它衡量的是“在等效数据扰动下
   是否反复被选中”，这正是稳定性的定义，但不是统计显著性。
5. **口径要一致**：replicate 的 RMSE 是在**重采样样本**上算的，共识 controller
   的 RMSE 是在**全量样本**上算的，两者只作量级对比，不能混为一谈。
6. **模型 Fmax 变更**（见 `MOCOINVERSE_TO_LASSO_WORKFLOW.md` §14）后要重跑；
   LASSO 与 CMA-ES 的数据/结果都不能跨版本复用。
7. **系数的数值对 BLAS 线程数敏感，结构不敏感。** 在同一份数据、同一份代码上
   实测：多线程与 `-singleCompThread` 相比，mask / support / selected lambda
   **完全相同**，但近共线的设计矩阵会让系数差到 ~1e-4、RMSE 差到 ~1e-9。
   这不是 bug，而是“平坦方向上的系数本来就不可辨识”。它恰好从反面支持了
   本方法的前提：**值得冻结的是 0-1 结构，不是数值**；跨机器比较系数时要
   固定 `maxNumCompThreads`。

---

## 11. 常见问题

| 现象 | 原因 / 处理 |
|---|---|
| 20 个种子的 mask 完全相同（`modeFraction = 1`） | 说明重采样没生效，或数据量远大于结构不确定性——检查 `resampleMode`、`D.phaseIdxMap` 是否真的有多段 run |
| `dominantClusterFraction` 很低（< 0.5） | 结构不可辨识：数据太少/相位划分不稳/`tolSelect` 太宽。先别拿去优化 |
| `rmseDelta` 很大（> 0.15） | 共识 mask 太稀。检查 `nAmbiguous` 高的条目，考虑 `tieBreak='include'` 或更多种子 |
| 某个种子失败 | 日志会给 `#k: <message>`；自动重试 3 次（seed+1000·attempt） |
| `resumeFromCheckpoint` 没有跳过任何种子 | checkpoint 的 `nSeeds`/`baseSeed` 与本次不一致，会被判定为不兼容并忽略（会 warning） |
| 想让结构更保守/更稀疏 | `consensusRule='stability'` + `piThreshold=0.7`，并对比 `mode`/`medoid` 的 nnz |

---

## 12. 与旧路径的关系

| 项目 | 旧路径 | 本次新增 |
|---|---|---|
| 单次拟合入口 | `Adaptive_LASSO_MocoInverse.m` | 不改（输出文件名不变） |
| 拟合函数 | `fit_adaptive_lasso_phase_model.m` | 变为薄 wrapper，签名/默认值/输出**完全不变** |
| 求解器 | `sgl_fit.m` | 只把无条件打印的 summary 改成受 `opts.verbose` 控制（数值不变） |
| 新 controller | — | `lasso_controller_result_moco_consensus.mat` |
| 运行时 | `load_lasso_reflex_controller` → `unpack` → CMA-ES | **完全不动** |

重构等价性由一次性基线对比保证：重构前用旧代码在
`moco_lasso_dataset.sto` 上跑一次并存档，重构后用同一套 options 重跑，
断言 `mask / beta / bias / RMSE_refit_full / S / selectedIdx` 全部逐项相同
（容差 1e-12）。步骤见
[`check_lasso_refactor_equivalence`](<../tests/check_lasso_refactor_equivalence.m>)。

> **实测结果**：同一次会话、同样单线程（`-singleCompThread`）下，重构前后
> 的 stance nnz = 111、swing nnz = 101、total = 230，**所有字段逐位相同**
> （最大 |ΔRMSE| = 0）。
> （该对比是在当时的数据集上做的；结论是“代码重构不改变数值”，与用哪版
> Moco 解无关。）
>
> ⚠️ 对比时两边必须用相同的线程设置。多线程 vs 单线程在本数据上会让系数
> 差 ~1e-4（mask 不变），原因见 §10 第 7 条。

---

## 13. 本数据上的参考结果（K=20, baseSeed=2026）

输入是 MocoInverse 外层交替迭代**收敛后**的解
`experiments/moco/outputs/B30_cont2_final.sto`（与
`moco_inverse_solution_B30_cont2_u021_merged.sto` **MD5 相同**）。
`moco_lasso_dataset.sto`：3001 帧 @100 Hz，有效 2991；stance 1468 帧 / 26 段，
swing 1532 帧 / 27 段。

| 项目 | stance | swing |
|---|---|---|
| 共识 nnz | 123 / 234 | 106 / 234 |
| 稳定的激活条目（π≥0.9） | 45 | 50 |
| 稳定不激活（π≤0.1） | 31 | 62 |
| 模糊条目（0.1<π<0.9） | 158 | 122 |
| 种子两两 Jaccard 均值 | 0.528 | 0.588 |
| 每个种子的 mask 覆盖率（共识保留了多少） | 0.78 | 0.82 |
| 众数 mask 出现次数 | 1 / 20 | 1 / 20 |

- **总参数**：247（`123+9+106+9`）；未冻结的全量单次拟合是 234。
- **保真度**：冻结 mask 的全量重拟合 RMSE 比未冻结全量拟合 **低 0.26%**
  （逐肌肉：+0.23, −0.13, −0.34, −1.35, −3.16, +0.05, −2.30, +1.09, +3.59 %，
  见诊断文件 `rmseDelta`）。
- **耗时**：20 个种子 + 1 次参照拟合，共 **262 s**（16 核 `parfor`，
  `workerThreads=1`）。

**怎么读这组数**：

1. `modeFraction = 1/20`：没有任何两个种子给出**逐位相同**的 234 维 mask。
   这是正常的——逐位相等是个极苛刻的判据；这不代表方法失败。
2. `jaccardMean ≈ 0.53–0.59`：不同重采样选出的结构有大约一半以上重叠，
   属于“中等稳定”。没有触发 `minStabilityJaccard = 0.35` 的警告。
3. `stable-on` 45/50 + `ambiguous` 158/122：真正**稳定激活**的条目只有一小半，
   其余在种子之间摇摆。这些摇摆条目就是 π 图里灰色的一带。
4. 冻结后 RMSE 反而略降：因为多数投票把“半数以上种子都选的条目”都纳入了，
   共识 mask 是**超集**，所以拟合不会变差；代价是 CMA-ES 的维度从 234 涨到 247。

**与旧目标（非收敛解）的对比**：同一流程跑
`_history_files/outputs/moco_inverse_solution_B30_right_final.sto` 得到
121/99 个条目、Jaccard 0.532/0.610、238 个参数、RMSE −0.53%。
换成收敛后的目标后这些数字基本在同一量级（Jaccard 略降、参数略增），
说明**结构不确定性主要来自 LASSO 的离散 lambda 选择，而不是 Moco 解是否收敛**。

> 换句话说：**这份数据的 0-1 图不是唯一可辨识的，但它的“稳定核心 + 多数
> 条目”足以构成一个不劣于单次拟合的控制器。** π 图应当和控制器一起被保存，
> 后续若要收紧结构（更少参数），可以把 `consensusRule` 改成 `'stability'`
> 并用 `piThreshold = 0.7`（本数据 → 79+78+18 = 175 个参数），再比较 nnz 与 RMSE。

---

## 14. 稀疏度：247 个参数是不是太多？

不是“错”，但**约三成是冗余**。三个容易误会的地方先澄清：

1. **冻结不会让 RMSE 变低**——冻结=少用参数，训练误差只会升。之前看到的
   “冻结后低 0.26%”是因为 π≥0.5 的共识是**超集**（比单次拟合多 13 个格子）。
2. **UN.sto 的 22/17 不可比**：那条路径的 `minRefitRMSE` 是 `1.2e-06` 量级，
   数据本身由同类线性反射律生成；Moco 的 excitation 是最优控制解，26 个线性
   特征的最优 RMSE 只能到 0.019。前者是“恢复生成系数”，后者是“近似富信号”。
3. **247 ≠ 全开**：它是每肌肉-相位平均 13 项（26 项里的一半）。

实测的三条曲线（工具：[`analyze_lasso_sparsity.m`](<../Sparse Group LASSO Validation/analyze_lasso_sparsity.m>)）：

| 选择方式 | 参数 | 训练 RMSE | 留出 RMSE |
|---|---|---|---|
| tolSelect=0.05 单次拟合 | 234 | 0.01937 | — |
| π≥0.5（当前 majority 共识） | 247 | 0.01934 | 0.02076 |
| π≥0.6 | 206 | 0.02031 | 0.02156（+4.4%） |
| **π≥0.7（推荐）** | **175** | 0.02147 | 0.02240（+8.5%） |
| π≥0.8 | 144 | 0.02537 | 0.02532（+22.6%） |
| tolSelect=0.2 单次拟合 | 155 | 0.02203 | — |
| tolSelect=0.3 单次拟合 | 117 | 0.02330 | — |

关键证据：**过拟合 gap（留出 − 同折训练）从 π≥0.3 的 +20.9% 单调降到 π≥1.0 的
+2.0%**；π≥0.5 时仍有 +14.5%，说明稠密 mask 有相当一部分在拟合噪声。留出误差
最低点在 π≥0.4，π≥0.5 只差 0.5%，砍到 π≥0.7 也只贵 8.5%。（“训练/留出在
π≈0.8 处交叉”是早前的误读，已按同折比较修正。）

**推荐**：要“简洁且每项都有 ≥70% 种子支持”，用
`consensusRule='stability'; piThreshold=0.70; resumeFromCheckpoint=true`
（复用已算好的 20 个种子，约 2 min），得到
`results/lasso_controller_result_moco_consensus_pi70.mat`（175 个参数）。
更狠的 `tolSelect=0.2`（155 参数）需要重跑一次 ensemble 才自洽。

另外，**i.i.d. 重采样会夸大稳定性**（Jaccard 0.67/0.75 vs run-level 0.53/0.59，
模糊条目几乎减半），但不改变 nnz——所以稳定性只能按 run-level 报告。
工具：[`compare_bootstrap_resampling.m`](<../Sparse Group LASSO Validation/compare_bootstrap_resampling.m>)。
完整推导与表格见
[`STABILITY_IDENTIFICATION_EXPLAINED.md`](STABILITY_IDENTIFICATION_EXPLAINED.md) §7.5–7.6。
