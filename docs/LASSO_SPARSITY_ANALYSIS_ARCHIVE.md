# LASSO 稀疏度分析归档（论文素材）

> 目的：回答“这个反射框架要多少参数才够”，并给出可写进论文的证据链。
>
> 相关工具：[`analyze_lasso_sparsity.m`](<../Sparse Group LASSO Validation/analyze_lasso_sparsity.m>)、
> [`plot_lasso_fit_quality.m`](<../Sparse Group LASSO Validation/plot_lasso_fit_quality.m>)、
> [`compare_bootstrap_resampling.m`](<../Sparse Group LASSO Validation/compare_bootstrap_resampling.m>)
>
> 方法细节见 [`STABILITY_IDENTIFICATION_EXPLAINED.md`](STABILITY_IDENTIFICATION_EXPLAINED.md)；
> 符号与 `.mat` 格式见 [`LASSO_STABILITY_IDENTIFICATION.md`](LASSO_STABILITY_IDENTIFICATION.md)。

---

## 0. 一句话结论

在本数据上，**247 个参数的 majority 共识里有约三成冗余**；把共识阈值提到 π≥0.7
可降到 **175 个参数**（−29%），代价是留出 RMSE +8.5%、视觉上几乎不可分辨。
若按“精度-稀疏度”曲线的膝点（2–3 项/肌肉）选，可以压到 ~90 个参数，但 RMSE
要约 1.5 倍——那已经接近“线性模型解释不了的部分”，不适合作为唯一方案。
**参数不能无限压缩，这本身可能反映了下肢 Ia 异源投射的解剖学事实**（见 §5）。

---

## 1. 研究问题与约束

目标**不是**精确拟合 MocoInverse 的 excitation，而是得到一个**尽量简洁、可解释**
的肌肉反射框架（0-1 结构），作为 CMA-ES 的结构先验。两条硬约束：

1. **MocoInverse 的 excitation 必然含噪声**（最优控制解 + sEMG 驱动的权衡），
   拟合精度超过某个点就没有生理意义；
2. **参数越多，方法越失去意义**：muscle-reflex-based predictive simulation 的
   立论是“用少量反射参数编码运动”；参数过多就退化为回归问题。

因此评价标准不是“RMSE 越小越好”，而是**在可接受的拟合下参数最少**。

---

## 2. 数据与评估协议

| 项 | 内容 |
|---|---|
| 输入 | `moco_lasso_dataset.sto`（3001 帧 @100 Hz，210–240 s），来自收敛解 `B30_cont2_final.sto` |
| 目标 | 清洗后的 excitation（Hampel + 6 Hz 低通），floor = 0.01 |
| 特征 | 26 个：9 `lCEN` + 9 `fATN` + 8 骨盆/关节角度与角速度 |
| 结构 | 2 组（stance = 相位 0/1，swing = 相位 2/3/4），每组 26×9 |
| 辨识 | censored adaptive LASSO（λ path 60 点 + post-selection refit + Huber IRLS） |
| 稳定性 | 20 个种子 × run-level bootstrap；π = 入选频率 |
| 参数计数 | `nnz(stance) + 9 + nnz(swing) + 9`（9 个偏置永远存在） |

> **记号**：设计矩阵在论文与
> [`REFLEX_NEUROPHYSIOLOGY_DESIGN_AND_VALIDATION.md`](REFLEX_NEUROPHYSIOLOGY_DESIGN_AND_VALIDATION.md)
> 中记作 `X`，在代码/旧文档中记作 `A`（`fit.A`；旧式 `A = X*M`，`M` 为单位阵，已移除）。
> 模型即 `y = Xβ + b`，每个相位组一套 `β, b`。

**评估口径（2 折，按步态段）**：

- `train(full)`：全部帧拟合、全部帧打分（偏乐观）；
- `train(fold)`：一折拟合、**同一折**打分（与留出同规模，用于算过拟合）；
- `held-out`：一折拟合、**另一折**打分（泛化）；
- 过拟合 gap = `(held-out − train(fold)) / train(fold)`；
- `相对最优留出` = 相对所有 π 阈值中最小的 held-out 的百分比差。

> 为什么按整段分：相邻帧高度相关，随机分帧会信息泄漏。详见 EXPLAINED §7.5。

---

## 3. 结果

### 3.1 λ path 的边际收益（每块肌肉、每个相位组平均）

| 容忍误差（相对该肌肉-相位最优 RMSE） | stance 项数 | swing 项数 |
|---|---|---|
| +2% | 14.1 | 12.6 |
| **+5%（当前 tolSelect）** | **12.7** | **11.3** |
| +10% | 10.1 | 10.1 |
| +20% | 7.7 | 7.6 |
| +50% | 4.0 | 4.6 |
| +100% | 2.7 | 3.3 |

**Kneedle 膝点（逐 muscle-group，排除 `support=0`）**：对 18 条路径各自的
`(support, 相对 RMSE)` 曲线，在非支配子集上取膝点：

- stance：nnz = `[6 5 2 5 2 2 3 3 3]`，合计 **31**；
- swing ：nnz = `[3 2 2 3 5 8 2 3 7]`，合计 **35**；
- 加 18 个偏置 = **84 个参数**；膝点处的 RMSE 为各自最优的 **1.14–2.13 倍**。

**为什么排除 `support=0`**：λ 最大时 β 全零，预测只剩 9 个偏置，这不是一个可用的
反射结构；它的 RMSE 很大（>3×），会把 Kneedle 的弦拉向原点，使膝点过度偏左。
实测：不排除时合计 63 参数（如 `vasti` 的 swing 膝点只有 1 项），排除后 84 参数。

**选中点在哪里**：规则保证每个 muscle-group 的选中点**紧贴 `y = 1+tolSelect`
的下方**——实测选中点 1.000–1.049、其前一个 λ 为 1.052–1.072（全部跨过 1.05 线）。

→ 也就是说：**“再加一项”的边际价值在 2–3 项/肌肉之后迅速衰减**；当前 5% 容差
（12–13 项）落在曲线的平段末端，确实保守。但这条曲线**只看训练精度、不看稳定性**：
按膝点选择稀疏化得到 84 参数，平均 RMSE 比未冻结参考高 **+40.4%**（逐肌肉 1.14–2.13 倍），
仍明显比“推荐”的 175 参数激进。

### 3.2 π 阈值扫描（冻结 mask、全量重拟合；2 折 CV）

| π | 参数 | stance | swing | train(full) | train(fold) | held-out | 过拟合 gap | 相对最优留出 |
|---|---|---|---|---|---|---|---|---|
| 0.30 | 325 | 173 | 134 | 0.01873 | 0.01729 | 0.02091 | +20.9% | +1.3% |
| 0.40 | 283 | 147 | 118 | 0.01891 | 0.01762 | **0.02065** | +17.2% | 0 |
| 0.50 | 247 | 123 | 106 | 0.01934 | 0.01813 | 0.02076 | +14.5% | +0.5% |
| 0.60 | 206 | 97 | 91 | 0.02031 | 0.01925 | 0.02156 | +12.0% | +4.4% |
| **0.70** | **175** | 79 | 78 | 0.02147 | 0.02042 | 0.02240 | +9.7% | **+8.5%** |
| 0.80 | 144 | 59 | 67 | 0.02537 | 0.02396 | 0.02532 | +5.7% | +22.6% |
| 0.85 | 132 | 53 | 61 | 0.02597 | 0.02448 | 0.02576 | +5.2% | +24.7% |
| 0.90 | 113 | 45 | 50 | 0.02752 | 0.02612 | 0.02714 | +3.9% | +31.4% |
| 0.95 | 101 | 38 | 45 | 0.02961 | 0.02812 | 0.02916 | +3.7% | +41.2% |
| 1.00 | 79 | 26 | 35 | 0.03267 | 0.03101 | 0.03163 | +2.0% | +53.2% |

**Kneedle 膝点**：**175 个参数（π≥0.70）**，held-out 0.02240。

三个信号：
1. 过拟合 gap 单调下降（+20.9% → +2.0%）：稠密 mask 有相当一部分在拟合噪声；
2. held-out 最低在 π≥0.4，但 π≥0.5 只差 0.5%——当前 247 不算错，只是**冗余**；
3. π≥0.8 之后 held-out 跳涨（+22.6%），开始删真信号。

### 3.3 tolSelect 扫描（全量单次拟合）

| tolSelect | 参数 | stance | swing | train RMSE | 相对 0.05 |
|---|---|---|---|---|---|
| 0.05 | 234 | 114 | 102 | 0.01937 | — |
| 0.10 | 200 | 91 | 91 | 0.02013 | +3.9% |
| 0.20 | 155 | 69 | 68 | 0.02203 | +13.7% |
| 0.30 | 117 | 49 | 50 | 0.02330 | +20.3% |

### 3.4 膝点汇总

| 方案 | 参数 | 怎么实现 | 备注 |
|---|---|---|---|
| λ path，**逐 muscle-group 膝点**（排除 nnz=0） | **84** | `selectionMode='knee'; minSupport=1`（已实现） | 31+35 项；RMSE 1.14–2.13× |
| **膝点规则 + π 共识（v2，K=20）** | **58 / 63 / 72（π=0.7/0.6/0.5）** | `selectionMode='knee'` + `consensusRule='stability'` | Jaccard 0.42/0.46；RMSE +74–84% → **过于激进** |
| π 阈值（总参数 vs 留出 RMSE） | **175** | `consensusRule='stability'; piThreshold=0.70` | **推荐**：留出 +8.5%，每项 ≥70% 种子支持 |

> **注**：`tolSelect` 是逐肌肉的停止规则，与“取曲率最大点”不等价，所以两者给出的
> 参数数不同。若最终采用膝点选择，`tolSelect` 可以完全不用；π 仍作为事后的一致性
> 统计（每个 replicate 各自取膝点 → 20 个 mask → 多数投票）。该方案需要先实现
> `selectionMode='knee'` 并重跑一次 ensemble（约 4–5 min）。

### 3.4b 膝点规则 + π 共识（v2 实验，K=20, π≥0.7, minSupport=1）

按用户要求正式运行了一次“膝点规则”的 ensemble：

- 每个 replicate 各自对 18 条路径取 Kneedle 膝点（排除 nnz=0），得到 20 个 mask；
- **replicate 间一致性明显低于 tolSelect 规则**：平均 Jaccard
  stance **0.415** / swing **0.463**（tolSelect 是 0.528 / 0.588），
  modal 1/20（0.05），coverage 0.48 / 0.56；
- π≥0.7 共识：stance **19** + swing **21** = 40 项 + 18 偏置 = **58 参数**；
- 全量数据膝点拟合（不冻结）是 87 参数；共识相对它 RMSE **+26.5%**，
  相对 234 参数参考 **+83.7%**（平均 RMSE 0.03558 vs 0.01937）；
- 逐肌肉 RMSE：`[0.0167 0.0684 0.0300 0.0331 0.0257 0.0132 0.0255 0.0781 0.0295]`
  （soleus、bifemsh 最差）。

**结论**：膝点规则确实把参数压到最少（58），但代价是精度大幅下降、稳定性也更差。
原因是**膝点是路径上的一个“拐点”，而拐点附近曲线很陡**：重采样让它左右滑动 1–2 个
λ，support 就换一批，π 因此被稀释；tolSelect 的 5% 容差带更宽，选中的 λ 在
重采样下更稳定。若要坚持膝点规则，建议 `minSupport≥2–3` 或把 π 降到 0.5。

**π=0.6 的对照（同一批 replicate，只换阈值）**：共识变为 stance 21 + swing 24 =
**63 参数**，平均 RMSE **0.03512**（相对 87 参数全量膝点拟合 +25.1%，相对 234 参考
**+81.3%**）。与 π=0.7 相比只多 5 个参数、RMSE 改善 1.3%——说明**瓶颈不是 π 阈值，
而是膝点本身的不稳定**：π 落在 [0.6, 0.7) 的条目很少，放宽阈值也捞不回精度。

**π=0.5 的对照（多数边界）**：共识变为 stance 28 + swing 26 = **72 参数**，平均 RMSE
**0.03372**（相对 234 参考 **+74.1%**，相对全量膝点拟合 +18.9%）。整条趋势：
π = 0.7/0.6/0.5 → **58 / 63 / 72 参数**，RMSE +83.7% / +81.3% / +74.1%。
改善存在但离 175 参数版本的 +11% 仍很远；π 低于 0.5 就不再是“多数支持”，
所以 **0.5 是这条路线的实用上限**。

### 3.5 拟合效果（图）

- `results/lasso_fit_quality.png`：9 块右腿肌肉的全程 30 s（247 vs 175 两条预测叠加）；
- `results/lasso_fit_quality_zoom.png`：3 s 放大，单个步态周期；
- `results/lasso_fit_rmse_comparison.png`：逐肌肉 RMSE 柱状对比（234 / 247 / 175）；
- `results/lasso_fit_quality_ctrl.png` / `_zoom_ctrl.png` / `_rmse_comparison_ctrl.png`：
  由 `plot_controller_fit_quality.m` 从**任意 controller `.mat`** 直接生成
  （此处是 84 参数的 knee controller，见 §3.4）。

肉眼结论：**175 参数的预测与 247 几乎重合**，主要差异在少数 burst 的峰值
（如 `tib_ant_r`、`gastrocnemius`），这也解释了 +8.5% 的 RMSE 差。

---

## 4. 推荐配置

| 方案 | 参数 | 生成方式 | 留出代价 | 适用 |
|---|---|---|---|---|
| A. majority 共识 | 247 | `consensusRule='majority'` | +0.5% | 保守、想尽量保住信号 |
| B. stability π≥0.7 | 175 | `consensusRule='stability'; piThreshold=0.7; resumeFromCheckpoint=true` | +8.5% | 需要“每项 ≥70% 支持”的稳定性叙事时 |
| **D. 膝点规则，无稳定性筛选（当前采用）** | **84** | `selectionMode='knee'; minSupport=1`（全量单次拟合） | **0.0303 held-out** | **CMA-ES 直接输入 + 初值**；参数最少且泛化最好 |
| C. 膝点 + π（消融） | 58 / 63 / 72 | `consensusRule='stability'` + π=0.7/0.6/0.5 | 0.0418 / 0.0403 | 仅作“π 会过度删减”的对照 |

产物：

- **D（当前）**：`results/lasso_controller_result_moco_consensus_knee.mat`（84 参数；
  `prepare_fresh_start.m` 的 `lassoFile` 已指向它）
- B：`results/lasso_controller_result_moco_consensus_pi70.mat`
- C：`results/lasso_controller_result_moco_kneerule{,_pi60,_pi50}.mat`

> 选择建议：**当前以 D 作为 CMA-ES 的输入**（跳过 stability selection；out-of-fold
> 0.0303，与 nested 膝点选择 0.0302 一致）。π 只保留为“每条连接可重复性”的诊断量，
> 不再用来裁剪 mask——因为实测 π 投票会删掉“有用但不稳定”的可替代特征
> （bifemsh 的膝角、soleus 的踝角），使 held-out 从 0.030 恶化到 0.040。
> B 作为“稳定性叙事”的备选；C 作为消融。

---

## 5. 解释一：为什么这个 mask 可能注定是稠密的（事后解释，不是先验）

> **定位**：与 [`REFLEX_NEUROPHYSIOLOGY_DESIGN_AND_VALIDATION.md`](REFLEX_NEUROPHYSIOLOGY_DESIGN_AND_VALIDATION.md)
> 一致——神经生理学文献用于**设计动机**与**结果的事后解释**，**不作为优化先验**
> 去剪枝或约束符号；结构由 LASSO 稳定性选择决定，增益由 CMA-ES 优化。

如果进一步压缩参数会明显损失拟合，这不一定是方法失败——**下肢 Ia 的异源
投射本来就是广泛分布的**：

> Meunier S, Pierrot-Deseilligny E, Simonetta M. *Pattern of monosynaptic
> heteronymous Ia connections in the human lower limb.*
> Exp Brain Res (1993) 96:534–544.

该文用 PSTH 方法在 10 名健康受试者上系统测了股神经、胫神经、腓神经、比目鱼肌
神经、腓肠肌神经对 Sol / GM / Per brev / TA / Q / Bi / ST 运动神经元池的单突触
异源 Ia 投射（其 Table 1）。可直接引用的定量事实：

- **股神经（膝伸肌/髋屈肌来源）的 Ia 传入投射到所有被探查的小腿运动神经元池**
  （Sol 10/12 = 83%，TA 55/63 = 79%）——即**跨关节**投射；
- **所有被探查的小腿肌肉都投射到 ST**；小腿三头肌（Sol、GM）还投射到 **Bi 和 Q**；
- **同一刺激可同时兴奋跨关节的拮抗肌**（小腿三头肌 → Q 与腘绳肌；VL → TA 与 Sol）；
- 作者的结论句：**“transjoint monosynaptic connections were found in all explored
  combinations, except pretibial flexors and peroneal muscles to Q”**；
- 但同时：“These monosynaptic projections are widely distributed, **but generally
  quantitatively weak**（异源峰一般 < 同源峰的 20%）”，例外是 GM→Bi、Q→Sol 较强；
- 还有“**近协同肌之间反而可能缺失**”的反例：GM→Sol 为 0/9，Sol→GM 仅 12%；
- 机制上，作者提出这些“预接线”的通路可以通过**下行控制突触前抑制**在不同任务
  间被“聚焦/选择”（其 Discussion 明确提到 presynaptic inhibition 对 Ia 传入的
  门控）。

**推论（事后解释）**：如果 β 的非零项在生理上对应这些单突触/寡突触通路，那么
一个 9×26 的反射矩阵**本来就不该是稀疏的**——非零项多、但多数增益小。这与我们
观察到的“参数数量降不下去、且多数条目增益较小”的图景在**定性上一致**。注意这只是
**解释**，不是“因此必须稠密”的约束。

**必须写清楚的边界**（避免审稿人误读）：

1. LASSO 的 mask 是**统计/功能**意义上的支持集，不等于解剖学连接矩阵：它混合了
   同源 lCEN/fATN（肌梭、GTO）与关节角度项，还受任务、噪声、特征共线性影响；
2. “稠密”不等于“每一项都真实”：π 分析显示相当一部分条目只有 40–60% 的
   可重复性，稳定性选择正是在剔除它们；
3. 1993 年文章测的是**电刺激下的单突触 EPSP**，与行走中功能性的反射增益不是
   同一层面（后者受突触前抑制、运动单位募集等调制）。

因此论文里应表述为：“**解剖学证据给出一个广泛连接、单个增益较小的预期；我们的
辨识结果与这一解剖学预期在定性上一致，且稳定性分析指出其中哪一部分是可重复的
核心。**”（这是**事后解释**，不是“用先验约束优化”。）

### 5.1 证据强度：哪些能用来说、哪些不能

| 结论 | 证据强度 | 本文用法 |
|---|---|---|
| Ia 异名投射**广泛、跨关节、但多数弱**（Meunier 1993） | 中（解剖实测；功能增益未知） | 解释“为什么稠密”；**不**用于设定非零 |
| `GM→Sol` 单突触 Ia 兴奋缺失（0/9，多研究重复） | 中–高 | 仅作事后核对 |
| Ib 从非单突触；安静时**自体抑制**（负力反馈） | 较硬（经典回路） | 解释“负力反馈”项合理 |
| 行走时**伸肌 Ib 反转**为正力反馈 | 较硬（猫）/中（人间接） | 解释 stance 组 `fATN` 正项 |
| 交互抑制：双突触 Ia、甘氨酸能、**跖屈→背屈更强** | 中–高 | 解释拮抗肌对负项与不对称 |
| 反射增益的**相位/任务调制**（站立 > 步行 > 跑步） | 高（人体 H-reflex 实测） | 支撑 stance/swing 两套 `β` |
| “Ib 贡献 SOL 支撑相 EMG ~50%” | **不成立**（Sinkjaer 是传入总贡献） | **不写** |
| 跨关节踝→膝 Ib | **无证据** | 若出现只作探索性讨论 |

> 更完整的对照清单（Ia 长度反馈 / Ib 力反馈 / 相位）见神经生理学文档 §7.3；
> 表述边界（7 条）见其 §7.4 与
> [`REFLEX_CONTROLLER_INTERPRETATION.md`](REFLEX_CONTROLLER_INTERPRETATION.md) §9。

### 5.2 本文如何使用这些证据（三步）

1. **设计动机**（Introduction/Model）：长度反馈 ← Ia；力反馈 ← Ib；负项 ←
   交互抑制/自体抑制；相位调制 ← H-reflex 的相位依赖。这套动机与 Geyer & Herr
   (2010) 的做法一致（力学/反射原理驱动设计，再用数据验证），而**不是**从
   神经生理参数赋值；
2. **数据驱动辨识**（Methods）：LASSO 稳定性选择定结构、CMA-ES 定增益；
3. **事后核对**（Discussion）：把辨识出的 `β` 与 §5.1 的定性预期对照——
   一致则加强生理合理性，不一致则讨论“新发现 / 模型局限 / 文献证据弱”，
   不强行附会。

**明确不做的**：不用文献设定 `mask`、不施加符号约束、不用文献值当初值。

---

## 6. 解释二：y = Xβ + b 的生理含义

见 [`REFLEX_CONTROLLER_INTERPRETATION.md`](REFLEX_CONTROLLER_INTERPRETATION.md)
（β = 感觉反馈通道；b = 反馈参考 + 脊髓 CPG；二者互补；反射增益本身受相位/任务
调制，正好对应 CMA-ES 的优化对象）。

---

## 7. 局限

1. **留出评估的 mask 来自全量 ensemble**：偏乐观，只能横向比较阈值；
2. 只有 26–27 个步态段/组，2 折 CV 的方差较大；
3. 两个相位组共享同一套 26 特征，未区分同源/异源、单/寡突触；
4. `tolSelect` 与 π 是两个独立的稀疏度旋钮，本文没有穷举其组合；
5. 所有结论基于单一受试者/单次试验（B30 窗口 210–240 s），跨受试者需重跑。

---

## 8. 复现

```matlab
cd('Sparse Group LASSO Validation')

% 1) 数据集（若尚未生成）
moco_inverse_to_lasso_dataset

% 2) 20 个种子的结构辨识 + majority 共识（约 4.5 min）
Adaptive_LASSO_MocoInverse_Stability

% 3) 稀疏度/膝点分析（约 7 min；输出 .mat + 图）
analyze_lasso_sparsity

% 4) π>=0.7 的稀疏 controller（复用 checkpoint，约 2 min）
%    consensusRule='stability'; piThreshold=0.7; resumeFromCheckpoint=true

% 5) 拟合效果图（从 diagnostics 文件）
plot_lasso_fit_quality('results/lasso_stability_diagnostics.mat', ...
    'Sparse Group LASSO Validation/moco_dataset/moco_lasso_dataset.mat', 'results');

% 5b) 从任意 controller .mat 直接画（不需要 diagnostics）
plot_controller_fit_quality( ...
    'results/lasso_controller_result_moco_consensus_knee.mat', ...
    'Sparse Group LASSO Validation/moco_dataset/moco_lasso_dataset.sto', ...
    'results', [], '_ctrl');

% 6) run-level vs iid 对照（约 5 min）
compare_bootstrap_resampling
```

产物清单：

| 文件 | 内容 |
|---|---|
| `results/lasso_stability_diagnostics.mat` | 20 个种子的 mask、π、Jaccard、参考拟合 |
| `results/lasso_controller_result_moco_consensus.mat` | 247 参数 majority controller |
| `results/lasso_controller_result_moco_consensus_pi70.mat` | 175 参数 π≥0.7 controller |
| `results/lasso_sparsity_analysis.mat` | 三条曲线 + 膝点 |
| `results/lasso_sparsity_tradeoff.png` | 4×3 布局：9 个肌肉子图（stance/swing 各一条路径 + 各自膝点）+ π 扫描 + tolSelect 扫描 + 图例 |
| `results/lasso_fit_quality*.png` | 9 肌肉拟合曲线（全程 / 3 s 放大 / RMSE 柱） |
| `results/lasso_bootstrap_mode_comparison.mat` | run vs iid 稳定性对照 |

---

## 9. 论文写法建议

- 表 1（正文或附录）：λ path 边际收益 + 膝点 → 说明 5% 容差保守；
- 表 2（正文）：π 扫描（参数、train、held-out、gap）→ 主结果与推荐；
- 图：`lasso_sparsity_tradeoff.png`（4×3：9 个肌肉子图 + π 扫描 + tolSelect 扫描）；
- 图：`lasso_fit_quality_zoom.png` 展示 175 vs 247 视觉上不可分辨；
- 讨论段：把“参数无法进一步压缩”接到 Meunier 1993 的广泛异源投射；
- 局限段：§7。
