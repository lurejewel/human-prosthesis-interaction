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

**Kneedle 膝点**（对 9 块肌肉的中位曲线）：stance **2.0 项**（RMSE/RMSE_best
= 1.67，隐含 tolSelect ≈ 0.67），swing **3.0 项**（1.50，≈ 0.50）。

→ 也就是说：**“再加一项”的边际价值在 2–3 项之后迅速衰减**。当前 5% 容差
（12–13 项）落在曲线的平段末端，确实保守。

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

| 曲线 | 膝点 | 对应该选择什么 | 备注 |
|---|---|---|---|
| λ path（每肌肉精度 vs 项数） | 2–3 项 | tolSelect ≈ 0.5–0.67 | **只看精度**，不考虑稳定性；过于激进 |
| π 阈值（总参数 vs 留出 RMSE） | 175 参数 | π ≥ 0.70 | **兼顾精度与稳定性**，推荐 |

两条膝点不同不矛盾：path 曲线问“再加一项值不值”，π 曲线问“总参数-泛化
的性价比”，后者才是论文里要的“最简框架”。

### 3.5 拟合效果（图）

- `results/lasso_fit_quality.png`：9 块右腿肌肉的全程 30 s（247 vs 175 两条预测叠加）；
- `results/lasso_fit_quality_zoom.png`：3 s 放大，单个步态周期；
- `results/lasso_fit_rmse_comparison.png`：逐肌肉 RMSE 柱状对比（234 / 247 / 175）。

肉眼结论：**175 参数的预测与 247 几乎重合**，主要差异在少数 burst 的峰值
（如 `tib_ant_r`、`gastrocnemius`），这也解释了 +8.5% 的 RMSE 差。

---

## 4. 推荐配置

| 方案 | 参数 | 生成方式 | 留出代价 | 适用 |
|---|---|---|---|---|
| A. majority 共识 | 247 | `consensusRule='majority'` | +0.5% | 保守、想尽量保住信号 |
| **B. stability π≥0.7（推荐）** | **175** | `consensusRule='stability'; piThreshold=0.7; resumeFromCheckpoint=true` | **+8.5%** | 论文主结果：简洁 + 每项 ≥70% 支持 |
| C. tolSelect=0.2 | 155 | 重跑 ensemble（`lambdaPathConfig.tolSelect=0.2`） | 训练 +13.7%（留出待测） | 想更狠，但需重算 π 才自洽 |

产物：B 已生成 `results/lasso_controller_result_moco_consensus_pi70.mat`。

> 选择建议：**以 B 作为主结果**，A 作为“上界/敏感性”，C 作为“极限稀疏”消融。
> 这样论文可以同时讲“冗余有多少”和“最简能到多少”。

---

## 5. 解释一：为什么这个 mask 可能注定是稠密的

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

**推论**：如果 β 的非零项在生理上对应这些单突触/寡突触通路，那么一个 9×26 的
反射矩阵**本来就不该是稀疏的**——非零项多、但多数增益小。这正好对应我们观察到的
“参数数量降不下去、且大多数条目增益较小”的图景。

**必须写清楚的边界**（避免审稿人误读）：

1. LASSO 的 mask 是**统计/功能**意义上的支持集，不等于解剖学连接矩阵：它混合了
   同源 lCEN/fATN（肌梭、GTO）与关节角度项，还受任务、噪声、特征共线性影响；
2. “稠密”不等于“每一项都真实”：π 分析显示相当一部分条目只有 40–60% 的
   可重复性，稳定性选择正是在剔除它们；
3. 1993 年文章测的是**电刺激下的单突触 EPSP**，与行走中功能性的反射增益不是
   同一层面（后者受突触前抑制、运动单位募集等调制）。

因此论文里应表述为：“**解剖学证据预测一个广泛连接、单个增益较小的矩阵；我们的
辨识结果与这一先验一致，且稳定性分析指出其中哪一部分是可重复的核心。**”

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

% 5) 拟合效果图
plot_lasso_fit_quality('results/lasso_stability_diagnostics.mat', ...
    'Sparse Group LASSO Validation/moco_dataset/moco_lasso_dataset.mat', 'results');

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
| `results/lasso_sparsity_tradeoff.png` | 三面板图（含膝点星标） |
| `results/lasso_fit_quality*.png` | 9 肌肉拟合曲线（全程 / 3 s 放大 / RMSE 柱） |
| `results/lasso_bootstrap_mode_comparison.mat` | run vs iid 稳定性对照 |

---

## 9. 论文写法建议

- 表 1（正文或附录）：λ path 边际收益 + 膝点 → 说明 5% 容差保守；
- 表 2（正文）：π 扫描（参数、train、held-out、gap）→ 主结果与推荐；
- 图：`lasso_sparsity_tradeoff.png` 三面板（膝点星标直接落在 π=175 上）；
- 图：`lasso_fit_quality_zoom.png` 展示 175 vs 247 视觉上不可分辨；
- 讨论段：把“参数无法进一步压缩”接到 Meunier 1993 的广泛异源投射；
- 局限段：§7。
