# MocoInverse + sEMG tracking 论文可行性评估

> 评估对象：MocoInverse + sEMG tracking + 逐肌 scale `s` + 外层交替优化；1.5 s 与 30 s 右腿数据；A-form/no-gain 基线；reserve 与 s 弱可辨识性。
>
> 检索截止：2026-11。文献结论以已核验的 DOI/官方文档为准；主观概率与期刊判断单独标注。

---

## 0. 一句话结论

**以“新 pipeline / 首次在 MocoInverse 做 EMG tracking / 逐肌 scale”为主张：不够创新。**  
**但 B-form 不是无差别重参数化：A-form 在增益消掉后有 \(||u||^2\sin^2\theta\) 的幅值偏置；B-form 是 \(||e||^2\sin^2\theta\)，不含 \(||u||^2\)。B-form 确实去掉了 A 的“把 u 往下压”的 tracking 偏置；它仍然留下尺度不可辨识（由 effort 决定幅值），但这是另一件事。**

原因不是数据少，而是核心已经被公开工作占位：

1. **OpenSim Moco 官方示例 `exampleEMGTracking`** 已经在 `MocoInverse` 上用 `MocoControlTrackingGoal` 跟踪 EMG，并用 `addScaleFactor` 给每块肌肉加尺度参数（2022 年 API 文档已存在）：  
   https://opensim-org.github.io/opensim-moco-site/docs/1.2.0/html_user/exampleEMGTracking_answers_8m-example.html
2. **Miller et al. 2025, PLoS Comput Biol** 在 EMG-informed Moco 步行中做逐肌 scale，在 calibration 阶段求解、跨多个 stance 平均后共享，使用右腿、reserve，并有植入式膝关节假体验证：  
   https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012833
3. `(e - s u)^2 = s^2 (u - e/s)^2` 是恒等式；但**标准 `addScaleFactor` 的 A-form 是固定权重、参考乘 g**，增益消掉后 tracking 残差 ∝ \(||u||^2\sin^2\theta\)，会把 u 的幅值往下压；B-form 的增益消掉后是 \(||e||^2\sin^2\theta\)，与 \(||u||\) 无关。  
   所以 B-form 不是“新残差”，但它是 **amplitude-unbiased tracking**；它仍有一个尺度规范方向（s 与 u 的分配），需要 effort/归一化来定。

**目前最可能站得住、且最值得做的核心是：**

> **A-form 的 gain 消元后存在 tracking 幅值偏置；B-form 去掉了该偏置，但 s 仍只通过 \(s\cdot u\) 被识别。把“amplitude-unbiased tracking + s 的可辨识性/归一化”写成主结果。**

把它做成“理论 + A/B/CEINMS 对照 + 可辨识性诊断”，可以形成一篇 **methods / negative-result 技术札记**。这比继续包装 pipeline 更有希望。

---

## 1. 先修正事实：哪些已经存在

### 1.1 官方示例已经做了 MocoInverse + EMG tracking + 逐肌 scale

`exampleEMGTracking_answers.m` 的关键代码：

```matlab
inverse = MocoInverse();
study = inverse.initialize();
tracking = MocoControlTrackingGoal('emg_tracking');
tracking.setReference(TableProcessor(emgReference));
tracking.addScaleFactor('gastroc_factor', '/forceset/gasmed_l', [0.01 1.0]);
tracking.addScaleFactor('tibant_factor',  '/forceset/tibant_l', [0.01 1.0]);
tracking.addScaleFactor('bifem_factor',   '/forceset/bfsh_l',  [0.01 1.0]);
tracking.addScaleFactor('gluteus_factor', '/forceset/glmax2_l',[0.01 1.0]);
solver.resetProblem(problem);
solver.set_parameters_require_initsystem(false);
```

这已经覆盖：

- MocoInverse；
- `MocoControlTrackingGoal` 跟踪 EMG；
- 每块肌肉一个 scale factor；
- `addScaleFactor` 的 MocoParameter；
- `parameters_require_initsystem(false)` 加速。

因此不能写“首次在 MocoInverse 中做 EMG tracking / 逐肌 scale”。

### 1.2 Miller 2025 更接近你的完整设定

Miller et al. 2025（PLoS Comput Biol 21(9):e1012833）做 EMG-informed Moco 步行：

- prescribed kinematics + GRF（MocoInverse 式）；
- 代价包含 excitation effort、膝/踝 reserve、其他 reserve、EMG tracking；
- calibration 阶段对每块肌肉引入一个 scale factor（约 [0.5, 1]），与 excitation/reserve 在**同一个优化里联合求解**；
- 在无负重条件的 10 个 stance 上各求一组，跨 stride 取平均，得到每肌肉一个 scale，再固定用于 execution；
- 只用右腿 9 路 EMG（另 3 路映射到相邻肌肉），膝/踝 reserve 重罚，EMG 有约 53 ms 平移；
- 有植入式膝关节假体验证。

这几乎覆盖了“逐肌 s + 跨段共享 + 右腿 + reserve + Moco”的组合。  
你的工作仍有差异（B-form、**外层交替而不是单层联合**、30 s 长窗、并行分段、可辨识性诊断），但**不能把组合本身当主创新**。Miller 的 s 是单层联合估计，不是你的外层交替；这是唯一还可守的“估计层级”差异之一。

### 1.3 CEINMS/EMG-driven 的标定对象不同，但概念相近

Pizzolato et al. 2015（CEINMS）标定的是：

- 最优纤维长度、肌腱松弛长度（±5%）；
- 全局 activation dynamics shape factor、C1/C2；
- 按功能分组的 strength coefficient（约 0.5–2.5）。

也就是说，CEINMS 用**机械参数**吸收“EMG 幅值与力矩不一致”，不直接用 excitation 的乘性 scale。  
你的 `s` 与 strength coefficient 在力矩通道上共线，但在 EMG 通道上不同；这是可以写的差异，但必须证明 `s` 不是 Fmax 的重新参数化。

Sartori et al. 2014 的 hybrid/EMG-assisted 目标已经是：

\[
\alpha(\bar M-M)^2+\beta e^2+\gamma(\bar e-e)^2
\]

即“固定增益 1 的 A-form + 最小努力正则”。所以 A-form 不是新思想；**但标准 `addScaleFactor` 的 A-form 与 B-form 不是同一个优化问题**，这正是需要写清楚的地方。

### 1.4 数学事实：B-form 去掉幅值偏置，但不让 s 可辨识

对单块肌肉：

\[
J_B=\int(e-su)^2dt=s^2\int(u-e/s)^2dt.
\]

这个恒等式成立，但要注意：**它对应的是 A-form 取 \(g=1/s\) 并且每肌肉权重乘 \(s^2\)**。  
标准 `MocoControlTrackingGoal::addScaleFactor` 的做法是“参考乘 g、权重固定”，目标近似：

\[
J_A^{\text{standard}}=\int(u-g e)^2dt.
\]

标准 A-form 消掉增益 \(g\) 后：

\[
g^*=\frac{u\cdot e}{\|e\|^2},\qquad
J_A^{\min}=\|u\|^2\sin^2\theta.
\]

残差依赖 \(\|u\|^2\)：即使形状（夹角 \(\theta\)）不变，缩小 u 也能降低 tracking；这就是 A 的幅值下压偏置，g 会倾向 lower bound。

B-form 消掉 \(s\) 后：

\[
s^*=\frac{e\cdot u}{\|u\|^2},\qquad
J_B^{\min}=\|e\|^2\sin^2\theta.
\]

残差**不含 \(\|u\|^2\)**，所以 B 的 tracking 只惩罚形状，不在“整体放大/缩小 u”的方向上出力。  
C-form（\(1-\cos\theta\)）也一样不含增益，梯度与 u 本身正交。

统计上可以这样理解：A 是把 \(u\) 回归到含噪的 \(e\) 上，存在经典 **attenuation/regression dilution**：e 的噪声会把斜率 g 压小；B 是把含噪的 \(e\) 回归到模型 \(u\) 上，是 errors-in-variables 的“正确方向”。  
因此 B 的 amplitude-unbiased 不只是代数技巧，而是与“EMG 是含噪观测、u 是模型量”的统计结构一致。  
C 更彻底：完全丢弃幅值，只保留形状；如果下游只需要 EMG 形状，C 是更干净的选择；B 保留 s 作为幅值 readout，但这 s 不是可辨识的生理量。

**B 仍然有一个尺度规范方向**：如果 u 恰好平行 e，取 \(u=e/s\)，任意 s 都使 tracking=0。  
但 **B 不会因此额外奖励负荷转移**：s 消元后的 tracking 是 \(\|e\|^2\sin^2\theta\)，不含 \(\|u\|\)；把 u 压小并不会让 B 的 tracking 更好。  
u 也不能无限小，因为 MocoInverse 的力矩/动力学是硬约束，其他肌肉或 reserves 必须补位。  
所以 B 的 s 只是**被动读取**当前的 u 幅值/负荷分配；如果 effort 或动力学把负荷分到别处，s 会跟着变，但 B 没有额外激励去制造这种转移。A 才有 \(\|u\|^2\) 的 tracking 下压项，会主动奖励“把 tracked 肌肉的 u 压低”。

因此更准确的表述是：

> **A-form 的 gain 消元后存在 tracking 幅值偏置；B-form 去掉了该偏置，但 s 仍只通过 \(s\cdot u\) 被识别。B-form 是 amplitude-unbiased，但不是 scale-identifiable。**
> 要真正让 s 有生理意义，还需要归一化/外部锚（固定几何平均、MVC/力量锚、尺度不变损失或先验）；但这是另一层问题，不否定 B 对 A 偏置的修正。

A/B 的两个闭式尺度估计也不是互为倒数：

\[
\hat s_B=\frac{\int eu}{\int u^2},\qquad
\hat g_A=\frac{\int eu}{\int e^2},
\]

含噪时 \(\hat s_B\hat g_A\approx R^2\)（经典 regression dilution）。  
这个关系本身可以作为一个方法学结果，而不是“两种基线差不多/不一样”。

---

## 2. 你的工作逐条拆解

| 你的主张 | 文献状态 | 可怎么写 |
|---|---|---|
| MocoInverse + EMG tracking | 官方示例已有 | 只能作为 methods 实现，不能当 novelty |
| 每肌肉一个 scale | 官方示例 `addScaleFactor`；Miller 2025 | 只能作为 setting，不能当 first |
| 跨段共享 s | Miller 2025 跨 stance 平均 | 需要证明你的跨段估计/方差处理有额外价值 |
| B-form `(e-su)^2` | 不是新残差；A-form 下界退化，B-form 上界退化，二者对偶 | **对偶退化 + 实际表现是可辩护 delta**；需 A/B/CEINMS head-to-head |
| 外层闭式 s 更新 | 交替最小二乘/坐标下降 | 数学上是标准工具，靠诊断取胜 |
| 1.5 s vs 30 s、分段并行 | 未检索到系统研究 | 可作为第二贡献，但需偏差-方差/收敛分析 |
| A-form/no-gain 基线 | A-form 是 EMG-assisted 谱系；no-gain 是 s=1 特例；标准 addScaleFactor 会压 lower bound | 必须复现 official addScaleFactor 并展示退化，不能只用 no-gain |
| s 弱可辨识 | 现象已知；本 pipeline 的量化是开放 | **最有潜力的核心** |
| reserve 与 s 共线 | 文献有参数补偿，但 reserve 作为诊断量未系统做 | 可做成机制实验 |

---

## 3. 创新性判定

### 3.1 作为“方法/流水线”论文：核心不够，但 B-form 的 amplitude-unbiased 是可辩护 delta

- 官方示例 + Miller 2025 会直接被审稿人指出；
- “交替优化 + 闭式最小二乘”是标准 alternating least squares/coordinate descent；
- 当前只有 1 名受试者、右腿、30 s、无独立肌肉力验证；reserve 仍大；
- A/B 的 objective 数值不可直接比较，且形状指标 hamstrings/gastroc 反而变差。

但 **A/B 的 tracking 偏置不是一回事**：A-form 消掉 g 后残差 \(\|u\|^2\sin^2\theta\)，会压 u 的幅值；B-form 消掉 s 后是 \(\|e\|^2\sin^2\theta\)，与 \(\|u\|\) 无关。  
所以 B-form 确实去掉了 A 的 tracking 幅值偏置；它仍然有尺度不可辨识（s 与 u 的分配由 effort 决定），但这是另一层问题，不是 A 的下压偏置。  
真正要证明的是：B 在实际数据上是否比 A 的 `addScaleFactor` 产生更合理的 u/力矩/reserve，而不只是“不会撞下界”。

### 3.2 作为“amplitude-unbiased tracking + s 可辨识性”论文：够一个技术札记的雏形

可辩护的核心：

1. **A-form 幅值偏置**：标准 `addScaleFactor` 的固定权重 A-form 消掉 g 后残差 \(\propto\|u\|^2\sin^2\theta\)，会把 u 往下压；解释为什么 scale 会趋向 lower bound。
2. **B-form 无幅值偏置**：B-form 消掉 s 后残差 \(\propto\|e\|^2\sin^2\theta\)，不含 \(\|u\|\)；所以 tracking 只约束形状。
3. **B 仍有尺度不可辨识，但不会额外 reward hacking**：若 u 平行 e，任意 s 都能让 tracking=0；但 B 的 tracking 不含 \(\|u\|\)，压小 u 不会降低 B 的 tracking，因此没有“为了 B 而转移负荷”的额外激励。s 只是读取 effort/动力学已经决定的 u 幅值；A 才有 \(\|u\|^2\) 的下压项。  
   当前 30 s 检查：左腿肌肉 RMS 变化 <0.1%，右腿未跟踪肌肉约 ±3%，右踝 reserve RMS +20%（绝对值小）；尚未看到大规模负荷转移，但必须查肌肉力/力矩贡献，不能只看 u。
4. **可检验预言**：跟踪越好、激励能量越低、reserve 越强，s 的 profile 越平、CI 越宽、越容易漂移。
5. **诊断**：profile likelihood、multi-start、effort 权重扫描、reserve 强度扫描、u 饱和/边界统计。
6. **归一化方案**：固定几何平均 s=1、外部 MVC/力量锚、尺度不变损失、或显式先验；比较能否把 s 从 bounds 拉回内部。
7. **报告建议**：主结果报告 \(s\cdot u\) 与拟合残差；s 只附 CI 和条件（权重、bounds、窗长），不报“生理增益”。

这可以成为：

> **Amplitude-unbiased EMG tracking with weakly identifiable scale: B-form vs A-form in inverse optimal control**

但必须补解析、合成 ground truth、A/B/CEINMS head-to-head、profile、multi-start；仅凭当前 6 次更新不够。

### 3.3 作为“算法比较”论文：有潜力，但需重做

如果你能证明：

- 标准 A-form `addScaleFactor` 在 [0.01,1] 下确实压 lower bound、跟踪退化；
- B-form 在相同数据/目标/weights/bounds 下去掉 A 的幅值偏置；s 仍不可辨识，但它是 u 幅值的 readout，不应把 load transfer 归因于 B；比较 tracking/力矩/reserve 的实际表现；
- 原生 fixed-s 外层交替 vs 官方 `addScaleFactor` 联合优化在速度、收敛、参数可辨识性、边界行为上更好；

那么可以是 numerical methods 论文（CMBBE/IJNMBE/JBE）。  
但当前没有这个对比；而且当前 s 弱辨识，不能证明交替更好。

---

## 4. Idea 可以落在哪里

### Idea A（推荐）：amplitude-unbiased B-form + s 可辨识性

**核心问题**：标准 A-form `addScaleFactor` 的 gain 消元后是否把 u 的幅值往下压？B-form 是否去掉该偏置？去掉之后，per-muscle s 还能被数据识别到什么程度？如何用归一化/外部锚让 s 有生理意义？

**贡献**：

- 解析：A 消 g 后 \(\|u\|^2\sin^2\theta\)；B 消 s 后 \(\|e\|^2\sin^2\theta\)；A 有幅值偏置，B 没有；
- 合成 ground truth：已知 s_true，测试 A 的幅值下压、B 的纯形状跟踪；
- 真实数据 head-to-head：A-form `addScaleFactor` [0.01,1] vs B-form native；报告 scale 轨迹、下界命中、tracking、u、s·u；
- profile likelihood：每块肌肉的 ΔJ(s) 与 95% CI；
- multi-start：不同 s0 是否到同一 s；
- effort 权重/reserve 强度扫描：s 是否系统移动；
- 归一化方案：固定几何平均 s=1、外部 MVC/力量锚、尺度不变损失；比较能否把 s 从 bounds 拉回内部；
- 与 CEINMS/EMG-assisted 对比；
- 报告建议：报 s·u/CI，不报单点生理 s。

**卖点**：不是“又一个 EMG pipeline”，而是“A-form 的 gain 消元会压 u 的幅值；B-form 是 amplitude-unbiased tracking；但 s 仍弱可辨识，需要归一化/报告 s·u”。

**适合**：J Biomech 短篇/technical note；JEK；Ann Biomed Eng；CMBBE；PLoS ONE 保底。

### Idea B：B-form vs A-form `addScaleFactor` 的 head-to-head

**核心问题**：在 Moco 中，逐肌 scale 应该作为 A-form 的参考乘数联合优化，还是 B-form 的模型侧尺度/外层交替？

**贡献**：

- 相同数据/目标/权重/bounds；
- 目标值、scale、u、s·u、耗时、收敛、边界、初值依赖；
- 证明 A-form 是否压 lower bound；B-form 是否去掉幅值偏置；B 的 s 只是 readout，不应把 load transfer 归因于 B；
- 比较两者在实际 tracking/力矩/reserve 上的表现；
- 证明交替是否更稳/更快，或只是近似。

**风险**：如果 B-form 没有实际优势，就不是论文；当前只有初步现象，需正式复现。

**适合**：CMBBE、IJNMBE、J Biomech Eng。

### Idea C：窗长/分段/并行的偏差-方差

**核心问题**：1.5 s 瞬态 vs 30 s 稳态，s 估计和 reserve 使用如何变化？

**贡献**：

- 窗长、overlap、分段数对 s、s·u、力矩残差、reserve 的 bias-variance；
- 并行分段的 seam/收敛；
- 给出长 trial 的分段建议。

**风险**：偏工程；需要多数据集。

**适合**：CMBBE、JBE、IJNMBE。

### Idea D：临床应用（当前不建议）

例如“EMG-informed MocoInverse 估计膝载荷/康复”。需要多受试者、临床问题、外部验证；当前单受试者右腿不够。

---

## 5. 最小可辩护论文 vs 强论文

### 最小（8–12 页 technical note）

- 1 名受试者右腿步行（明确外推限制）；
- 合成 ground truth（已知 s_true，展示不可辨识/恢复）；
- 解析 + profile likelihood（5–8 块代表性肌肉，含低激励/近边界）；
- multi-start（≥5 个 s0）；
- effort 权重/reserve 强度扫描；
- 与官方 `addScaleFactor` 联合优化同数据对比；
- A/B 估计量关系与 R² 衰减；
- 主结论：报 s·u 与残差，s 附 CI，不报生理值。

**可投**：J Biomech 短篇 → JEK → CMBBE → PLoS ONE。

### 强（full paper）

- 5–10 名受试者，双腿，2 个任务/速度；
- held-out trial / leave-one-trial-out；
- 重测信度；
- CEINMS/EMG-assisted 同数据对比；
- 独立肌肉力/接触力验证（植入假体、测力台、超声）；
- 分段/并行收敛性与段数无关；
- 下游结论敏感性：报 s 还是 s·u 会不会改变临床/力学结论；
- 开源代码/数据。

**可投**：J Biomech、Ann Biomed Eng、TBME（验证足够时）、JEK。

---

## 6. 期刊判断

主观概率，不是计量证据。

| 期刊 | 契合 | 门槛 | 当前判断 |
|---|---|---|---|
| Journal of Biomechanics | 高 | 中高 | **最佳目标**，但必须改成可辨识性短篇；现状 <15% |
| J Electromyogr Kinesiol | 高（EMG 增益/归一化） | 中 | 强替代；补最小实验后 50–60% |
| Ann Biomed Eng | 高（EMG-informed 先例） | 中高 | 合适第二选择 |
| CMBBE | 中高（算法） | 中 | 现实；数学贡献偏薄 |
| J Biomech Eng (ASME) | 中高 | 中 | 现实备选 |
| IJNMBE | 中（数值/收敛） | 中 | 备选 |
| IEEE TBME | 中 | 高 | 当前不建议；单受试者几乎必拒 |
| JNER | 中 | 中高 | 只有加康复场景才考虑 |
| PLOS Comput Biol | 中低 | 很高 | 不建议；Moco 主文与 Miller 2025 已在此 |
| PLoS ONE | 中 | 低 | 最稳保底；负结果可发 |
| Scientific Reports | 中 | 低 | 现实保底 |
| Front Bioeng Biotechnol | 高（methods） | 中低 | 现实可投；APC/声誉较低 |
| Gait & Posture | 中 | 中 | 当前不现实（单受试者+纯方法） |
| J Sports Sciences | 低 | 中 | 不建议 |
| BMMB | 低中 | 高 | 易判 out of scope |

**排序**：J Biomech ≈ JEK > Ann Biomed Eng ≈ CMBBE > JBE/IJNMBE > PLoS ONE ≈ Front Bioeng ≈ Sci Rep > TBME > JNER > Gait & Posture > PLOS CB.

---

## 7. 建议的 90 天路线

### 第 1–2 周：停止盲跑，做解析和合成实验

- 写清 B-form/A-form 估计量关系；
- 合成 e/u，已知 s_true，测试 \(s^*=C/U\)、尺度平坦方向；
- 加 profile likelihood 原型；
- 计算每块肌肉的 \(||u||\)、激励能量、s 的渐近 SE。

### 第 3–4 周：多起点 + 权重/reserve 扫描

- 5 个起点 × 1 受试者；
- effort 权重 ×0.5/1/2；
- reserve 强度扫描（加强 reserve 惩罚，s 离散度应下降或不可行）；
- 报告 s 与 s·u 的分离。

### 第 5–6 周：与官方 joint scale 对比

- 同数据、同目标、同 bounds；
- `addScaleFactor` 联合优化 vs 外层交替；
- 目标、s、u、s·u、耗时、收敛、边界；
- 如果交替无优势，论文重心完全转到可辨识性。

### 第 7–10 周：补数据/验证（按目标期刊）

- 最小：1 受试者 + 合成 + profile + 多起点；
- 强：3–10 受试者、双腿、2 任务、held-out；
- CEINMS/EMG-assisted 对照；
- 如做临床/载荷，补力/接触力验证。

### 第 11–12 周：写作

- 标题围绕 identifiability，不围绕 MocoInverse pipeline；
- 主图 4 张：尺度平坦性、profile CI vs ||u||、reserve 扫描、joint vs alternating；
- 主结论：s·u 可识别，s 仅条件可报。

---

## 8. 审稿人最可能的问题

1. **“这不就是 Moco 官方 exampleEMGTracking？”**  
   必须在引言明确引用官方示例，并把贡献转成可辨识性/诊断。
2. **“B-form 只是 A-form 重参数化。”**  
   主动写出恒等式，不要声称新误差；把 A/B 的估计量差异（R² 衰减）作为结果。
3. **“s 与 Fmax/strength coefficient 共线。”**  
   解释 s 作用在 excitation 通道，不是力矩通道；但必须承认在纯力矩意义上不可分。
4. **“excitation ≠ activation，EMG 到 excitation 还有 activation dynamics/EMD。”**  
   报告滤波、EMD 敏感性；或明确 s 吸收这些因素，是 nuisance scale。
5. **“reserve 吸收了力矩误差，s 才有空间漂移。”**  
   做 reserve 强度扫描；报告 reserve 激活/峰值；必要时做无 reserve 对照。
6. **“只有 1 名受试者/右腿。”**  
   技术札记可接受但必须明确外推限制；full paper 需多受试者/双腿。
7. **“A/B 目标值不可比。”**  
   不要用 J 直接证明生理优越；报告 s·u、残差、CI。
8. **“没有肌肉力验证。”**  
   若主张生理，需植入假体/测力台/超声；否则只主张 EMG-informed trajectory。
9. **“s 的 0.5–2 是生理范围吗？”**  
   明确这是算法/先验约束，不是生理边界；做 bounds 敏感性。
10. **“为什么不用 CEINMS/EMG-assisted 对照？”**  
    这是事实标准；至少讨论，强论文补同数据对比。

---

## 9. 当前数据能支撑什么、不能支撑什么

### 能支撑

- 一个可复现的 MocoInverse + EMG tracking + 逐肌 s 的实现；
- 1.5 s/30 s 右腿的初步结果；
- A/no-gain 基线的 tracking 指标对比；
- s 漂移、raw u 变化、s·u 稳定、reserve 大的现象；
- 一个明确的“弱可辨识”信号。

### 不能支撑

- “首次”任何 pipeline；
- s 是生理增益或可跨任务迁移；
- 肌肉力已经验证；
- B-form 幅值无偏 ⇒ s 一定可辨识；
- 交替优化优于官方 joint scale；
- 临床/载荷结论。

---

## 10. 参考文献（核心）

1. Dembia CL, Bianco NA, Falisse A, Hicks JL, Delp SL. OpenSim Moco: Musculoskeletal optimal control. *PLoS Comput Biol.* 2020;16(12):e1008493. https://doi.org/10.1371/journal.pcbi.1008493
2. OpenSim Moco official example: `exampleEMGTracking_answers.m`. https://opensim-org.github.io/opensim-moco-site/docs/1.2.0/html_user/exampleEMGTracking_answers_8m-example.html
3. Miller DE, Brown AE, Bianco NA, Bhise R, Delp SL, Collins SH. How peak knee loads are affected by changing the mass of lower-limb body segments during walking. *PLoS Comput Biol.* 2025;21(9):e1012833. https://doi.org/10.1371/journal.pcbi.1012833
4. Pizzolato C, Lloyd DG, Sartori M, Ceseracciu E, Besier TF, Fregly BJ, Reggiani M. CEINMS: A toolbox to investigate the influence of different neural control solutions on the prediction of muscle excitation and joint moments during dynamic motor tasks. *J Biomech.* 2015;48(14):3929–3936. https://doi.org/10.1016/j.jbiomech.2015.09.021
5. Sartori M, Farina D, Lloyd DG. Hybrid neuromusculoskeletal modeling to best track joint moments using a balance between muscle excitations derived from electromyograms and optimization. *J Biomech.* 2014;47(15):3613–3621. https://doi.org/10.1016/j.jbiomech.2014.10.009
6. Lloyd DG, Besier TF. An EMG-driven musculoskeletal model to estimate muscle forces and knee joint moments in vivo. *J Biomech.* 2003;36(6):765–776. https://doi.org/10.1016/S0021-9290(03)00010-1
7. Buchanan TS, Lloyd DG, Manal K, Besier TF. Neuromusculoskeletal modeling: estimation of muscle forces and joint moments and movements from measurements of neural command. *J Appl Biomech.* 2004;20(4):367–395. https://doi.org/10.1123/jab.20.4.367
8. Hambly MJ, Worsey MTO, Lloyd DG, Pizzolato C. Improving calibration of EMG-informed neuromusculoskeletal models through differentiable physics and muscle synergies. *IEEE Trans Biomed Eng.* 2025;72(12):3453–3463. https://doi.org/10.1109/TBME.2025.3569682
9. Burden A. How should we normalize electromyograms obtained from healthy participants? What we have learned from over 25 years of research. *J Electromyogr Kinesiol.* 2010;20(6):1023–1035. https://doi.org/10.1016/j.jelekin.2010.07.004
10. Ackland DC, Lin YC, Pandy MG. Sensitivity of model predictions of muscle function to changes in moment arms and muscle-tendon properties: a Monte-Carlo analysis. *J Biomech.* 2012;45(8):1463–1471. https://doi.org/10.1016/j.jbiomech.2012.02.023
11. De Groote F, Kinney AL, Rao AV, Fregly BJ. Evaluation of direct collocation optimal control problem formulations for solving the muscle redundancy problem. *Ann Biomed Eng.* 2016;44(10):2922–2936. https://doi.org/10.1007/s10439-016-1591-9
12. Hicks JL, Uchida TK, Seth A, Rajagopal A, Delp SL. Is my model good enough? Best practices for verification and validation of musculoskeletal models and simulations of movement. *J Biomech Eng.* 2015;137(2):020905. https://doi.org/10.1115/1.4029304
13. Raue A, Kreutz C, Maiwald T, et al. Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood. *Bioinformatics.* 2009;25(15):1923–1929. https://doi.org/10.1093/bioinformatics/btp358
14. Gutenkunst RN, Waterfall JJ, Casey FP, et al. Universally sloppy parameter sensitivities in systems biology models. *PLoS Comput Biol.* 2007;3(10):e189. https://doi.org/10.1371/journal.pcbi.0030189
15. Reinbolt JA, Schutte JF, Fregly BJ, Koh BI, Haftka RT, George AD, Mitchell KH. Determination of patient-specific multi-joint kinematic models through two-level optimization. *J Biomech.* 2005;38(3):621–626. https://doi.org/10.1016/j.jbiomech.2004.03.031

---

## 11. 最终建议

**不要把当前工作写成“新的 MocoInverse + sEMG pipeline”。**  
那会被官方示例和 Miller 2025 直接压住。

**把它写成：**

> “在 EMG-informed inverse optimal control 中，per-muscle EMG scale 只通过 s·u 进入目标，因此 s 弱可辨识；我们用 profile likelihood、multi-start、reserve 扫描和 joint-vs-alternating 对照量化这一点，并给出报告规范。”

这条路线：

- 有真实数据 + 合成 ground truth；
- 有明确机制（尺度规范、effort/reserve 打破）；
- 有可检验预言；
- 有临床/方法学警示价值；
- 不依赖“首次”声明。

**最小可投**：J Biomech 短篇 / JEK / CMBBE / PLoS ONE。  
**强论文**：多受试者、双腿、CEINMS 对照、力/接触力验证后，再冲 J Biomech/Ann Biomed Eng。
