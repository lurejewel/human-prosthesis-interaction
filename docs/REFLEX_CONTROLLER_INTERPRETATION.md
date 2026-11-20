# 反射控制器 `y = Xβ + b` 的生理学解释

> 用途：为论文的 Discussion / 模型假设部分提供一套自洽、可引用的解释框架。
> 与 [`LASSO_SPARSITY_ANALYSIS_ARCHIVE.md`](LASSO_SPARSITY_ANALYSIS_ARCHIVE.md)
> 配合使用（后者给出“mask 为什么可能注定稠密”的解剖学证据）。
>
> **定位**：下面引用的神经生理学证据用于**设计动机**（为什么框架里要有
> 长度/力/抑制/相位项）与**结果的事后解释**（辨识出的结构与已知脊髓反射
> 是否自洽），**不是优化先验**；`β` 的数值由数据（LASSO 稳定性选择 + CMA-ES）
> 确定。证据强度分级见
> [`REFLEX_NEUROPHYSIOLOGY_DESIGN_AND_VALIDATION.md`](REFLEX_NEUROPHYSIOLOGY_DESIGN_AND_VALIDATION.md)。
>
> 符号：`X` 26 个特征（设计矩阵旧版记作 `A = X*M`，`M` 为单位阵，已移除，现统一用 `X`），`β` 26×9 反馈增益，`b` 1×9 偏置，`y` 9 块肌肉的 excitation。
> 相位组：stance（相位 0/1）、swing（相位 2/3/4）；每组一套 `β, b`。

---

## 0. 核心命题

**在稳定的周期性运动（平地行走）中，肌肉 excitation 由两件事共同决定：
脊髓节律发生器（CPG）给出“什么时候、按什么顺序、以多大强度”激活的基本模板，
感觉反馈在此基础上做实时修正。** 本模型的 `β` 对应前者之外的反馈通道权重，
`b` 对应“反馈参考点 + CPG 的偏置输出”，相位划分则是对二者做时序调制。
二者不是替代关系，而是**互补**：CPG 负责生成与协调，反馈负责适应负载、接触与扰动。

---

## 1. 模型与两个自由度

```
y_m(t) = Σ_i β(i,m) · X_i(t) + b(m)          m = 1..9, 每组一套 β、b
```

- 每个相位组内，`β` 与 `b` 是**常数**：这一组相位共享同一个“反射增益表”；
- 组间不同：stance 与 swing 各有一套，等价于一个**相位调制的增益调度**
  （phase-dependent gain scheduling）；
- 9 个 `b` 永远存在（不受 0-1 mask 约束），`β` 的 0-1 结构才是被辨识/优化的对象。

这个式子能且只能表达两类东西：**对当前状态（`X`）的瞬时反应**（`β` 项）和
**与当前状态无关的偏置驱动**（`b`）。两类恰好对应反馈与 CPG。

---

## 2. `β`：感觉反馈的通道权重

`X` 的 26 个特征分三族，每一族对应一类本体感受器：

| 特征族 | 生理对应 | 主要感受器 | 信号含义 |
|---|---|---|---|
| `lCEN`（9 个，归一化肌纤维长度） | 肌肉长度的变化与变化率（经 `fATN`、关节速度间接体现） | **肌梭 Ia** | 肌肉被拉伸 → Ia 放电增强 → 同源/异源 α 运动神经元兴奋 |
| `fATN`（9 个，归一化 MTU 力） | 肌肉-肌腱张力 | **高尔基腱器官 Ib** | 张力升高 → Ib 放电 → 抑制同源运动神经元、兴奋拮抗肌（功能性反转见 §5） |
| 关节角度/角速度（8 个：骨盆 tilt、髋、膝、踝） | 关节位置与运动 | **关节囊/韧带感受器（Ruffini、Pacinian 类）**，以及跨关节的肌梭长度变化 | 关节处于什么构型、运动多快；为跨关节协调提供状态量 |

**为什么说这组特征“基本完备”**：它覆盖了本体感觉的三大来源——**长度、力、关节
运动学**。对平地行走这类任务，缺少的主要是**皮肤触觉**（足底接触/压力）、
**前庭**与**视觉**；因此更严谨的表述是“**对下肢本体感觉反射而言基本完备**”，
而不是“对所有感觉输入完备”。这也解释了为什么模型在支撑相早期（足底触觉最
重要）的拟合误差偏大。

**逐项解释 `β` 的每一行**：

- `β(lCEN_j, m)`：肌肉 j 的长度变化对肌肉 m 兴奋性的影响。对角项（j = m）是
  **同源** Ia 投射，非对角项是**异源** Ia 投射（肌肉 A 的 Ia 影响肌肉 B 的 α
  神经元）。1993 年的经典工作表明异源投射在下肢**广泛分布**且可跨关节——
  这正是 `β` 非对角项应当普遍非零的解剖学依据（见稀疏度归档 §5）。
- `β(fATN_j, m)`：肌肉 j 的张力对肌肉 m 的影响（Ib 通路）。它把“负载”引入
  反射：支撑相负载大时抑制屈肌、促进伸肌，是承重与稳定性的基础。
- `β(关节量, m)`：关节层面的状态对肌肉 m 的调制，承担**跨关节协调**与
  **姿势参照**（例如骨盆 tilt 影响髋伸肌的活动）。

### 2.4 三类反馈的现代证据与建模含义

> 以下证据用于**设计动机**与**事后解释**；证据强度分“较硬 / 中 / 推断”，
> 写作时据此选择措辞（详见
> [`REFLEX_NEUROPHYSIOLOGY_DESIGN_AND_VALIDATION.md`](REFLEX_NEUROPHYSIOLOGY_DESIGN_AND_VALIDATION.md) §3.3/§4.4）。

**A. Ia 长度反馈（兴奋 + 交互抑制）**

- 同源单突触兴奋（膝跳反射）在几乎所有运动单位都能测到；
- **异名 Ia 兴奋广泛分布**，且可跨关节（Meunier 1993：股神经 → Sol 83%、→ TA 79%；
  GM→Bi 21/21），但“**widely distributed, generally quantitatively weak**”
  （异源峰一般 < 同源峰 20%）；
- 对拮抗肌则是**双突触抑制**（交互抑制，经 Ia 抑制性中间神经元，甘氨酸能，
  多约 1 ms）；**不对称**：跖屈肌→TA 的抑制强于反向（Crone 1987；Petersen 1999）。

*建模含义*：`β(lCEN)` 的非对角正项可以**广泛但小**；拮抗肌对（如
`gastroc/soleus → tib_ant`）**预期出现负项**。这为“稠密但多数很弱”的 mask
提供了动机，但**不规定**哪一项一定非零。

**B. Ib 力反馈（安静抑制、行走反转）**

- GTO 与肌纤维**串联**，编码肌腱张力；Ib **从不单突触**（至少一个中继）；
- **安静/站立**：力↑ → 抑制自己（自体抑制，负力反馈，保护性）；
- **行走（伸肌）**：同一 Ib 网络被运动状态重新配置，**翻转为兴奋**（正力反馈，
  帮助承重）；证据：脊髓猫（Pearson & Collins 1993）、清醒行走猫（Donelan 2009），
  人体为间接证据（Grey 2007：仅 SOL、晚支撑相）；
- **正反馈为何不发散**：环路增益 = 神经增益 K × 肌肉增益 M，M 随肌肉缩短而下降
  （力–长度曲线），使 K·M 自限；延迟反而帮助稳定（Prochazka 1997）。

*建模含义*：`β(fATN)` 的**符号应当与相位有关**——伸肌（尤其 SOL/GAS）在支撑相
以正为主，摆动相可能为负；这正是“stance/swing 各一套 `β`”的生理依据。

**C. 交互抑制与“抑制不止一种”**

抑制至少三种来源，写作时不要混为一谈：

1. **长度 Ia → 抑制拮抗肌**（交互抑制，双突触）；
2. **力 Ib → 安静时抑制自己**（自体抑制，负力反馈）；
3. **运动神经元回返抑制**（Renshaw，激活驱动）。

> 常见误解澄清（引自神经生理学文档 §5）：
> (i) **异名力反馈在生理上是可能的**（Ib 中间神经元多进多出），只是缺少像 Meunier
> 那样的人体全图；跨关节踝→膝 Ib 基本未证。
> (ii) **力负反馈真实存在**，Geyer 用负力反馈本身没有错；只是他那条
> “跖屈肌力→抑制 TA”的具体跨肌肉项是**功能等效抽象**（真实交互抑制是 Ia 长度驱动）。
> (iii) 协同肌之间的单突触连接**可能缺失**（GM→Sol = 0/9），所以“近协同”不等于“强连接”。

> 注意：`β` 是**功能性增益**，不是解剖学连接强度。它把突触前抑制、突触后
> 调制、运动单位募集阈值等全部“折合”成一个等效系数（详见 §5.5）。
> 因此文献里的百分比（如 Crone 的 14.9%）只能作为“通路存在/不对称”的证据，
> **不是可以移植到模型里的增益数值**。



---

## 3. `b`：反馈参考点 + 脊髓节律发生器

`b` 有两层含义，二者并不冲突：

### 3.1 反馈的参考/工作点（setpoint）

线性反馈 `β·X` 需要一个零点：当 `X` 处于什么水平时该肌肉“不激活”。`b` 就是
这个偏移。例如若 `lCEN ≈ 1`（静息长度）时希望肌肉静息，则 `b` 要抵消
`β(lCEN)·1`。从这个角度，`b` 编码的是“**在什么长度/张力下开始产生反射**”。
（旧版 UN.sto 路径把它打印为 `prestimulation`，是同一个意思。）

### 3.2 脊髓 CPG 的偏置输出

把 `b` 只当“参考点”会漏掉它的另一半功能：在**没有相应感觉输入变化**时，
肌肉仍然按节律激活——这就是 CPG。`b` 是常数项，它不依赖瞬时状态，因此天然
适合承载“**与相位绑定、与瞬时反馈无关**”的那部分驱动：

- **节律生成**：决定周期运动的基本时序（左右腿何时交替、支撑相/摆动相如何
  组织）；
- **肌群协调**：决定不同肌肉激活的**先后、持续时间、相对强度**；
- **相位模板**：在我们的实现里，相位划分（5 个相位 → 2 个组）给出**时序调制**，
  而每组的 `b` 给出**幅度/偏置调制**。二者合起来就是一张“按相位展开的激活模板”。

### 3.3 一个更准确的表述

> `β` 描述“**此刻的状态如何修正激活**”，`b` 描述“**这个相位本来应该激活多少**”。
> 前者是反馈，后者是节律与前馈模板；`β` 与 `b` 的分离，正是“反射”与“CPG”
> 在同一个线性模型里的分离。

---

## 4. CPG 与反馈：互补，而非替代

三条来自运动控制的经典观察支持这种分工（引用见 §9）：

1. **去传入（deafferented）动物仍能在跑台上产生节律性行走**——基本节律由
   脊髓 CPG 产生，不依赖感觉反馈；
2. **但失去感觉反馈后无法适应扰动与负载变化**（不能根据地面/障碍调整），
   说明反馈负责“适应性修正”；
3. **反射增益随任务/相位系统性变化**（站立 vs 行走 vs 跑），说明反馈不是
   固定增益，而是被节律网络调度（§5）。

映射到模型：`b`（+相位划分）提供第 1 条的基本节律；`β` 提供第 2 条的适应性；
`β`、`b` 的**分组相位依赖**提供第 3 条的调度。因此本模型不是“用反射替代 CPG”，
而是“**用一组稀疏的时变增益同时表达两者**”。

---

## 5. 反射通路本身受调制：增益不是常数

这是本框架最容易被忽略、但对论文最重要的一点：**同一份 Ia 输入，在不同运动
状态、不同相位、甚至不同身体条件下，产生的效果并不相同。**

### 5.1 相位依赖

踝关节 H-reflex（可视为 Ia 通路的整体增益）在行走的支撑相与摆动相幅度差异显著，
且与背景 EMG 解耦——说明增益调制发生在反射通路本身，而不只是运动神经元兴奋性
变化。定量上（Capaday & Stein）：**站立时 soleus H-reflex 可比行走支撑相大至
3.5 倍**（1986）；**跑步时 EMG 是行走的 2.4 倍、而 H-reflex 反而更小**（1987），
作者归因于突触前抑制，并把它解释为“自动增益补偿”（防止反射饱和/失稳）。
**我们按 stance/swing 分组各拟合一套 `β`，正是这一现象的最小实现。**

### 5.2 任务依赖

同一块肌肉的 H-reflex 在站立、行走、跑步中幅度不同，**Ia 增益大致为
站立 > 步行 > 跑步**；增益随任务“重新配置”。负荷也会缩放力反馈的贡献
（猫 ±25° 坡度：下坡 ~35% → 平地 ~39% → 上坡 ~54%，Donelan 2009）。
这意味着 `β` 不可能是一个跨任务的普适常数——**它是 task-dependent 的**。

### 5.3 方向反转（reflex reversal）

某些通路在特定相位/任务下从兴奋变为抑制（如猫行走时来自跖肌的 Ib 传入对
内侧腓肠肌的作用反转）。这说明矩阵的**符号**也受状态调制，而不只是幅度。
交互抑制同样有相位门控：**背屈肌→跖屈肌的抑制在摆动相大、支撑相缺如**
（Petersen 1999）。

### 5.4 机制：突触前抑制（presynaptic inhibition）

上述调制的共同机制之一是**对 Ia 终末的突触前抑制**：同一组 Ia 纤维在不同
运动状态下被不同程度地“门控”，从而改变其到各运动神经元池的有效增益。
Meunier 等（1993）在讨论异源投射时明确指出：这些“预接线”的通路可以在不同
任务间通过**下行控制突触前抑制**被“聚焦/选择”——这与我们的 `β` 在组间不同、
且需要按任务辨识完全一致。

### 5.5 对辨识的含义

- `β` 是**功能增益**（functional gain）：它把突触前抑制、突触后调制、运动单位
  募集等折合为一个等效系数。因此它**不能**直接与解剖学的 EPSP 幅度一一对应；
- 功能增益是 **task-dependent 与 subject-dependent** 的，所以：
  - 换任务/换受试者必须重新辨识 `β`；
  - 在仿真里，“让模型走起来”所需的增益与实测辨识值之间存在 real2sim gap，
    由 **CMA-ES 在 `β` 上继续优化**是合理的——它优化的正是这个功能增益；
  - 结构（mask）可以跨任务迁移，**数值不行**——这正是我们把 0-1 结构冻结、
    只让 CMA-ES 调数值的理由。

### 5.6 “谁在拧增益旋钮”：至少四只手

把反射增益想成反馈回路上的音量旋钮，至少四只手同时在拧：

1. **脊髓自己的运动模式（CPG/节律网络）**：进入“走路模式”就按相位调音量、甚至
   翻转 Ib 符号。**脊髓猫（去大脑）照样有 Ib 反转与相位调制** → 很大程度是脊髓
   自带的，但它是**动态网络状态**，不是固定查表；
2. **大脑下行**：按任务/意图调整；反射变化**先于 EMG 约 50 ms**，说明是“提前
   命令”而不是反馈的副产品；
3. **突触前抑制**：直接作用在 Ia 终末（传感器“输出端子”），相当于旋钮本身，
   脊髓和大脑都能拧；
4. **任务力学**：负载/长度/速度不同 → 传感器放电不同 → 有效反馈不同。

> 结论：**不是纯 CPG、也不是纯大脑，而是分层叠加**；每个任务里两者各占多少，
> 目前尚未完全分清（Petersen 等明确称 unresolved）。这正好解释了为什么我们
> 把一个“功能增益”交给 CMA-ES 去拟合，而不试图从解剖学推出它。

---

## 6. 可检验的预测

如果上述解释成立，应当能在仿真/数据中观察到：

1. **去掉反馈（`β = 0`）**：仍应产生基本节律（若 `b` 与相位划分合适），但抗扰动、
   变速、变坡度能力显著下降；
2. **stance 与 swing 的 `β` 应系统性不同**（而不是仅仅数值抖动）——可直接比较
   两组的 mask/β；
3. **增益矩阵应“广泛但弱”**：非对角项普遍非零、但多数幅值小（与 Meunier 1993
   的“widely distributed but quantitatively weak”一致）；
4. **跨受试者**：mask 的重叠应高于 `β` 数值的重叠（结构比增益更保守）；
5. **扰动实验**：在支撑相加短暂外力，模型应当通过 `fATN`（负载）项产生与
   实测同向的修正；若模型只靠 `b` 响应，则说明反馈项权重不足。

---

## 7. 与 mask 稀疏性/稠密性的关系

- 若 `β` 对应的是**广泛的异源 Ia 投射**（1993），那么**稠密是预期**，不是失败；
- 但“稠密”不等于“每一项都可重复”：π 分析显示相当一部分条目只有 40–60% 的
  可重复性，稳定性选择（π≥0.7 → 175 项）剔除的是这部分；
- 因此论文可以这样表述：
  > 解剖学证据给出一个**广泛连接、单个增益较小**的预期；本研究的辨识结果与
  > 这一**解剖学预期**在定性上一致（稠密、多数增益小），并通过多随机种子
  > 一致性分析区分出其中**可重复的核心**（π≥0.7）与**任务/噪声相关的外围**。
  > 注意这是**事后解释**，不是把解剖学当作优化先验。

---

## 8. 论文段落结构建议

1. **Model** 段：给出 `y = Xβ + b`，说明 `β`/`b` 分别对应反馈与 CPG；
2. **Interpretation** 段：用 §2–§3 的材料解释 26 个特征的三族与 `b` 的双重含义；
3. **Phase-dependent modulation** 段：§5.1–5.4，落到“为什么按相位分组、为什么
   增益要由 CMA-ES 优化”；
4. **Why dense** 段：§7 + Meunier 1993 的具体数字；
5. **Limitations** 段：`β` 是功能增益、缺少皮肤/前庭输入、单受试者。

---

## 9. 参考文献

> 完整、带 DOI/PMCID 且已对照原文核对的清单见
> [`REFLEX_NEUROPHYSIOLOGY_DESIGN_AND_VALIDATION.md`](REFLEX_NEUROPHYSIOLOGY_DESIGN_AND_VALIDATION.md) §9；
> 逐字引文见 `references/_review_Ib_force_feedback.md` 与
> `references/_review_reciprocal_inhibition.md`。

**核心（本仓库已核对）**

- Meunier S, Pierrot-Deseilligny E, Simonetta M. Pattern of monosynaptic
  heteronymous Ia connections in the human lower limb. *Exp Brain Res.*
  1993;96:534–544.

**反射的相位/任务调制**

- Capaday C, Stein RB. Amplitude modulation of the soleus H-reflex in the human
  during walking and standing. *J Neurosci.* 1986;6:1308–1313.
- Capaday C, Stein RB. Difference in the amplitude of the human soleus H reflex
  during walking and running. *J Physiol.* 1987;392:513–522.
- Zehr EP, Stein RB. What functions do reflexes serve during human locomotion?
  *Prog Neurobiol.* 1999;58:185–205.
- Yang JF, Stein RB, James KB. Contribution of peripheral afferents to the
  activation of the soleus muscle during walking in humans.
  *Exp Brain Res.* 1991;87:679–687.

**反射反转与 Ib 通路**

- Forssberg H, Grillner S, Rossignol S. Phase dependent reflex reversal during
  walking in chronic spinal cats. *Brain Res.* 1975;85:103–107.
- Pearson KG, Collins DF. Reversal of the influence of group Ib afferents from
  plantaris on activity in medial gastrocnemius muscle during locomotor
  activity. *J Neurophysiol.* 1993;70:1009–1017.
- Prochazka A, Gillard D, Bennett DJ. Positive force feedback control of
  muscles. *J Neurophysiol.* 1997;77(6):3226–3236.
- Grey MJ, Nielsen JB, Mazzaro N, Sinkjaer T. Positive force feedback in human
  walking. *J Physiol.* 2007;581(Pt 1):99–105.
- Sinkjaer T, Andersen JB, Ladouceur M, Christensen LO, Nielsen JB. Major role
  for sensory feedback in soleus EMG activity in the stance phase of walking
  in man. *J Physiol.* 2000;523(Pt 3):817–827.
- Donelan JM, McVea DA, Pearson KG. Force regulation of ankle extensor muscle
  activity in freely walking cats. *J Neurophysiol.* 2009;101(1):360–371.

**突触前抑制与脊髓回路**

- Hultborn H, Meunier S, Pierrot-Deseilligny E, Shindo M. Changes in presynaptic
  inhibition of Ia fibres at the onset of voluntary contraction in man.
  *J Physiol.* 1987;389:757–772.
- Pierrot-Deseilligny E, Burke D. *The Circuitry of the Human Spinal Cord:
  Spinal and Corticospinal Mechanisms of Movement.* Cambridge University Press;
  2012.

**交互抑制**

- Crone C, Hultborn H, Jespersen B, Nielsen J. Reciprocal Ia inhibition between
  ankle flexors and extensors in man. *J Physiol.* 1987;389:163–185.
- Petersen N, Morita H, Nielsen J. Modulation of reciprocal inhibition between
  ankle extensors and flexors during walking in man. *J Physiol.*
  1999;520(Pt 2):605–619.
- Pierrot-Deseilligny E, Morin C, Bergego C, Tankov N. Pattern of group I fibre
  projections from ankle flexor and extensor muscles in man. *Exp Brain Res.*
  1981;42:337–350.

**反射框架的设计来源**

- Geyer H, Herr H. A muscle-reflex model that encodes principles of legged
  mechanics produces human walking dynamics and muscle activities.
  *IEEE Trans Neural Syst Rehabil Eng.* 2010;18(3):263–273.

**CPG**

- Grillner S. Neurobiological bases of rhythmic motor acts in vertebrates.
  *Science.* 1985;228:143–149.
- McCrea DA, Rybak IA. Organization of mammalian locomotor rhythm and pattern
  generation. *Brain Res Rev.* 2008;57:134–146.
- Dietz V. Do human bipeds use quadrupedal coordination? *Trends Neurosci.*
  2002;25:462–467.

**待补充的支撑点**

- “脊髓网络调节反射通路本身”的直接人体证据（除 H-reflex 调制外）；
- 反射增益的**受试者间变异**（subject-dependent）的定量文献；
- 关节感受器（Ruffini/Pacinian）在步行中的具体贡献——目前 `β` 的关节项只能
  说“功能性”，不宜过度归因到某一种感受器。

**表述边界（写作时务必注意，摘自神经生理学文档 §7.4）**

1. 不要写“Ib 力反馈贡献 SOL 支撑相 EMG ~50%”：Sinkjaer 2000 的 ~50% 是**传入
   总贡献**（作者倾向 group II）；Ib 特异数值来自猫（Donelan 约 35–54%），
   人体力相关证据是 Grey 2007（仅 SOL、晚支撑、间接）。
2. 不要照抄 Geyer 对 Christensen 2001 的描述：原文是**支撑相大、摆动相小**
   （Geyer 写反了）；可辩护的是“短潜伏期 M1/Ia 成分摆动有、支撑无”。
3. Geyer 的 `F−`（跖屈肌力→抑制 TA）是**功能抽象**：真实交互抑制是 Ia 长度驱动，
   经典 Ib 符号相反；但**力负反馈本身真实存在**，不要说它不存在。
4. HAM→HFL 的抑制是**假设**（髋屈肌运动神经元池从未被采样），宜作消融候选。
5. 跨关节踝→膝 Ib **无证据**；异名 Ib 仅在踝伸肌群内有证据（猫）。
6. 支撑相“跖屈肌→背屈肌大抑制”是**推断**（Petersen 明言无直接证据）。
7. `β` 是功能增益，文献百分比（如 Crone 14.9%）是“通路存在/不对称”的证据，
   **不是可移植的模型增益**。
