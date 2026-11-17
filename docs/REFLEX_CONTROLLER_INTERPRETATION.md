# 反射控制器 `y = Xβ + b` 的生理学解释

> 用途：为论文的 Discussion / 模型假设部分提供一套自洽、可引用的解释框架。
> 与 [`LASSO_SPARSITY_ANALYSIS_ARCHIVE.md`](LASSO_SPARSITY_ANALYSIS_ARCHIVE.md)
> 配合使用（后者给出“mask 为什么可能注定稠密”的解剖学证据）。
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

> 注意：`β` 是**功能性增益**，不是解剖学连接强度。它把突触前抑制、突触后
> 调制、运动单位募集阈值等全部“折合”成一个等效系数（详见 §5.5）。

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
变化。**我们按 stance/swing 分组各拟合一套 `β`，正是这一现象的最小实现。**

### 5.2 任务依赖

同一块肌肉的 H-reflex 在站立、行走、跑步中幅度不同；增益随任务“重新配置”。
这意味着 `β` 不可能是一个跨任务的普适常数——**它是 task-dependent 的**。

### 5.3 方向反转（reflex reversal）

某些通路在特定相位/任务下从兴奋变为抑制（如猫行走时来自跖肌的 Ib 传入对
内侧腓肠肌的作用反转）。这说明矩阵的**符号**也受状态调制，而不只是幅度。

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
  > 解剖学证据预测一个**广泛连接、单个增益较小**的反射矩阵；本研究的辨识结果
  > 与这一先验一致，并通过多随机种子一致性分析区分出其中**可重复的核心**
  > （175/234 项，π≥0.7）与**任务/噪声相关的外围**。

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

> ⚠️ 除 Meunier 1993（本仓库 `references/BF00234121.pdf`，页码已核对）外，
> 其余条目的卷期/页码请在投稿前用数据库核对一遍。

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

**突触前抑制与脊髓回路**

- Hultborn H, Meunier S, Pierrot-Deseilligny E, Shindo M. Changes in presynaptic
  inhibition of Ia fibres at the onset of voluntary contraction in man.
  *J Physiol.* 1987;389:757–772.
- Pierrot-Deseilligny E, Burke D. *The Circuitry of the Human Spinal Cord:
  Spinal and Corticospinal Mechanisms of Movement.* Cambridge University Press;
  2012.

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
