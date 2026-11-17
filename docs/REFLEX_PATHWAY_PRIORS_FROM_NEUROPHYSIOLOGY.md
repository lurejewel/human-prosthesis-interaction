# 反射通路存在性的神经生理学先验：Meunier 1993、Ib 力反馈与交互抑制

> **用途**：为 muscle-reflex 肌骨预测仿真（Geyer & Herr 2010 风格，OpenSim/Moco + CMA-ES）
> 提供“**哪些反射通路存在 / 不存在**”的解剖学与电生理证据，用于约束反射增益矩阵 `β` 的
> **结构（mask）与符号**，从而降低优化维度。
>
> **本文不讨论增益大小**。增益由 LASSO 稳定性选择 + CMA-ES 在仿真中优化获得（见
> [`LASSO_STABILITY_IDENTIFICATION.md`](LASSO_STABILITY_IDENTIFICATION.md)）。本文只回答
> “**这条通路该不该留在优化空间里、该允许什么符号**”。
>
> 配套文档：[`REFLEX_CONTROLLER_INTERPRETATION.md`](REFLEX_CONTROLLER_INTERPRETATION.md)
> （`y = Xβ + b` 的生理学解释）。本文与其 §7“稠密预期”互补：那里论证“为什么 β 倾向稠密”，
> 这里论证“**哪些具体条目可以可靠地剪掉**”。

---

## 0. 符号与模型对齐

本文所有讨论对齐本仓库的反射控制器：

```
y_m(t) = Σ_i β(i,m) · X_i(t) + b(m)
```

> 记号：设计矩阵直接用 `X`（旧版写作 `A = X*M`，`M` 为单位阵，已移除）。

- `β` 为 **26 特征 × 9 肌肉** 的**有符号**增益矩阵，按相位组（stance / swing）各一套。
- 特征行：
  - **1–9：`lCEN`**（归一化肌纤维长度）→ 生理对应 **肌梭 Ia**（长度/拉伸反馈）。
  - **10–18：`fATN`**（归一化 MTU 力）→ 生理对应 **高尔基腱器官 Ib**（力/负荷反馈）。
  - **19–26：关节角度/角速度** → 关节感受器 + 跨关节肌梭状态量。
- 9 块肌肉（列）顺序：`hamstrings, bifemsh, glut_max, iliopsoas, rect_fem, vasti, gastroc, soleus, tib_ant`。
- `mask(i,m)=false` → **结构零**，`β(i,m)` 恒为 0，**不进入 CMA-ES**（真正降维）。
- `mask(i,m)=true` → 自由参数，可正可负。

**关键记号**：`β(lCEN_j, m)` 表示“肌肉 j 的长度对肌肉 m 兴奋性的影响”。
- `> 0`：j 被拉伸 → 兴奋 m（**兴奋性长度反馈**；j≠m 时即**异名 Ia 兴奋**，Meunier 1993 所测）。
- `< 0`：j 被拉伸 → 抑制 m（**抑制性长度反馈**；功能上接近**交互抑制**，由 Ia 抑制性中间神经元介导）。
- `= 0`：无长度耦合。

> ⚠️ 因为 `β` 是**单个有符号系数**，“把正反馈硬编码为零、但保留负反馈”在实现上
> **不是结构零**，而是**符号约束 `β ≤ 0`**（保留参数、限制取值）。这一点在 §2.3 展开。

---

## 1. Meunier 1993 与 muscle-reflex 仿真的关系

### 1.1 两者属于同一问题，但证据层级不同

- **Geyer 2010**：用肌肉感觉反馈（长度、力、关节、负荷）直接驱动肌肉刺激，让步行从
  神经—肌肉—地面交互中自组织。它是一个**功能性 / 现象学**反射控制模型。
- **Meunier 1993**（`references/BF00234121.pdf`，*Exp Brain Res* 96:534–544）：用单运动单位
  PSTH 系统绘制人类下肢 **肌梭 Ia 传入到各肌肉运动神经元的单突触兴奋性投射**。它是
  **解剖/电生理连接组**层面的证据。

二者关系：

```
Geyer muscle-reflex control
    ├── Ia 长度反馈（同名 + 异名兴奋）   ← Meunier 直接相关
    ├── Ib / 力反馈（F+）               ← Meunier 不相关（见 §3）
    ├── 交互抑制 / 负反馈               ← Meunier 不相关（见 §4）
    ├── 负荷 / 接触相位切换             ← Meunier 不相关
    └── 上层 / 相位增益调制             ← Meunier 不相关
```

**结论**：Meunier 不是 Geyer 模型的完整神经基础，也不提供增益数值；它只为
**“Ia 兴奋性长度反馈”这一子矩阵**提供**存在性、方向、符号**的先验。

### 1.2 Geyer 2010 各反射规则与 Meunier 的对应

| Geyer 2010 规则 | 生理机制 | Meunier 能否支撑 |
|---|---|---|
| SOL/GAS/VAS 正力反馈 `F+` | Ib / 负荷反馈 | **否**（Meunier 测 Ia，不测 Ib） |
| TA 正长度反馈 `L+` | 同名 Ia stretch reflex | **间接**（同名 Ia 单突触兴奋普遍存在） |
| SOL 力反馈抑制 TA（`F−`） | 双突触 Ia 交互抑制 | **否**（Meunier 不测抑制） |
| HFL 正长度反馈 `L+` | 髋屈肌同名 stretch reflex | **有限**（未直接研究髋屈肌靶） |
| HAM 长度反馈抑制 HFL（`L−`） | 跨肌肉抑制（Geyer 自述缺乏 Ia 交互抑制证据） | **否** |
| GLU/HAM 正力反馈 `F+` | 力/负荷反馈 | **否** |
| VAS 膝角阈值抑制、对侧负荷抑制 VAS、躯干 PD | 关节/负荷/前庭 | **否** |

### 1.3 Meunier 的核心发现（与建模相关）

1. 人类下肢存在**广泛**的单突触异名 Ia 兴奋，**可跨踝—膝关节**。
2. 连接**多数弱**（异名峰多为同名峰的 7%–22%），少数强（GM→Bi ≈54%，F→Sol ≈39%）。
3. **方向不对称**：Sol→GM 有（弱），GM→Sol 未检出；Sol↔Per、Sol↔Q 双向。
4. **同关节协同肌不一定强耦合**（GM→Sol 缺失）；**看似拮抗的肌肉间也可有兴奋**（DP→GM）。
5. 这些“预接线”通路可被**下行突触前抑制**按任务选择/聚焦（与相位调制增益一致）。

> 对 mask 的含义：**一般预期是“广泛连接”（稠密），但存在一批可重复的“具体缺失”**，
> 后者正是可用来降维的剪枝点（§2）。

---

## 2. 核心问题：未检出阳性连接，能否硬编码为零？

> 用户问题：比起“某通路存在”，更关心“**某通路不存在**”，因为可直接减少优化维度。
> (a) 能否把未检出阳性连接的路路的**长度反馈**硬编码为 0？
> (b) 能否只把**长度正反馈**硬编码为 0、但保留**长度负反馈（抑制）**的可能？

### 2.1 先厘清：“未检出” ≠ “不存在”

Meunier 表 1 中的 `0/n` 含义是：**在该实验条件下，未检测到显著的单突触 Ia 兴奋**。
它**不是**“该通路解剖学上不存在”的证明。原因有四：

1. **统计功效有限**：每个组合只测 4–25 个运动单位（`0/4`、`0/6`、`0/9`、`0/12`）。
   弱连接可能漏检。`n` 越小，“缺失”的证据越弱。
2. **只测单突触兴奋**：判据是中枢潜伏期差 ≤0.2 ms。它**主动排除**了：
   - **少突触/多突触兴奋**（潜伏期更长的兴奋性效应）；
   - **抑制**（交互抑制是双突触，本就不会出现在单突触兴奋峰里）；
   - **group II** 长度反馈。
3. **状态受限**：坐位、<5% MVC、低阈值运动单位。某通路可能在步行、高阈值运动单位或
   不同相位下才表达。
4. **突触前门控**：通路可能存在但在测试条件下被突触前抑制压低到检测阈以下。

> 因此：**“未检出单突触 Ia 兴奋” 只能否定“快速单突触 Ia 兴奋”这一特定成分，
> 不能否定整条长度耦合（尤其是抑制性成分）。**

### 2.2 关键区分：正反馈 vs 负反馈，结构零 vs 符号约束

Meunier 的阴性结果**只约束兴奋性（正）Ia 长度反馈**，对抑制性（负）长度反馈**完全没有约束**。
而交互抑制（负长度反馈）恰恰是**拮抗肌之间预期存在**的通路（§4）。

由此得到两条不同强度的操作：

| 操作 | 含义 | 是否被 Meunier 阴性支持 | 实现 |
|---|---|---|---|
| **整条长度反馈 = 0** | `β(lCEN_j,m)=0`，正负都禁 | **不支持**（会误删可能的抑制） | `mask=false`（真降维） |
| **只禁正、允许负** | `β(lCEN_j,m) ≤ 0` | **较合理**（否定单突触 Ia 兴奋，保留交互抑制） | **符号约束**（保留参数，限制上界=0） |

**直接回答**：
- **(a) 不建议**仅凭 Meunier 阴性就把整条长度反馈设结构零——除非你另有理由相信该路路
  连抑制也没有。对拮抗肌对（如 gastroc/soleus ↔ tib_ant），抑制是预期存在的，设结构零会删错。
- **(b) 方向正确**：把**兴奋性单突触 Ia 长度反馈**禁为零、保留负反馈，是生理上更站得住的做法。
  但在你的有符号 `β` 模型里，这应实现为 **`β ≤ 0` 的符号约束**，而**不是 `mask=false`**。
  二者对 CMA-ES 维度的影响不同（§2.6）。

### 2.3 推荐策略（分级，按证据强度）

把 Meunier 阴性当作**软先验 / 符号约束**，而非硬解剖真理：

1. **最强剪枝**：`gastroc → soleus` 的正 Ia 长度反馈。
   - Meunier `GM→Sol = 0/9`，且被 Pierrot-Deseilligny 1981、Mao 1984 **独立重复**。
   - 这是全文最可靠的“缺失”。可对该路路的**正成分**设强先验/符号约束。
2. **中等剪枝**（`n≥8` 或机制清晰）：`tib_ant → vasti`（DP→Q 0/9）、
   `tib_ant → bifemsh`（DP→Bi 0/8）。
3. **弱剪枝**（`n≤6`，谨慎）：`gastroc → tib_ant`（0/6）、`soleus → tib_ant`（0/6）、
   `tib_ant → soleus`（0/4）。这些多为拮抗肌对，**预期负（交互抑制）**，故应**禁正留负**，
   而非整条置零。
4. **不要剪**：所有 Meunier 阳性路路（§2.4 表），以及 Meunier **未覆盖**的路路（§2.5）。
5. **敏感性检验**：优化后，释放少数被约束路路，看 CMA-ES 是否把它们推离零。若始终贴近零，
   约束合理；若显著非零，说明该“缺失”先验过强，应回退。

### 2.4 映射到你的 9 肌肉 `β(lCEN)`：可剪枝 / 应保留

> 源神经 → 模型源肌肉的近似映射：`Inf Sol≈soleus`，`GM≈gastroc(内侧头)`，
> `DP≈tib_ant(胫前肌群)`，`F≈vasti+rect_fem+iliopsoas(混合，含髋屈肌，归属模糊)`，
> `SP≈腓骨肌(模型无此肌肉)`。靶肌肉：`Sol→soleus`，`GM→gastroc`，`TA→tib_ant`，
> `Q(VL)→vasti`，`Bi(短头)→bifemsh`，`ST→hamstrings`，`Per brev→模型无`。

**A. Meunier 阴性 → 候选“禁正留负”路路（兴奋性 Ia 长度反馈）**

| Meunier 路路 | 检出 | 模型 `β(lCEN_j, m)` 正成分 | 建议 | 置信度 |
|---|---|---|---|---|
| GM→Sol | 0/9 | `gastroc → soleus` | 禁正（可强约束） | **高**（多研究重复） |
| DP→Q | 0/9 | `tib_ant → vasti` | 禁正 | 中 |
| DP→Bi | 0/8 | `tib_ant → bifemsh` | 禁正 | 中 |
| GM→TA | 0/6 | `gastroc → tib_ant` | 禁正留负（拮抗） | 中–低 |
| InfSol→TA | 0/6 | `soleus → tib_ant` | 禁正留负（拮抗） | 中–低 |
| DP→Sol | 0/4 | `tib_ant → soleus` | 禁正留负（拮抗） | 低 |
| GM→Per | 0/12 | （腓骨肌，模型无） | 不适用 | — |
| SP→TA / SP→Q | 0/6 / 0/9 | （腓骨肌源，模型无） | 不适用 | — |

**B. Meunier 阳性 → 应保留（不要剪）的兴奋性路路**

| Meunier 路路 | 检出 / 相对强度 | 模型对应（正成分保留） |
|---|---|---|
| GM→Bi | 21/21 / 54% | `gastroc → bifemsh`（最强） |
| F→Sol | 49/70 / 39% | `vasti/rect_fem/iliopsoas → soleus`（F 归属模糊） |
| F→TA | 55/63 / 21% | `vasti/rect_fem/iliopsoas → tib_ant` |
| F→GM | 10/12 / 7% | `vasti/rect_fem/iliopsoas → gastroc` |
| InfSol→Bi | 10/16 / 22% | `soleus → bifemsh` |
| InfSol→GM | 11/16 / 12% | `soleus → gastroc`（弱；与 GM→Sol 缺失不对称） |
| InfSol→Q | 17/33 / 11% | `soleus → vasti` |
| GM→Q | 5/19 / 13% | `gastroc → vasti`（弱、低频） |
| InfSol→ST | 6/20 / 14% | `soleus → hamstrings` |
| GM→ST | 11/25 / 12% | `gastroc → hamstrings` |
| DP→ST | 8/23 / 11% | `tib_ant → hamstrings` |

### 2.5 Meunier 的覆盖盲区（切勿过度剪枝）

Meunier 表 1 的覆盖是**不对称**的：它测“小腿神经 / 股神经 → 小腿或大腿靶”，
但**没有把大腿肌肉当作异名源**去测。因此下列模型路路**几乎没有 Meunier 先验**，
**不能**因为“表里没有”就置零：

- 以 `hamstrings`、`bifemsh` 为**源**的异名长度反馈（Meunier 只把它们当靶，未当源测）；
- 以 `glut_max`、`iliopsoas`、`rect_fem` 为**靶**的路路（Meunier 未记录髋肌靶）；
- `rect_fem` 作为独立源（Meunier 的 Q=股外侧肌，股神经混合了髋屈肌，归属模糊）；
- 任何涉及**腓骨肌**的路路（模型无此肌肉，Meunier 的 Per brev 数据无法映射）。

> 这些路路的取舍应交给 **LASSO 稳定性选择 + 仿真优化**，而不是 Meunier。

### 2.6 实现注意：结构零 vs 符号约束（对 CMA-ES 维度的影响）

- **`mask=false`（结构零）**：`β=0` 恒定，**从 CMA-ES 向量中删除** → 真正降维。
  仅当你相信该路路**正负都不存在**时才用。Meunier 阴性**不足以**单独支持这一点。
- **符号约束 `β≤0`（或 `β≥0`）**：**保留参数**，只限制取值范围 → **不降维**，但缩小搜索空间、
  注入生理先验。当前 `load_lasso_reflex_controller.m` 的 0-1 mask **不支持符号界**，
  需要在 CMA-ES 的 `bounds`/`initPara` 上额外加上界（例如把该条目上界设为 0）。
- **务实折中**：若你的 `β(lCEN)` 子矩阵**只打算表达 Ia 兴奋性 stretch reflex**，而把抑制
  交给 `fATN` 负项或专门的抑制项，那么对 Meunier 阴性路路设 `mask=false` 是**可接受的建模约定**，
  但必须在论文中**显式声明**“此处忽略了这些路路的抑制性长度耦合”。否则对拮抗肌对会失真。

### 2.7 与“稠密预期”的 reconcile

[`REFLEX_CONTROLLER_INTERPRETATION.md`](REFLEX_CONTROLLER_INTERPRETATION.md) §7 论证
“β 倾向稠密”（Meunier 的广泛弱异名投射）。本文不矛盾：

> **一般规律是广泛连接（稠密），但存在一批可重复的具体缺失（洞）。**
> 剪掉这几个well-supported的“兴奋性洞”，与“其余广泛连接”完全相容。
> 而且**负（抑制）子矩阵由交互抑制文献支配（§4），不由 Meunier 支配**——
> 所以完整的 `β(lCEN)` 仍可在符号意义上稠密，只是**正成分有特定的洞**。

---

## 3. Ib 力反馈证据（Prochazka / Sinkjaer / Grey / Pearson / Donelan）

> 完整逐字引文、方法细节与页码见 [`references/_review_Ib_force_feedback.md`](../references/_review_Ib_force_feedback.md)。
> 本节只提炼**对 `β(fATN)` 子矩阵存在性/符号先验**有用的结论。所有引文已对下载原文逐字核对。

### 3.0 访问状态与本地文件

| 文献 | 访问 | 本地文件（`references/`） |
|---|---|---|
| Prochazka 1997 | **全文 + PDF** | `Prochazka-1997-...pdf` / `_FULLTEXT.txt` |
| Grey 2007 | **全文**（PMC HTML，无 PDF） | `Grey-2007-..._FULLTEXT.txt` / `Grey-2007-PMC2075215.html` |
| Sinkjaer 2000 | **全文**（PMC HTML，无 PDF） | `Sinkjaer-2000-..._FULLTEXT.txt` / `Sinkjaer-2000-PMC2269822.html` |
| **Pearson & Collins 1993** | **仅摘要（付费墙）** | `Pearson-Collins-1993-..._ABSTRACT-ONLY.txt` |
| Donelan 2009 | **全文 + PDF** | `Donelan-2009-...pdf` / `_FULLTEXT.txt` |
| Prochazka 2002（可选） | **全文 + PDF**（机构访问） | `Prochazka-2002-...pdf` / `_FULLTEXT.txt` |

> ⚠️ Pearson & Collins 1993 全网封闭（OpenAlex `any_repository_has_fulltext:false`），仅得截断摘要；
> 涉及它的细节均标注为“摘要”或来自 Prochazka 1997 / Donelan 2009 的**二手佐证**。

### 3.1 经典 Ib 通路：从来不是单突触

- 经典回路：**高尔基腱器官 → Ib 传入 → Ib 抑制性中间神经元 → 同源（及协同肌）α 运动神经元**，
  双突触“自体抑制”（autogenic inhibition），静息时为**负力反馈**。
- **Ib 作用完全经中间神经元，无单突触 Ib→运动神经元通路**：
  > "tendon organ Ib reflex action, which is **entirely mediated through interneurons**"（Prochazka 1997, p3226）
- Ib 中间神经元**汇聚/发散**：接受多块肌肉的 Ib（也有 Ia）传入 + 下行（如红核脊髓）输入
  → **符号与增益可被中枢调制**（Prochazka 1997, p3226）。
- **行走中反转**：静息抑制 → 运动时兴奋；有效通路**比双突触更长**（猫 30–40 ms，人卸载反应 ~64 ms），
  与脊髓伸肌“半中枢”节律网络整合（Sinkjaer 2000 明确认为双突触太快，解释不了 64 ms）。

> **建模含义**：任何 `F+` 块都是**集总、延迟、经中间神经元**的通路，不应宣称为单突触或简单双突触 Ib 反射；
> 数十毫秒的长延迟有生理依据，不是凑数。

### 3.2 正力反馈 F+ 与“抑制→兴奋”反转的证据链（强→弱）

| 证据 | 来源 | 确立了什么 |
|---|---|---|
| **状态依赖符号反转，排除 Ia** | Pearson & Collins 1993（**仅摘要**） | 跖肌→MG 在**运动时兴奋**、同一刺激串**静息时抑制** tonic MG；跖肌振动（Ia 选择性）**两者都不产生** → 效应属 Ib。潜伏期 30–50 ms，靶为**伸肌半中枢**。异名、脊髓猫。 |
| **清醒猫因果 + 定量** | Donelan 2009（全文） | MG 活动跟随**力**而非长度/速度；预测 **Ib** 匹配、Ia/II 不匹配；移动 peg 扰动给因果支持；loop gain 平地 0.26/0.51，恒 <1；**神经增益 K 长度无关、肌肉增益 M 长度依赖**。**同源 MG，异名通路已被去神经排除**。 |
| **控制论可行性** | Prochazka 1997（全文） | 真实人体肌肉 + **人工反馈回路**（hybrid）：F+ 稳定补偿负载；**Gf>1 仍稳定**（力–长度曲线随肌肉缩短自动降增益）；**40 ms 延迟反而稳定**。**非解剖学证明**。 |
| **明确开关 + 可实现模型** | Prochazka 2002（全文） | "when locomotion starts, there is a switch from inhibition … to excitation"；给出 Ib(s) 传递函数、35 ms 延迟、激活门控。 |
| **人体步行相关** | Grey 2007（全文） | 晚支撑相 SOL 卸载反应随**跟腱力**增大、随踝运动学**反向**变化（±4% 坡度）。**间接**。 |

**为何正反馈不发散**（仿真关键）：神经增益 K 大致长度无关，但肌肉增益 M 随缩短下降（力–长度特性），
故 **loop gain K·M 自限趋于 1**（Donelan 2009 实测此分解）；且**正反馈回路中的延迟是稳定化的**
（Prochazka 1997）。警告：对**不动/等长**负载，F+ 可能失稳——“可能只存在于腿伸肌这类始终惯性负载的肌肉”
（Prochazka 1997, p3234）。

### 3.3 人体步行证据：Grey 2007 与 Sinkjaer 2000 必须分开引用

**⚠️ 最易误引的一点：Sinkjaer 2000 并不证明正力反馈。**

- **Sinkjaer 2000**：20 人，卸载使 SOL EMG 降 **~50%（最大 51±7%）**， onset **63.6±3.5 ms**。
  该反应**非交互抑制**（CPN 阻滞不变）、**非 Ia**（缺血阻滞不变，而 Ia 短潜伏期反射已消失）、
  **非大皮肤传入**。但作者明确指出：反应随支撑相力增大而**减小**——**与 Ib 力反馈预测相反**；
  倾向 **group II**；摘要结论为 "**Group II and/or group Ib**"。
  → 它确立的是**传入总贡献的量级**，**不是 Ib/正力反馈**。
- **Grey 2007**：21 人，±4% 坡度，**晚支撑相（60–80%）**扰动 + **直接测跟腱力（ATF）**。
  SOL 卸载反应 **2.0/2.2/2.7 μV·s**、ATF 卸载反应 **0.009/0.021/0.047**（下坡/平地/上坡），
  而踝位移 **3.6/3.1/2.6 deg**、速度 **39/35/31 deg/s** **反向**变化 → 力信号而非长度信号。
  **MG 反应未分析**（多数受试者缺如/极短）→ **这是比目鱼肌结果，不是腓肠肌结果**。
  结论："compelling, albeit **indirect**, evidence … excitatory group Ib pathway … **late stance**"。
  → 它确立的是**力相关性**（且仅限 SOL、晚支撑相）。

> **正确引用**：
> - “传入反馈贡献 SOL 支撑相 EMG 约 50%” → Sinkjaer 2000 ✅
> - “Ib 力反馈贡献约 50%” → ❌ **无证据**。Ib 特异数值是**猫**的 ~30%（Donelan 摘要/讨论）至 ~39%（Donelan 结果）。
> - Grey 2007 与 Sinkjaer 2000 **互补**：后者提出“力依赖”质疑（早/中支撑相反力增大而减小），
>   前者移到晚支撑相 + 直接测腱力来回应。成对阅读才完整。

### 3.4 异名 Ib 连接：哪些存在、哪些不存在

**有证据（强）**：
- **踝伸肌协同肌之间的异名 Ib**：Pearson & Collins 1993 本身就是异名（**跖肌→MG**，两块不同肌肉，
  但**都是踝伸肌**）；Prochazka 1997 p3234 称之为“selective Ib reflex action from one muscle to another”。
- Ib 中间神经元汇聚/发散；fictive locomotion 中**分级异名 Ib 兴奋**（Gossard 1994，经 Prochazka 1997 转述）。

**无证据 / 含糊**：
- **踝伸肌 → 膝肌（跨关节 Ib）：这六篇均未证明。** 跖肌与 MG 同属踝伸肌群，是**群内异名**，非跨关节。
  唯一指向踝/趾→膝的是**参考文献条目**（非结果）：**Bonasera & Nichols 1994**（J Neurophysiol 71:1096–1110，
  未下载）；另有 Harrison/Jankowska 1983。注意 Nichols 传统的 “heterogenic” 指**力学（长度/力）反射耦合**，
  不自动等于 Ib 特异。
- Grey 2007 **不区分**同源/异名（"homonymous **and/or** heteronymous"）；Donelan 2009 **刻意排除异名**（去神经）。
- **无任何人体异名 Ib 图谱**。

> **对 `β(fATN)` 的含义**：若你想在模型里加“踝伸肌力 → 膝伸肌（VAS）”的 Ib 路路，
> **这六篇不支持**；需要 Bonasera & Nichols 1994 / Jankowska 中间神经元文献，且应标注为**外推**。

### 3.5 对 `β(fATN)` 子矩阵的存在性 / 符号先验

| 路路 | 证据强度 | 建议 |
|---|---|---|
| **SOL、GAS 同源正力反馈（F+）** | **最强**（猫反转 + 清醒猫定量 + 人晚支撑相力依赖 + 稳定性理论） | 保留正力反馈；同源/自体是可辩护默认 |
| **VAS 正力反馈（F+）** | **控制原理外推**，非已证明的人体 Ib 通路 | 可保留（Geyer 靠它使分段腿柔顺、防膝塌），但按**控制原理**而非神经生理学论证 |
| **踝伸肌力 → 膝（跨关节 Ib）** | **本组文献无证据** | 不要据此声称存在；若加入须标注外推 |
| **异名 Ib（踝伸肌群内）** | 有（跖肌→MG，猫） | 可考虑，但人体证据缺如 |
| **符号** | **状态依赖反转** | `β(fATN)` 符号可相位调制：支撑相兴奋、静息/摆动相可能抑制 |

**稳定性先验**（供 CMA-ES 边界/ sanity check，非增益赋值）：
- 稳态 **loop gain < 1**（Donelan 2009；猫 MG 平地 0.26/0.51 可作 sanity 范围，但是猫 MG 非人 SOL）。
- 依赖 Hill 型肌肉**力–长度曲线**自动降增益来稳定（Prochazka 1997）；等长/不动负载会失稳。
- **长延迟生理可行且助稳**（猫 Ib 30–40 ms，人卸载 ~54–64 ms；Geyer 的 20 ms 是通用 H 反射延迟，非 Ib 延迟）。

**⭐ 可直接实现的升级**（Prochazka 2002, p363）：
```
Ib(s) = KIb · F · (s+0.15)(s+1.5)(s+16) / [(s+0.2)(s+2)(s+37)]
```
35 ms 延迟、**仅在受体肌肉激活时门控**、平均反射贡献 ~15%。比 Geyer 的“原始力增益”更生理，可直接引用。

### 3.6 Geyer 的 F+ 是现象学控制律，不是字面 Ib 回路

已对本地 Geyer 2010 PDF 核对：F+ = **prestimulation + 延迟、增益化的肌肉自身力**（同源 SOL/GAS/VAS）。
- 增益初值来自 **Geyer 2003 bouncing gaits**（ref [29]），非 Ib 生理；
- 延迟 **20/10/5 ms** 来自 M 波–H 波间隙（通用神经延迟），**比猫 Ib 30–40 ms、人 ~64 ms 都短**；
- 引用恰为 [35] Pearson 1993、[28] Prochazka 1997、[36] Grey 2007、[62] Sinkjaer 2000、[61] Donelan 2009。

| 生理成分 | Geyer F+ 是否包含 |
|---|---|
| GTO 力转导动态（lead-lag Ib(s)） | 否（直接反馈力） |
| 静息 Ib 自体抑制 / 符号开关 | 否（无状态依赖） |
| 双突触 vs 长潜伏期节律网络通路 | 否（塌缩为单一延迟） |
| 相位/激活门控 | 部分（构造上仅支撑相） |
| 长度依赖自动降增益 | **隐式有**（来自 Hill 肌肉力–长度曲线） |
| 异名扩散到协同肌 | 否（各肌用自身力） |
| 长 Ib 延迟 | 否（用 20/10/5 ms） |

> **诚实表述**：F+ 作为**控制律**有据且稳定机制正是 Prochazka 1997 提出、Donelan 2009 实测的那个；
> 作为**特定解剖 Ib 回路**的主张则是抽象——Geyer 自己的措辞是“observed in cats … **suggested** in humans”。
> 另：Geyer 的 **F−**（SOL 力→抑制 TA）是**另一个现象学构造**，也不是经典 Ib 自体抑制通路。

### 3.7 待补 / 下一步可追

- **Pearson & Collins 1993 全文**（最重要缺失，需机构/馆际）；
- **Bonasera & Nichols 1994**（跨关节 heterogenic 反射的关键）；
- Gossard 1994、Guertin 1995、Angel 1996（双突触 vs 长潜伏期 Ib 兴奋的一手来源）；
- Grey 2004 / Mazzaro 2005（排除皮肤/Ia/II 的对照）；Prochazka 1999（Ia/II/Ib 传递函数来源）；
- Geyer Seyfarth Blickhan 2003（F+ 控制律与初值的真正来源）。

---

## 4. 交互抑制证据（Crone / Petersen / Christensen / Tanaka / Shindo / Pierrot-Deseilligny）

> 完整逐字引文、方法、页码与检索日志见 [`references/_review_reciprocal_inhibition.md`](../references/_review_reciprocal_inhibition.md)。
> **六篇全部拿到全文**（三篇 J Physiol 为 PMC 绿色 OA；三篇 Exp Brain Res 经机构订阅）。
> 本节提炼对 `β(lCEN)` **负成分**与 Geyer `F−`/`L−` 的存在性/符号/相位先验。

### 4.0 访问状态与本地文件

| 文献 | 访问 | 本地文件（`references/`） |
|---|---|---|
| Crone 1987 | 全文（PMC1192076，扫描 OCR） | `Crone-1987-...pdf` |
| Petersen 1999 | 全文（PMC2269590） | `Petersen-1999-...pdf` / `-fulltext-...txt` |
| Christensen 2001 | 全文（PMC2278473） | `Christensen-2001-...pdf` / `-fulltext-...txt` |
| Tanaka 1974 | 全文（机构订阅；PubMed 无摘要） | `Tanaka-1974-...pdf` |
| Shindo 1984 | 全文（机构订阅） | `Shindo-1984-...pdf` |
| Pierrot-Deseilligny 1981 | 全文（机构订阅） | `Pierrot-Deseilligny-1981-...pdf` |

### 4.1 经典 Ia 交互抑制通路：双突触、甘氨酸能、与单突触兴奋不同

- 回路：**激动肌肌梭 Ia → Ia 抑制性中间神经元 → 拮抗肌 α 运动神经元**，**双突触抑制**（一个中间神经元，甘氨酸能）。
- **比单突触 Ia 兴奋多 ~0.9–1.0 ms 中枢延迟**（Crone 1987, p179，三受试者）→ “one interneurone is interposed”。
- **同一跟腱叩击**既**易化** Sol H 反射（单突触）又**抑制** TA H 反射（多 ~1 ms）；弱叩击 + 弱电刺激
  **相互易化**（10%+5%→24%，而非 15%）→ 汇聚于同一中间神经元通路（Crone 1987, p169）。
- **甘氨酸能**：甘氨酸受体 α1 亚基突变患者既无 H 反射也无 EMG 抑制（Petersen 1999, p615）。
- **符号规则**（Pierrot-Deseilligny 1981, p349）：**Ia 兴奋协同肌 / 抑制拮抗肌；Ib 相反**。
- 回路特征：**对侧 Ia 中间神经元相互抑制**、**Renshaw 细胞抑制 Ia 中间神经元**、**下行 α–γ 联动**
  （Tanaka 1974；Shindo 1984；Petersen 1999）。

> **建模含义**：抑制项应比单突触长度反馈**多 ~1 ms 中枢延迟 + 拮抗肌传入传导延迟**；
> 它是**脊髓、甘氨酸能、符号反转**连接，不是皮层的。

### 4.2 踝部肌对：两个方向都存在，但不对称

| 方向 | 神经配对 | 静息状态 | 强度 |
|---|---|---|---|
| 背屈肌 → 跖屈肌 | CPN（TA+腓骨肌）→ Sol | 多数有，但**变异大** | Crone 均值 **14.9%**（0–40%, n=60）；Shindo 4/5，max 14.1%；Pierrot **2/9**；Tanaka 1/6 |
| 跖屈肌 → 背屈肌 | PTN（Sol+腓肠肌）→ TA | 能引出 TA H 反射时**很强** | TA H 反射静息仅 **9/60** 可引出；Pierrot 表 1：**Sol→TA `−`，GM→TA `−−`** |

- **不对称明确**（Crone 1987, pp179–180）：**跖屈肌→TA 抑制强于 TA→Sol**（与猫“屈肌 Ia IPSP 大于伸肌”一致）。
- **GM→TA（`−−`）强于 Sol→TA（`−`）** → 模型把交互抑制项**从 GAS 引出至少和从 SOL 引出一样有据**。
- 腓骨肌：CPN 是混合神经（TA+腓骨肌），无纯腓骨肌→Sol 的分离证据。
- 人踝拮抗肌 Ia 效应**弱、受试者变异大**；真正强的是**踝→股四头肌**群 I 投射（Pierrot 1981, p349）。

> **对 `β(lCEN)` 负成分**：`SOL/GAS → TA` 抑制项有据且应**更强**；`TA → SOL` 也有据但**更弱、更变异**。

### 4.3 步行相位调制（Petersen 1999）——本组最硬的证据

| 通路 | 摆动相 | 支撑相 |
|---|---|---|
| **背屈肌 → 跖屈肌**（CPN→Sol） | **大**（摆动 onset 显著） | **缺如**（H 反射 0/4+0/5；EMG 抑制减小，早支撑最明显） |
| **跖屈肌 → 背屈肌**（PTN→TA） | **小**（显著低于匹配 EMG 的 tonic 背屈，n=10） | **未测**（任何受试者支撑相都引不出 TA H 反射） |

- **任务特异**（匹配背景 EMG），**非外周反馈**（缺血阻滞大直径传入后差异仍在）；
  **突触前抑制**是候选机制（未证实）；对侧 Ia 中间神经元相互抑制、Renshaw 也被讨论。
- **⚠️ 支撑相“跖屈肌→背屈肌大抑制”是推断，不是测量**：作者明言
  "**but we have at present no direct evidence of this**"（p617）。
- **刺激强度会饱和通路**（>~1.2–1.3×MT），解释了 Capaday 1990“支撑相与 tonic 相同”的表面矛盾。

> **建模含义**：交互抑制项需**相位依赖增益**。最有据的调度：**背屈肌→跖屈肌 摆动 ON / 支撑 OFF**；
> 跖屈肌→背屈肌 摆动 OFF（小）；其支撑相值是**假设**。

### 4.4 ⚠️ Christensen 2001：Geyer 的描述方向反了

- **实际发现**：**支撑相**（TA 静息）出现**单个很大**的拉伸反应（**全部 17 人**，onset 76±5 / 92±4 ms，
  peak ~110 ms，group I 驱动、部分经皮层）；**摆动相**（TA 激活）只有**小而变异**的反应。
  （摘要 p545、结果 p546、讨论 p554；Zuur 2009 独立佐证。）
- **Geyer [37] 的原话**“a large TA stretch response is present in **swing**, but suppressed … in **stance**”
  **与原文相反**。
- **但**：短潜伏期 **M1/Ia 成分（41±5 ms）**确实**摆动相有、支撑相无**；Christensen 明确把支撑相缺如归因于
  **双突触交互 Ia 抑制**（引 Petersen 1999，p554）。
- **正确表述**：把 TA 正长度反馈 `L+` 门控到摆动相**仍可辩护**，但依据是
  “**短潜伏期单突触 Ia 成分在支撑相被抑制**”，**不是**“TA 拉伸反应摆动大、支撑被抑制”。
- 模型**无法重现**支撑相那个大晚反应（onset ≥75 ms、部分经皮层）；实测拉伸反应有**速度阈 ~100 deg/s**、
  >~250 deg/s 饱和——纯线性长度反馈没有这种阈值。

### 4.5 “负长度反馈”是否是正确抽象？——传入对、增益错

- **传入支明确是 Ia（拉伸/肌梭）**：Crone 1987（p169）**纯跟腱叩击**（阈下、零激动肌激活）就抑制 TA H 反射
  → **纯激活驱动的抑制是错的**（丢失被动拉伸证据）。
- **但增益与激活相关**：γ 驱动使 Ia 放电 ≠ 长度（Vallbo）；**α–γ 联动**（Tanaka/Shindo）；
  下行门控**先于 EMG 50 ms**（Crone）；突触前抑制；收缩平行分级（Shindo）。
  **注意**：Crone **未能重复** Tanaka/Shindo 的 tonic 背屈增强，且认为 Ia 中间神经元**不是**收缩平行
  “自然”交互抑制的主源（主源是长潜伏期 group I 通路，p183）。
- **结论**：负长度反馈是**正确的传入抽象、错误的增益抽象**。
  **推荐**：**Ia 样长度/速度负反馈 + 相位依赖增益**，可加小量激活依赖项近似 α–γ；
  加 ~1 ms 额外延迟；两方向**耦合**（对侧 Ia 中间神经元相互抑制），不要设成彼此独立。

### 4.6 对 Geyer `F−` / `L−` 与 `β(lCEN)` 负成分的先验

| 模型项 | 通路存在性 | 方向/符号 | 相位调制 | 信号模态 | 总评 |
|---|---|---|---|---|---|
| **SOL/GAS→TA 抑制**（Geyer `F−`） | **有据**（Crone；Pierrot） | 正确，强于反向 | 支撑 ON 为**推断**，摆动 OFF 为**实测** | **错配**：模型用**力**（Ib 样），真实通路是 **Ia（长度）**；且 Ib 符号实际相反（GM→TA Ib `+`） | **保留**，但诚实重述、考虑改用 Ia 样输入 |
| **TA→Sol 抑制**（反向） | **有据**，弱/变异 | 正确 | **摆动 ON / 支撑 OFF 实测**（Petersen） | Ia 样 | **保留** + 相位增益（全组最有据的调制） |
| **HAM→HFL 抑制**（Geyer `L−`） | **人体无据**；猫同源 plausible；池未采样 | 假设 | 摆动特异，假设 | 长度（Ia 样，模态至少一致） | **作为显式建模假设保留**，或在 CMA-ES 中**消融**检验必要性 |
| **TA `L+` 门控到摆动** | 仅 **M1/Ia 成分**有据（Christensen） | 正确 | M1 摆动 ON / 支撑 OFF | 长度（Ia，一致） | **保留**，但正确引用；大晚支撑反应未被捕获 |

**两个关键纠正**：
1. **Geyer `F−` 模态错配**：它用 **SOL 力**驱动，却引用 **Ia** 通路；而经典群 I 模式中
   **GM→TA 的 Ib 是易化（`+`）**、Sol→TA Ib 无显著效应。所以“跖屈肌**力**抑制 TA”**既不是已证明的 Ib 通路，
   也不是字面 Ia 通路**，而是因为支撑相 SOL 力与 SOL 激活/长度相关而产生的**功能性抽象**。
2. **HAM→HFL `L−` 是显式假设**：六篇 + Meunier 1993 **均未采样髋屈肌运动神经元**；Geyer 的增益/偏移来自
   **双摆力学**（Appendix I），非神经生理辨识。**未被证伪，只是未测**——是 CMA-ES 消融的好候选。

### 4.7 与 §2 mask 问题的衔接：Pierrot 1981 符号图

Pierrot-Deseilligny 1981 表 1（p348）给出一张**符号图**，直接支撑 §2 的“**禁正留负**”：

| 条件→测试 | Ia | Ib |
|---|---|---|
| Sol→TA | **−** | NS |
| GM→TA | **−−** | **+** |
| CPN→Sol | **−** | NS |
| GM→Sol | **NS（缺如）** | −− |
| GM→Bi | +++ | − |
| Sol→Q / GM→Q | +++ | （易化） |

- **拮抗肌对（Sol/GM→TA、CPN→Sol）**：**Ia 抑制（负）有据**，而 Meunier 的**单突触 Ia 兴奋（正）缺如/弱**
  （GM→TA 0/6、InfSol→TA 0/6、DP→Sol 0/4）。→ `β(lCEN_{gastroc,soleus}, tib_ant)`、
  `β(lCEN_{tib_ant}, soleus)` 应**禁正留负**（§2.3），**不是结构零**。
- **GM→Sol Ia NS** 与 Meunier GM→Sol 0/9 **互相印证** → §2 最强剪枝点。
- **Ib 符号与 Ia 相反**（GM→TA Ib `+`）→ 再次说明 `β(fATN)` 与 `β(lCEN)` 的符号先验**不能混用**。

> ⚠️ Pierrot 1981 表 1 来自旧排版扫描的 `-layout` 提取，Ia 关键项（Sol→TA `−`、GM→TA `−−`、
> CPN→Sol `−`、GM→Sol NS）清晰，但**定量使用前应对照 PDF 图像复核**。

### 4.8 待补 / 下一步可追

- **Zuur 2009**（J Physiol 587:1669–1676，仅摘要）：早支撑 TA 大拉伸反射被 rTMS 抑制——佐证 §4.4。
- 髋屈肌运动神经元的交互抑制证据（HAM→HFL）——本组与 Meunier 均缺如，需专门文献。
- Crone 1987 提到的“长潜伏期 group I 抑制”一手来源（收缩平行自然抑制的主源）。

---

## 5. 综合：三类证据分别约束 β 的哪个子矩阵

| β 子矩阵 | 生理通路 | 主要证据来源 | 存在性 / 符号先验用法 |
|---|---|---|---|
| `β(lCEN)` **正**（j≠m） | 异名 Ia 单突触**兴奋** | **Meunier 1993** | 阳性→保留；阴性→**禁正**（符号约束/强先验），§2 |
| `β(lCEN)` **负**（j≠m） | Ia **交互抑制**（双突触、甘氨酸能） | **Crone/Petersen/Pierrot**（§4） | 拮抗肌对→**允许负**；跖屈肌→TA 强于反向；**需相位依赖增益**；**不由 Meunier 定** |
| `β(lCEN)` 对角 | 同名 Ia stretch reflex（M1） | Meunier + Christensen（§4.4） | 一般保留（正）；TA 的 M1 成分**摆动有/支撑无** |
| `β(fATN)` | Ib 力反馈（含 F+ 反转） | **Prochazka/Donelan/Grey/Pearson**（§3） | **同源 SOL/GAS 正力反馈最有据**；VAS 为控制原理外推；**跨关节踝→膝 Ib 无证据**；符号状态依赖反转 |
| `β(关节量)` | 关节感受器 + 跨关节状态 | 功能性，证据较弱 | 不宜过度归因到单一感受器 |

**一句话总结**：

> Meunier 1993 是“**异名 Ia 兴奋性长度反馈子矩阵**”的存在性/符号先验，不是增益表，也不覆盖抑制与力反馈。
> 用它的**阴性结果剪“兴奋性洞”是合理的降维**，但应实现为**禁正留负的符号约束**而非一刀切结构零；
> **抑制子矩阵**交给交互抑制文献（双突触、相位调制、跖屈肌→TA 更强），
> **力反馈子矩阵**交给 Ib 文献（同源 SOL/GAS 最有据，跨关节无证据）。

### 5.1 关键纠正清单（写论文/引用时务必避免的误述）

1. **不要写“Ib 力反馈贡献 SOL 支撑相 EMG 约 50%”。** Sinkjaer 2000 的 ~50% 是**传入总贡献**
   （非 Ia/非交互/非皮肤），作者倾向 group II，且反应随支撑相力增大而**减小**（与 Ib 预测相反）。
   Ib 特异数值是**猫**的 ~30–39%（Donelan 2009）。力相关性的人体证据是 **Grey 2007（仅 SOL、晚支撑相、间接）**。
2. **不要照抄 Geyer 对 Christensen 2001 的描述。** 原文是**支撑相大、摆动相小**；Geyer 写反了。
   可辩护的是“**短潜伏期 M1/Ia 成分**摆动有、支撑无（被交互抑制压低）”。
3. **Geyer `F−`（SOL 力→抑制 TA）模态错配。** 真实交互抑制是 **Ia（长度）**驱动；经典 **Ib 符号相反**
   （GM→TA Ib `+`）。力驱动抑制 TA 是功能性抽象，不是已证明的群 I 通路。
4. **HAM→HFL `L−` 是显式假设**，人体无直接证据（髋屈肌池未被采样），增益来自双摆力学。宜作 CMA-ES 消融候选。
5. **跨关节踝伸肌→膝 Ib 通路**：本组 Ib 文献**不支持**；异名 Ib 仅在**踝伸肌协同肌群内**（跖肌→MG，猫）。
   若要声称需 Bonasera & Nichols 1994 等，并标注外推。
6. **交互抑制支撑相“跖屈肌→背屈肌大抑制”是推断**（Petersen 1999 明言无直接证据）；
   实测的是“背屈肌→跖屈肌 摆动大/支撑缺如”和“跖屈肌→背屈肌 摆动小”。
7. **Pearson & Collins 1993 仅摘要可得**（付费墙）；涉及其细节须标注“摘要”或二手佐证。

---

## 6. 参考文献

**本仓库已核对**

- Meunier S, Pierrot-Deseilligny E, Simonetta M. Pattern of monosynaptic heteronymous Ia
  connections in the human lower limb. *Exp Brain Res.* 1993;96:534–544.
  （`references/BF00234121.pdf`）
- Geyer H, Herr H. A muscle-reflex model that encodes principles of legged mechanics produces
  human walking dynamics and muscle activities. *IEEE TNSRE.* 2010;18(3):263–273.
  （`references/Geyer-2010-...pdf`）

**Ib / 正力反馈（§3，DOI/PMCID 已对下载原文核对）**

- Prochazka A, Gillard D, Bennett DJ. Positive force feedback control of muscles.
  *J Neurophysiol.* 1997;77(6):3226–3236. DOI `10.1152/jn.1997.77.6.3226`. **全文+PDF**。
- Grey MJ, Nielsen JB, Mazzaro N, Sinkjaer T. Positive force feedback in human walking.
  *J Physiol.* 2007;581(Pt 1):99–105. DOI `10.1113/jphysiol.2007.130088`. PMC2075215. **全文（HTML，无 PDF）**。
- Sinkjaer T, Andersen JB, Ladouceur M, Christensen LO, Nielsen JB. Major role for sensory
  feedback in soleus EMG activity in the stance phase of walking in man.
  *J Physiol.* 2000;523(Pt 3):817–827. DOI `10.1111/j.1469-7793.2000.00817.x`. PMC2269822. **全文（HTML，无 PDF）**。
- Pearson KG, Collins DF. Reversal of the influence of group Ib afferents from plantaris on
  activity in medial gastrocnemius muscle during locomotor activity.
  *J Neurophysiol.* 1993;70(3):1009–1017. DOI `10.1152/jn.1993.70.3.1009`. PMID 8229157. **仅摘要（付费墙）**。
- Donelan JM, McVea DA, Pearson KG. Force regulation of ankle extensor muscle activity in
  freely walking cats. *J Neurophysiol.* 2009;101(1):360–371. DOI `10.1152/jn.90918.2008`. **全文+PDF**。
- Prochazka A, Gritsenko V, Yakovenko S. Sensory control of locomotion: reflexes versus
  higher-level control. *Adv Exp Med Biol.* 2002;508:357–367. DOI `10.1007/978-1-4615-0713-0_41`. **全文+PDF**。

**交互抑制（§4，DOI/PMCID 已对下载原文核对；六篇均全文）**

- Crone C, Hultborn H, Jespersen B, Nielsen J. Reciprocal Ia inhibition between ankle flexors
  and extensors in man. *J Physiol.* 1987;389:163–185. DOI `10.1113/jphysiol.1987.sp016652`. PMC1192076.
- Petersen N, Morita H, Nielsen J. Modulation of reciprocal inhibition between ankle extensors
  and flexors during walking in man. *J Physiol.* 1999;520(Pt 2):605–619.
  DOI `10.1111/j.1469-7793.1999.00605.x`. PMC2269590.
- Christensen LO, Andersen JB, Sinkjaer T, Nielsen J. Transcranial magnetic stimulation and
  stretch reflexes in the tibialis anterior muscle during human walking.
  *J Physiol.* 2001;531(Pt 2):545–557. DOI `10.1111/j.1469-7793.2001.0545i.x`. PMC2278473.
- Tanaka R. Reciprocal Ia inhibition during voluntary movements in man.
  *Exp Brain Res.* 1974;21:529–540. DOI `10.1007/BF00237171`.（PubMed 无摘要）
- Shindo M, Harayama H, Kondo K, Yanagisawa N, Tanaka R. Changes in reciprocal Ia inhibition
  during voluntary contraction in man. *Exp Brain Res.* 1984;53:400–408. DOI `10.1007/BF00238170`.
- Pierrot-Deseilligny E, Morin C, Bergego C, Tankov N. Pattern of group I fibre projections from
  ankle flexor and extensor muscles in man. *Exp Brain Res.* 1981;42:337–350. DOI `10.1007/BF00237499`.

**佐证 / 待追（仅摘要或未下载）**

- Zuur AT, Christensen MS, Sinkjaer T, Grey MJ. Tibialis anterior stretch reflex in early stance
  is suppressed by repetitive transcranial magnetic stimulation. *J Physiol.* 2009;587:1669–1676.
  DOI `10.1113/jphysiol.2009.169367`.（仅摘要）
- Bonasera & Nichols 1994, *J Neurophysiol* 71:1096–1110（跨关节 heterogenic 反射，未下载）。
- Geyer H, Seyfarth A, Blickhan R. Positive force feedback in bouncing gaits?
  *Proc R Soc B.* 2003;270:2173–2183.（Geyer F+ 控制律与初值的真正来源，未下载）。

> 卷期/页码/DOI 已对下载原文核对；精确逐字引文见 `references/_review_Ib_force_feedback.md`
> 与 `references/_review_reciprocal_inhibition.md`。Pearson & Collins 1993 全文需机构/馆际获取。