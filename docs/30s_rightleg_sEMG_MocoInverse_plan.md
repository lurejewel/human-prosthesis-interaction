# 30 s（210–240）右腿 sEMG Tracking（per-muscle s，跨段共享）计划

> 只对 **右腿 5 块肌肉** 做 sEMG tracking；  
> **每块肌肉一个 s，整个 30 s 的所有段共享**；  
> 不是每段一套 s，也不是 5 块肌肉共用一个标量；分段并行。

## 0. 先明确三个约束

### 0.1 “只做右腿肌肉”不等于删掉左腿

建议保留左腿肌肉作为模型的 actuators，只是**不给左腿加 sEMG tracking**：

- 右腿 5 块：`rect_fem_r, vasti_r, hamstrings_r, tib_ant_r, gastroc_r`；
- 左腿 5 块仍然存在、仍可被 effort 优化；
- 如果把左腿肌肉从模型里删掉，左腿没有主动驱动，MocoInverse 的动力学通常不可行。

准确说法：

> **只对右腿 5 块加 EMG tracking；左腿肌肉保留但不加 EMG 约束。**

### 0.2 s：每块右腿肌肉一个值，所有段共享

采用：

\[
s_i \quad (i=\text{rect\_fem\_r, vasti\_r, hamstrings\_r, tib\_ant\_r, gastroc\_r})
\]

- 5 个 s；
- 每个 s 在整个 210–240 s 内不变；
- **不允许每段一个 s**；
- 这 5 个 s 是右腿的全局校准，不是 30 s 内随时间变化的参数。

### 0.3 右腿 tracking 不是“计算量减半”

只跟踪 5 块而不是 10 块，减少的是 tracking 项的 control 数和对应导数，不是整个 NLP：

- 动力学、states、18 块肌肉的 effort、reserves 仍然存在；
- 因此不会严格快一半；预计 tracking 部分省 20–35%，加上窗口从 70 s 减到 30 s（段数约 43%）和并行，整体 wall time 可降到原来的约 1/4–1/3。

## 1. 窗口与数据

### 1.1 窗口

```text
210–240 s，共 30 s
```

确认 IK/GRF/model 覆盖 210–240，且没有缺帧/NaN。

### 1.2 右腿 sEMG reference

从当前处理后的 sEMG 生成：

```text
experiments/moco/outputs/emg_excitation_210_240_right.sto
```

只保留 5 列：

```text
/forceset/rect_fem_r
/forceset/vasti_r
/forceset/hamstrings_r
/forceset/tib_ant_r
/forceset/gastroc_r
```

要求：

- 时间覆盖 210–240；
- 幅值沿用当前 0–1 归一化；
- 10 个原始通道中只取右腿 5 个；
- 检查窗口边缘 NaN；必要时裁到有效区间。

### 1.3 左腿

左腿肌肉保留在模型里，但不设 reference label、不设 tracking weight。
它们仍受 effort 和动力学约束，可以自由优化。

## 2. 分段与并行

### 2.1 分段

沿用已验证配置：

```text
retain segment = 2.0 s
overlap = 0.5 s（两侧 burn-in）
solve window = 3.0 s
mesh = 0.02 s（外层）
```

30 s 约 **15–16 段**（4 个 chunk 各 4 段）。

每段输出：

- solution；
- summary；
- activation report；
- residuals；
- 该段实际使用的 5 个 s（全局同一组）。

### 2.2 并行结构

推荐 **按窗口 chunk 并行**：

```text
210–217.5
217.5–225
225–232.5
232.5–240
```

4 个 chunk，每个约 7.5 s、4 段，4 个进程并行。
每个 chunk 内部仍按 retain 2 s / solve 3 s 分段。

建议 4 workers；如果每进程内存过大，降到 2–3。

### 2.3 Agent 层面

Agent 层面只适合并行做：

- 结果 QC；
- 画图；
- 逐段异常排查；
- 文献/参数核对。

不要让多个 agent 同时写同一个 solution/CSV；重计算仍用本地进程并行。

## 3. s 的更新

### 3.1 每块右腿肌肉一个 s

固定轨迹时，对每个肌肉 \(i\)：

\[
s_i=\mathrm{clip}\left(
\frac{\sum_{k}\int e_i^{(k)}u_i^{(k)}}
{\sum_{k}\int (u_i^{(k)})^2},
0.5,2
\right)
\]

其中 \(k\) 是所有 retain 的 2 s 段。

- 每个肌肉独立更新；
- 所有段共享同一个 \(s_i\)；
- 不使用 overlap/burn-in；
- 当前 tracking 权重为 1；若以后给通道不同权重，要写进分子/分母。

### 3.2 初始值

从 1.5 s 最终 s13 取右腿 5 个值：

```text
rect_fem_r   1.349567
vasti_r      1.623255
hamstrings_r 0.744262
tib_ant_r    1.107898
gastroc_r    0.862097
```

第一版直接用这 5 个值作为 30 s 初值；另外做一次全部 s=1 的初值对照。

### 3.3 外层停止

每轮记录：

- \(J\)、\(r_J\)；
- \(R_s=\max_i|\hat s_i-s_i|/s_i\)（5 个右腿肌肉的最大值）；
- \(G_s=J(z,s)-J(z,s_{\rm next})\)；
- \(R_u\)：相邻轮 controls/activation 相对变化；
- 段级 seam。

建议：

- \(r_J<0.1\%\) 连续两轮；
- \(R_s<1\%\)；
- \(G_s\) 接近重复求解噪声；
- \(R_u<1–2\%\)。

5 个 s 比 1 个 s 更容易出现弱辨识；如果 J/轨迹稳定但某个 \(s_i\) 仍漂移，
应报告该通道的弱辨识区间，而不是无限追小数位。

## 4. 代码改动

当前 `moco_inverse.py` 默认 10 块肌肉、`--emg-s-values` 要求 10 个值。需要加：

### 4.1 `--emg-leg right`

- 定义 `EMG_MUSCLES_R`（5 块）；
- `--emg-tracking` 的 reference labels/weights 只加右腿；
- effort 仍覆盖全部 18 块肌肉；
- `activation_report`/`residuals` 仍可包含全部 18 块，但 tracking report 只列 5 块。

### 4.2 `--emg-s-values` 支持 5 个值

- `--emg-leg right` + 5 个 s；
- 解析器按 active muscle 数校验；
- `moco_inverse_s_used_<tag>.csv` 写 5 行（每个右腿肌肉一行）；
- 不写 1 行全局 s。

### 4.3 `read_moco_results`

- 已支持 `moco_inverse_s_used_<tag>.csv`；
- 5 行 s 文件时，未跟踪肌肉为 NaN，图 1 只画右腿 5 条绿线；
- 标题显示 `s=...`，不再显示 g；
- 右腿 5 个 s 分别来自 CSV，不能共用一个值。

### 4.4 并行 wrapper

新增：

```text
experiments/moco/run_30s_right_parallel.py
```

职责：

1. 生成 4 个 chunk 的 `--start/--end`；
2. 每个 chunk 调用 `moco_inverse.py`，共享同一组 5 个 s；
3. 等待 4 个进程；
4. 合并 summary、solution、activation report；
5. 计算 seam；
6. 写 `s_next`（5 个值）。

## 5. 运行模板

### 5.1 第一个 chunk（smoke）

```powershell
python experiments\moco\moco_inverse.py `
  --start 210 --end 217.5 --segment 2.0 --overlap 0.5 --mesh 0.02 `
  --sagittal-cost --nonsagittal-weight 1e-3 `
  --emg-tracking experiments\moco\outputs\emg_excitation_210_240_right.sto `
  --emg-weight 20 --emg-leg right `
  --emg-native-b-fixed `
  --emg-s-values "1.349567,1.623255,0.744262,1.107898,0.862097" `
  --emg-sparsity-random --tol 1e-4 `
  --tag B30_right_s13_chunk0
```

### 5.2 4 个 chunk 并行

```text
chunk0 210–217.5
chunk1 217.5–225
chunk2 225–232.5
chunk3 232.5–240
```

### 5.3 外层

```text
round 0: s = s13 right 5 values
round 1: s = s_next
round 2: s = s_next
round 3: s = s_next
```

每轮用上一轮同 chunk/同段的 solution 作 guess（`setGuess` 已修复，原生路径也生效）。

## 6. 计算量估计

30 s ≈ 15–16 段；单段 retain 2 s / solve 3 s，约 3–6 min（0.02 mesh，5 个 s，无 guess）。

- 串行一轮：约 **60–90 min**；
- 4 chunk 并行：约 **20–35 min/轮**（实测 4 段/chunk 约 19–34 min）；
- 3 轮：约 **60–100 min**；
- 0.01 mesh 最终筛查：约 2–3 倍单轮，视段数选择性跑。

这比 70 s、10 s 的方案快很多；但“只右腿”不会把动力学/effort 减半，实际加速主要来自
窗口缩短和并行。

## 7. 验收

### 7.1 段级

- 每个 chunk/segment `solved=True`；
- 无 max_iterations；
- 约束/对偶残差在 tol 内；
- 段边界 seam 不恶化。

### 7.2 Tracking

右腿 5 块：

- cosθ；
- Pearson r；
- 同空间 nRMSE；
- excitation/activation 曲线；
- 峰值、饱和时间比例。

B 比较 `e` vs `s_i u`（EMG 空间）或 `u` vs `e/s_i`（control 空间）；A 比较 `e` vs `u`。

### 7.3 饱和与 reserves

- 右腿 5 块是否 saturation；
- 左腿/非跟踪肌肉的 saturation 单独列；
- reserves 峰值/持续时间；
- 与 A 基线比较。

### 7.4 稳健性

- 网格 0.02→0.01（固定最终 5 个 s，选 2–4 段）；
- 初值全部 s=1 的 30 s 外层 2 轮；
- 5 个 s 与 single-s 的对照只作为额外敏感性，不作为主方案。

## 8. 交付物

```text
experiments/moco/outputs/
  moco_inverse_solution_B30_right_s13.sto
  moco_inverse_summary_B30_right_s13.txt
  moco_inverse_s_used_B30_right_s13.csv   # 5 行
  B30_right_tracking_report.csv
  B30_right_seam_report.csv
  B30_right_reserve_summary.csv
```

以及：

- A no-gain 30 s 右腿基线；
- 4 chunk 的逐段 summary；
- 3 轮外层 summary；
- 网格/初值稳健性；
- 单页生物力学解读。

## 9. 推荐执行顺序

1. **P0**：生成 210–240 右腿 5 列 sEMG reference；确认 IK/GRF 覆盖；
2. **P0**：chunk0 1 段冒烟；再 chunk3 1 段（检查窗口末端）；
3. **P1**：A no-gain 30 s 右腿基线；
4. **P1**：B 5-s（s13 初值）全 4 chunk；
5. **P1**：QC、seam、tracking、reserve；
6. **P2**：外层 3 轮（更新 5 个 s）；
7. **P2**：网格/初值稳健性；
8. **P3**：如果某个 s 弱辨识，报告其区间；不要退回每段一套 s。

## 10. 风险

| 风险 | 预案 |
|---|---|
| 左腿被误删导致不可行 | 保留左腿肌肉，只取消其 EMG tracking |
| 5 个 s 中某个弱辨识 | 报告区间；不强行固定或加先验，除非研究需要 |
| sEMG 各通道归一化不同 | 每块肌肉一个 s，分别更新；报告归一化方式 |
| 30 s 仍然慢 | 4 chunk 并行；0.02 mesh；只在最终做 0.01 |
| seam 差 | 检查 chunk/segment 边界和 guess 传递 |
| reserves 大 | 报告绝对值和权重设置，不把 EMG 拟合当生理验证 |

## 11. 执行结果（2026-10-29）

### 11.1 已完成

- 210–240 s 右腿 5 列 sEMG reference（实际覆盖 209.5–240，含 burn-in）；
- `--emg-leg right`：只跟踪右腿 5 块，左腿保留但无 EMG tracking；
- 4 chunk × 4 段并行，全部 `solved=True`；
- A no-gain 30 s 右腿基线；
- B fixed-s 初值 + 2 轮外层 s 更新；
- 合并后的 30 s 解、summary、s_used、tracking/outer 报告；
- `read_moco_results` 默认 `RightOnly=true`。

### 11.2 最终候选

```text
solution: moco_inverse_solution_B30_right_final.sto
summary : moco_inverse_summary_B30_right_final.txt
s_used  : moco_inverse_s_used_B30_right_final.csv
tracking: B30_right_tracking_report_final.csv
outer   : B30_right_outer_summary.csv
```

最终使用 s（5 个右腿肌肉，全 30 s 共享；第 7 轮实际使用值）：

```text
rect_fem_r   1.662868
vasti_r      1.493015
hamstrings_r 0.851994
tib_ant_r    1.177731
gastroc_r    1.068418
```

外层 **6 次有效更新**（r3 因脚本重启把 `S` 初始化回 s2，重复了 s2，已排除）：

| round | J | mean cos | mean r | mean nRMSE | R_s→next |
|---|---:|---:|---:|---:|---:|
| 0 | 204.830963 | 0.8904 | 0.7424 | 0.4054 | 3.97% |
| 1 | 204.490759 | 0.8906 | 0.7419 | 0.4036 | 3.80% |
| 2 | 204.111790 | 0.8906 | 0.7419 | 0.4021 | 3.67% |
| 4 | 203.792004 | 0.8907 | 0.7419 | 0.4005 | 3.57% |
| 5 | 203.515685 | 0.8910 | 0.7420 | 0.3988 | 3.47% |
| 6 | 203.244337 | 0.8911 | 0.7422 | 0.3972 | 3.37% |
| 7 | **203.028772** | **0.8914** | **0.7424** | **0.3956** | **3.27%** |

- 6 次更新总 J 下降 **0.88%**；
- 每轮 \(r_J\) 约 **0.11–0.19%**；
- tracking 指标几乎不变（cos +0.001，nRMSE −0.010）；
- \(R_s\) 从 3.97% 降到 **3.27%**，每轮约 **0.10%**；
- s 仍在缓慢漂移：rect_fem_r 1.35→1.66，vasti_r 1.62→1.49，gastroc_r 0.86→1.07。

**不能据此断言“再跑必然撞 bounds”**。目前只能确认：6 次更新尚未收敛，趋势与尺度补偿/慢收敛相符。

补充检查：

- 相邻轮 r6→r7：五块 tracked 肌肉原始 u 的 RMS 只变 **0.3–3%**，相关系数≈1，`s·u` 几乎不变；
- 累计 r0→r7：原始 u 变化更大（rect_fem_r −18%，gastroc_r −17%，hamstrings_r −10%，tib_ant_r −4%，vasti_r +4%），但 `s·u` 只变 1–4%，形状相关 >0.995；
- 右腿未 tracking 肌肉的 RMS 变化 ±3% 以内；reserves 峰值基本不变。

所以：**EMG 拟合稳定不等于原始 u 稳定；累计幅值确实变了，但主要是 s 在补偿尺度。**
最终是否内部收敛仍未判定，不能仅凭这 6 轮下结论。

因此结论是：

> **J 和 tracking 已稳定；s 没有内部稳定点，属于弱可辨识/边界漂移。**
> 不应把 7 轮的 s 当作唯一生理增益；也不应为了 \(R_s\) 再盲目跑 20 轮。

### 11.3 A vs B（30 s 右腿）

A = A no-gain（scale=1，无 g）；B = 上述 5-s B。

| 指标 | A | B |
|---|---:|---:|
| J_total | 208.303850 | **203.028772** |
| mean cos | 0.8698 | **0.8914** |
| mean Pearson r | 0.7034 | **0.7424** |
| mean nRMSE | 0.4219 | **0.3956** |

逐肌肉：

| muscle | A cos/r | B cos/r | B nRMSE |
|---|---:|---:|---:|
| rect_fem_r | 0.987 / 0.972 | **0.994 / 0.988** | 0.114 |
| vasti_r | 0.584 / 0.089 | **0.733 / 0.339** | 0.681 |
| hamstrings_r | **0.868 / 0.768** | 0.833 / 0.727 | 0.553 |
| tib_ant_r | 0.958 / 0.810 | **0.960 / 0.817** | 0.280 |
| gastroc_r | **0.953 / 0.878** | 0.937 / 0.841 | 0.350 |

- 最大改善是 **vasti_r**（r 0.089→0.339，cos 0.584→0.733）；
- rect_fem_r、tib_ant_r 也改善；
- hamstrings_r、gastroc_r 的 B 拟合反而比 A 差；
- 所以 B 的总体 mean 更好，但不是每块肌肉都更好。

### 11.4 seam 与饱和

chunk 边界 seam（activation max）：

```text
chunk0->1 0.0051
chunk1->2 0.0114
chunk2->3 0.0157
chunk3 end 0.0086
```

最大 **0.0157**，低于 0.02；可接受，但比 1.5 s 单窗口的 seam 略大。

activation >0.99：

```text
soleus_r    1.000   2.8% 时间
glut_max_l  1.000   8.4% 时间
```

这两块都**不在右腿 5 块 tracking 列表**；右腿 tracked 肌肉没有 activation 饱和。

### 11.5 reserves

最终 30 s 的 reserve 峰值（控制量）最大的包括：

```text
hip adduction   L 62.5 / R 51.9
hip rotation    L 33.9 / R 24.1
knee extension  L 29.1 / R 12.0
back bending    L 30.8
back extension  24.4
hip flexion     L 6.7
ankle dorsiflex R 6.9
```

左腿虽然不跟踪 EMG，但仍参与 effort/动力学，因此左侧 reserves 仍然存在，部分比右侧大。
不能把“右腿 tracking 好”解释成“全身肌肉力已验证”。

### 11.6 计算量

- 4 chunk 并行，单 chunk 4 段约 **19–34 min**；
- 一轮 30 s 约 **34 min wall**（4 进程并发）；
- 2–3 轮约 **1.5–2 h**；
- 总计算量约 **118 min CPU**（4 chunk 相加）。

### 11.7 还没做

- 网格 0.02→0.01 的最终筛查；
- 初值 s=1 的 30 s 对照；
- 完整 30 s serial 交叉验证（已停掉，避免重复）；
- 左腿/右腿 reserves 的独立肌骨验证。

当前结果应表述为：

> **右腿 5-s 的 30 s EMG-informed 轨迹；J 和 tracking 指标稳定，但 s 弱可辨识（R_s≈3.7%），reserves 仍大，不能声称参数唯一或肌肉力已验证。**
