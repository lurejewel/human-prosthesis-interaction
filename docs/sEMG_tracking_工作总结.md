# sEMG Tracking 工作总结

## 0. 一句话概览

目标是在 MocoInverse 里让肌肉激活跟踪 sEMG，并比较两种跟踪写法：原来的 A 形式
\((u-ge)^2\)，和新的反向缩放 B 形式 \((e-su)^2,\ g=1/s\)。
B 形式的数学和插件实现已经跑通，0.3 s 结果已经可以收尾；1.5 s 的 fixed-s B 形式也收敛了，
但 **1.5 s 的 free-s joint 求解还没有收敛**，目前卡在“s 方向太平面、联合优化病态”这一点上。

---

## 1. 要解决什么问题

### 1.1 A 形式的问题

原来的 `MocoControlTrackingGoal` 误差是：

\[
(u - g e)^2
\]

其中 \(u\) 是模型激活，\(e\) 是 sEMG 参考，\(g\) 是可优化的增益。这个写法的问题是：

- 如果 \(u\) 偏大，优化器可以把 \(g\) 调小，让误差看起来变小；
- 误差里 \(u\) 和 \(g\) 的尺度可以互相补偿，容易把 \(g\) 推到下界；
- 我们真正想要的是“模型激活和 EMG 形状是否一致”，而不是让一个自由增益把幅值差吃掉。

### 1.2 B 形式

把误差改成：

\[
(e - s u)^2,\qquad g = 1/s
\]

这里 \(s\) 是“模型幅值相对 EMG 的缩放”。固定 \(u\) 时，最优 \(s\) 有闭式解：

\[
s_i^* = \mathrm{clip}\left(
\frac{\int e_i u_i\,dt}{\int u_i^2\,dt},
0.5, 2
\right)
\]

B 形式的好处是：它把“幅值尺度”显式地放在 \(s\) 上，而不是让 \(g\) 和 \(u\) 一起漂。
但后来发现，B 形式也有自己的尺度简并：\(s\) 变大时 \(u\) 可以变小，effort 项会下降，
所以 \(s\) 会被 effort 偏置往 2 推。这一点在结果里必须说明。

需要注意：\(g=1/s\) 不意味着 A/B 是同一个目标的简单倒数参数化，因为

\[
(e-su)^2=s^2(u-ge)^2.
\]

B 形式同时改变了误差的尺度，所以 A/B 的 objective 不能直接比较，应主要比较 \(\cos\theta\)、形状和同一空间的 RMS。
另外，`MocoControlTrackingGoal` 跟踪的是 **controls**；对肌肉通常是 **excitation**。如果没有忽略 activation dynamics，
不应把 control 与 activation 状态混称，分析和绘图也要用同一种变量。

---

## 2. 实现路线

### 2.1 插件

复制 `MocoControlTrackingGoal` 的源码，改名为：

- `MocoReverseScaledControlTrackingGoal`
- `MocoReverseScaledControlTrackingGoalReference`

只改误差行：

```cpp
double error = refValue - (scaleFactor * control);
```

也就是 \((e-su)^2\)。其他部分（参考信号、control 权重、scale factor 注册）尽量保持原样。

插件文件位置：

```text
experiments/moco/revgoal/
  MocoReverseScaledControlTrackingGoal.h
  MocoReverseScaledControlTrackingGoal.cpp
  RegisterTypes_revgoal.cpp
  CMakeLists.txt
  build_vs/Release/osimRevGoal.dll
```

### 2.2 Python 接入

`experiments/moco/moco_inverse.py` 增加了：

- `--emg-rev-goal`
- `--emg-rev-plugin`
- `--emg-initial-guess`
- `--emg-s-values`（每个肌肉固定 s）
- `--emg-sparsity-none`
- `--emg-sparsity-random`
- `--emg-params-initsystem`

流程：

1. 先建一个普通 `MocoControlTrackingGoal`，加 reference labels 和 scale factors；
2. `printToXML()`；
3. 把 XML 里的 goal/helper 标签改成插件类名；
4. `LoadOpenSimLibrary()` 加载 DLL；
5. `MocoGoal.makeObjectFromFile()` 读回；
6. `problem.addGoal()`。

这里有一处需要独立验证：原生 `MocoControlTrackingGoal` 的 `control → scale factor name` 关联保存在
`m_scaleFactorMap`，它是普通 C++ 成员，不是可序列化 property。仅保存 `<MocoScaleFactor>` 再改 XML 类名读回，
理论上不足以保留这层映射。当前插件用 `gain_<shortName>` 的命名回退来恢复映射，且日志里能看到 `found=true`，
但这仍是项目特定的约定；C++ 源码里也因此加了回退逻辑。下一步应补一个非 \(s=1\) 的 XML 参数绑定测试。

### 2.3 V1：s=1 时与 A 形式精确等价

用 no-gain 的 A 形式和 B 形式（不加 scale factor，s 默认为 1）在同一个 0.3 s 窗口比较：

```text
A 形式 no-gain: objective = 1.2127594, 33 iterations
B 形式 no-gain: objective = 1.2127594, 33 iterations
```

误差公式和插件实现是**精确等价**的。这是整个工作的数学基线。

---

## 3. 遇到的难点与解决过程

### 3.1 插件编译、加载、XML roundtrip

难点：

- `OSIMMOCO_API`、导出宏、`not` 关键字、`/permissive-`、`/utf-8`；
- helper 类 `MocoReverseScaledControlTrackingGoalReference` 没有注册，XML 读回时报“unrecognized type”；
- `MocoGoal.h` 路径、spdlog 头文件、Simbody 路径；
- CMake 的 NMake/U1054、MSBuild 临时目录权限问题。

解决：

- 用 VS 18 BuildTools + CMake VS generator；
- 把 `TMP/TEMP` 指到 workspace 的 `.tmp_nmake`；
- 注册 goal 和 helper 两个类型；
- 加 `/permissive-` 和 `/utf-8`；
- 在 SDK include 里补 `spdlog/include`、Simbody include。

### 3.2 `addScaleFactor` 导致 native crash

在 OpenSim 4.6 里，给 goal 加 scale factor 后第一次求解会出现 `0xC0000005`。
原因是 solver 没有正确重建 problem/parameter 初始化。

解决：

```python
solver.resetProblem(problem)
solver.set_parameters_require_initsystem(False)
```

这两句现在是所有 B 形式求解的固定配置。

### 3.3 最重要的问题：s 完全不更新（H2）

#### 现象

把 s 的 bounds 设成 `[0.5, 2]` 后，求解成功，但 10 个 s **全部停在 1.25000001**。
1.25 正好是 bounds 的中点，也就是 scale factor 的默认初值。

这说明：参数存在、cost 也确实用了 s，但优化器没有更新它。

#### 诊断

做了几组对照：

1. **E1：插件内打印读到的 s**
   - 在 `calcIntegrandImpl` 里打印 s 和组件地址；
   - 在 `none` 模式下能看到 `1e-5` 量级的有限差分扰动，说明 s 的更新能传到插件读到的组件；
   - 首轮 177 次调用 s=0，之后才出现非零值。

2. **E2：关闭 sparsity detection**
   - 本工作流中 MocoInverse 设成 `optim_sparsity_detection="random"`、`fd="forward"`；
   - 改成 `none` 后，s 离开 1.25，变成 `0.5 / 2.0 / 1.32 / ...`；
   - 同一个 0.3 s 配置，`none` 能学 s，`random` 不能。

3. **E4：A 形式对照**
   - 同样的 random 设置下，A 形式的 g 能明显离开 1.25；
   - 所以不是“整个 solver 配置变了”，而是 **random/initial-guess sparsity 对这个自定义 goal 不可靠**。

4. **initial-guess 模式**
   - s 初值 1.0：只漂到 1.24968，objective 6.366；
   - s 初值 1.25：完全不动，objective 6.294；
   - s 初值 1.6：收敛到 1.24817，objective 6.59；
   - bounds 改成 [0.5, 3]：s 跑到 1.765，正好是新区间的中点。
   - 这些结果说明：`initial-guess` 下 s 基本只被 barrier 拉向 bounds 中心，目标梯度没有真正进入 s。

5. **fixed-s 对照**
   - 固定 s、用 `random`：J=5.0976，复现不了 joint-none 的 2.30；
   - 固定 s、用 `none`：J=2.2742，和 joint-none 的 2.3015 只差 1.2%。
   - 所以 **fixed-s 路径本身没坏，random sparsity 才是 5.0 的来源**。

#### 结论

- `random` 和 `initial-guess` 对这个自定义 goal 都不可靠；
- 只有 `optim_sparsity_detection="none"` 能让 s 真正参与优化；
- `parameters_require_initsystem=False` 是必须的；设成 `True` 会直接 crash；
- **H2 的源码机制已经找到**：该版本 `CasOCFunction.cpp` 的 `get_jac_sparsity()` 在探测某个输入块时，会把**其他输入块置零**。controls 和 parameters 是两个不同输入块，于是 B 形式的乘积项
  \[
  (e-su)^2
  \]
  在探测 s 时 \(u=0\)，变成 \(e^2\)；在探测 u 时 \(s=0\)，也变成 \(e^2\)。也就是说，**s 和 u 的依赖都会被“另一个块置零”遮掉**。
- 这解释了为什么 `random`/`initial-guess` 下 s 基本停在 bounds 中心，也解释了为什么 **fixed-s + random 也可能错**：fixed-s 只是把 s 的 bounds 设成 \([s,s]\)，参数块仍然存在，探测 controls 时仍可能被置零。
- 这是 **Moco/CasOC 核心层的稀疏结构问题**，不是只重编译 goal DLL 能修的；换新版本也不能假定已解决（main 里仍有同样写法）。


### 3.4 `none` 太慢，1.5 s joint 求解不收敛

`none` 下 0.3 s 能收敛，但 1.5 s 的 free-s joint 非常难。这里要说明：`none` 使用的是**保守块稀疏**，不是把整个 NLP 变成稠密矩阵；它慢是因为块更大、有限差分更多，而不是“10 个参数让全矩阵稠密”。

- A+none 1.5 s：2563 iter / 25.6 min 收敛；
- B+none 1.5 s：3000 iter / 28.4 min 仍失败（早期从零开始）；
- 后来用 1.5 s A 解 + 0.3 s 的 s 做 guess，跑到 2586 iter / 约 40 min，objective 在 16–20 之间振荡，inf_pr≈1e-3，inf_du≈0.5–4，没有 EXIT。

初始点诊断：

```text
iter 0: objective=16.7677, inf_pr=4.37, inf_du=0.416
iter 1: objective=18.0225, inf_pr=3.90, inf_du=3.55
```

初始 guess 本身不可行（inf_pr 4.37），所以不能简单说“求解器从 16.8 上坡到 18”。
一种可能的解释是：B 形式在 s 方向非常平，Hessian 用 limited-memory 近似，
s 和 u 的尺度又可以互相补偿，使联合问题弱辨识，在 restoration 和振荡之间来回；
但这只是提示，尚未证明 Hessian 近奇异。

#### 3.4.1 自重启对照：先修 guess 调用顺序

之前从自己的解重启仍出现 `iter 0 objective=80.09, inf_pr=297`，原因已定位：
**`--emg-initial-guess` 只在 `--emg-rev-goal` 分支里调用 `solver.setGuess()`；
原生 fixed-s / 外层路径根本没有设置 guess**，所以每次都是 bounds 冷启动。

把 `setGuess()` 移到 rev 分支外之后，用 s5 的同一完整解自重启：

| 设置 | iter 0 objective | iter 0 inf_pr | 最终 objective | iterations | wall |
|---|---:|---:|---:|---:|---:|
| 原生 random，修复前（未设 guess） | 80.09 | 297 | 10.369664 | 1446 | 264.8 s |
| 原生 random，修复后 | **10.2298** | **7.11** | **10.325376** | **228** | **53.9 s** |

修复后迭代数从 1446 降到 **228**，时间从 265 s 降到 **54 s**，目标基本不变。
说明 **primal guess 本身是有效的**；之前的 80/297 不是缺少对偶变量造成的，而是 guess 根本没进 solver。

仍然没有传入 `lam_x0/lam_g0`；iter 0 的 inf_pr=7.11 也不小（来自 `.sto` 重采样、贴界变量初始化和
IPOPT 的 bound push），所以**完整原始—对偶热启动仍可能进一步减少迭代**，但已经不是当前的主要瓶颈。

### 3.5 fixed-s 1.5 s：保底结果

固定 s 为 0.3 s 的 10 个值，用 1.5 s A 解做状态/控制初值，`none`，tol=1e-4：

```text
solved=True
iterations=1555
wall=904.6 s (15.1 min)
objective=10.653509
saturated=3/18
```

这是一个**有效的 1.5 s B 结果**，但 s 是给定的，不是辨识出来的。
它比 free-s joint 的 18 附近低很多，说明困难主要来自 s 的平坦方向，而不是 B 形式本身。

用这个 fixed-s 解的 u 反算 s*：

| muscle | s_0.3s | s*_1.5s | 漂移 |
|---|---:|---:|---:|
| rect_fem_r | 0.7602 | 0.8234 | **+0.063** |
| vasti_r | 2.0000 | 2.0000 | 0 |
| hamstrings_r | 0.5365 | 0.5807 | +0.044 |
| tib_ant_r | 1.4820 | 1.3782 | **−0.104** |
| gastroc_r | 0.5559 | 0.5939 | +0.038 |
| rect_fem_l | 2.0000 | 1.8511 | **−0.149** |
| vasti_l | 2.0000 | 1.6219 | **−0.378** |
| hamstrings_l | 2.0000 | 1.9063 | **−0.094** |
| tib_ant_l | 1.2057 | 1.2318 | +0.026 |
| gastroc_l | 0.5000 | 0.5537 | +0.054 |

**0.3 s 的 s 不能直接推广到 1.5 s**：vasti_l 漂移 0.38，rect_fem_l 0.15，
tib_ant_r 0.10，hamstrings_l 也从 2.0 回到 1.91。

### 3.6 s 的可辨识性很弱

在 0.3 s 的 1e-5 解上，用**求解器自己的 objective** 做固定-s profile：

| 通道/扰动 | objective | Δ vs base |
|---|---:|---:|
| base | 2.25959 | — |
| rect_fem_r −10% | 2.25990 | +3.1e-4 |
| rect_fem_r +10% | 2.25984 | +2.5e-4 |
| hamstrings_r −10%（到 0.5） | 2.25965 | +6e-5 |
| hamstrings_r +10% | 2.25969 | +1e-4 |

- rect 的 ±10% 只差 6e-5，hamstrings_r 只差几 e-5；
- 之前用梯形后处理看到的“不对称”主要是求积方式与 goal 不一致，用 solver objective 后基本对称；
- **s 方向非常平，s 是弱可辨识的**，报告里应给区间/说明，而不是单点“最优 s”。

同时，5 个通道在 0.3 s 贴界：

```text
vasti_r = 2.0
rect_fem_l = 2.0
vasti_l = 2.0
hamstrings_l = 2.0
gastroc_l = 0.5
```

这些数值是 bounds 决定的，不是数据单独识别的；报告里必须单独说明。

### 3.7 固定 s 时其实不需要插件（已交叉验证）

固定 \(s_i\) 时，B 形式可以严格改写成原生 `MocoControlTrackingGoal`：

\[
w_i(e_i-s_i u_i)^2
=
w_i s_i^2\left(u_i-\frac{e_i}{s_i}\right)^2.
\]

因此可以：

- 把 reference 改成 \(e_i/s_i\)；
- 把对应 control 的权重改成 \(w_i s_i^2\)；
- **不调用 `addScaleFactor()`，不保留任何代表 s 的 `MocoParameter`**。

这条路线完全绕开插件的参数绑定和 CasOC 参数稀疏问题。已经实际跑了两个交叉验证：

| 配置 | 插件 fixed-s B | 原生 fixed-s B | 差异 |
|---|---:|---:|---:|
| 0.3 s, 同一 s | **2.25959** | **2.26058** | +0.044% |
| 1.5 s, 同一 s | **10.653509** | **10.653837** | **+0.0031%** |

1.5 s 的 s* 也几乎逐位一致（例如 vasti_l：插件 1.621896，原生 1.621652；rect_fem_l：1.851134 vs 1.851100）。
原生路径没有参数，gains CSV 全是 `nan` 是预期的。需要注意：不能只改 reference 漏掉 \(s^2\) 权重；每轮从原始 reference 重建，不能累积缩放；多个 control 共用 reference 列时要复制成独立列。

**这说明：固定 s 时，B 形式的数学、XML 参数绑定和 cost 求值都是对的；剩下的问题在 free-s joint 求解，而不是插件实现。**

原生无参数路径还可以重新使用 `random`（不再有跨 controls/parameters 的遮蔽结构）：

| 配置 | none | random | random 耗时 |
|---|---:|---:|---:|
| 0.3 s | 2.26058 | 2.26049 | 122.9 s vs 296.4 s |
| 1.5 s 自重启 | 10.369963 | 10.369664 | 264.8 s vs 889.3 s |

差异只有 9e-5 / 3e-4，但 random 快 2.4–3.4 倍；因此外层固定-s 求解可以用 `random` 加速，free-s 插件路径仍然必须用 `none`。

### 3.8 热启动：当前状态与完整方案

当前 `solver.setGuess(trajectory)` 只把轨迹变量作为 `x0` 传给 CasADi/IPOPT。该版本
`CasOCTranscription.cpp` 里调用 `nlpFunc` 时只传了：

```text
x0, lbx, ubx, lbg, ubg
```

没有传 `lam_x0`、`lam_g0`，也没有设置 IPOPT 的 `warm_start_init_point=yes`。
不过之前“自重启 80/297、1292 次”的主要原因是 **`--emg-initial-guess` 只写在 rev 分支里，
原生路径根本没调用 `setGuess()`**；修复后同一 s5 自重启是 `iter 0 10.23/7.11`、228 次迭代、53.9 s。
所以现在的状态是：

- **primal warm start 已经可用**；
- **dual warm start 仍缺失**，但不再是当前最大瓶颈。

完整热启动可以分三层：

#### 方案 A：原始变量热启动（当前已经做到）

- `setGuess(trajectory)` 现在对原生 fixed-s / 外层路径也生效；
- 实测同一 s5 自重启：iter 0 从 80.09/297 降到 **10.23/7.11**，迭代从 1446 降到 **228**，时间从 265 s 降到 **53.9 s**；
- 如果还要进一步改善，可以在 IPOPT options 里设置：
  ```text
  warm_start_init_point = yes
  warm_start_bound_push = ...
  warm_start_mult_bound_push = ...
  ```
- 但 Moco 目前没有暴露任意 IPOPT options 的 property；需要改 `MocoCasADiSolver` 的
  `solverOptions`，或通过 `ipopt.opt` 文件做独立实验。
- 这仍然没有提供对偶变量，只能算“primal warm start”，不能保证从原解立刻收敛。

#### 方案 B：完整原始—对偶热启动（推荐，但要改 Moco 核心）

需要同时保存和恢复：

```text
x, lam_x, lam_g
```

其中 `lam_x` / `lam_g` 是 CasADi/IPOPT 的 NLP 对偶变量，**不是** `MocoTrajectory.multipliers`
（后者是运动学约束乘子，不是完整 NLP 对偶）。

修改点：

1. `CasOCTranscription::solve()` 从 `nlpResult` 取出 `lam_x`、`lam_g`，存到 `Solution`；
2. 给 `Solution` / `MocoSolution` 增加对偶字段，或写一个 sidecar 文件（例如 `.npz`）与 `.sto` 配套；
3. 下一次求解时，在 `nlpFunc` 的 `DMDict` 里加入：
   ```text
   lam_x0 = previous_lam_x
   lam_g0 = previous_lam_g
   ```
4. 同时在 `solverOptions` 里设置 `warm_start_init_point = yes`，以及合适的
   `warm_start_bound_push` / `warm_start_mult_bound_push`；
5. 保证变量/约束顺序、缩放、mesh、bounds 完全一致；如果 s 改变，约束不变但梯度变了，
   旧对偶只能作为近似初值。

伪代码：

```cpp
casadi::DMDict inputs{
    {"x0", x},
    {"lam_x0", lam_x},
    {"lam_g0", lam_g},
    {"lbx", lbx}, {"ubx", ubx},
    {"lbg", lbg}, {"ubg", ubg}
};
nlpFunc(inputs);
```

验收标准：同一 s、同一解重启时，primal 已经做到 228 次；完整 dual 应在此基础上进一步降低迭代数，
并把 iter 0 的 `inf_pr` 从 7.11 降向原解约束违背。如果 primal 已经足够快，dual 的收益可能有限。

#### 方案 C：不追 IPOPT 对偶，改用 continuation / 外层分块

如果不想改 Moco 核心，可以：

- 粗 mesh/低权重先解，再插值到细 mesh；
- 或继续现在的 fixed-s + 外层闭式 s 更新；
- 这些不是 IPOPT 意义的热启动，但能提供可检查的下降路径。

当前建议：primal guess 已经修好并验证（228 次）；先继续用外层分块，不必立即改核心。
如果后面确实需要反复求解同一结构、且 228 次仍太慢，再按方案 B 做 dual checkpoint。

---

## 4. 目前的主要结果

### 4.0 A 形式基线（B 形式的出发点）

- 1.5 s、w=20 的 A 形式 no-gain 主结果：mean r = **0.700**，RMSE = **0.097**，2/18 肌肉饱和，wall 约 1084.8 s。
- w = 1 / 5 / 10 / 20 的 1.5 s 扫描：mean r = 0.348 / 0.582 / 0.640 / 0.698；饱和数 2 / 4 / 3 / 2。
- A 形式 free-gain 在 w=20 时没有收敛（1000+ iter 后 objective 仍在 24 附近），这也是后来转向 B 形式的原因之一。
- A 的增益 g 也受 effort 偏置影响：g 越小，u 可以越小，所以大量通道的 g 被推到 0.5 下界。

### 4.1 0.3 s B 形式（当前最好已知点）

```text
joint-none, tol=1e-5
J = 2.25958
iterations = 2284
wall = 569.6 s
```

s：

```text
rect_fem_r   0.7602
vasti_r      2.0000
hamstrings_r 0.5365
tib_ant_r    1.4820
gastroc_r    0.5559
rect_fem_l   2.0000
vasti_l      2.0000
hamstrings_l 2.0000
tib_ant_l    1.2057
gastroc_l    0.5000
```

注意：2.25958 是“best known @ tol=1e-5”，不是严格证明的局部最优；
tol 从 1e-3 收到 1e-4、1e-5 时 J 还在下降，说明对容差敏感。

### 4.2 1.5 s fixed-s B 形式（有效结果，s 给定）

```text
s 固定为 0.3 s 的值
插件 fixed-s: J = 10.653509, 1555 iter, 904.6 s, 3/18 saturated
原生 fixed-s: J = 10.653837, 1505 iter, 893.8 s, 3/18 saturated
两者差 0.0031%
```

原生 fixed-s 用 reference \(e/s\) + per-control weight \(w s^2\)，没有任何 `MocoScaleFactor` 参数；
它和插件的 s* 几乎逐位一致。这一步同时验证了固定-s 的数学、XML 参数绑定和 cost 求值。

#### 4.2.1 1.5 s 外层闭式 s 更新（5 轮）

从 1.5 s fixed-s 解出发，固定 s 求 u，再用
\(s_i=\int e_i u_i/\int u_i^2\) 更新。前 3 轮用 `none`，第 4/5 轮改用原生无参数路径的 `random`（已验证，见 3.7）：

| 轮次 | objective | iterations | wall | 饱和 |
|---|---:|---:|---:|---:|
| 初始 s（0.3 s 值） | 10.653837 | 1505 | 893.8 s | 3/18 |
| 第 1 轮 | **10.465979** | 1626 | 1219.5 s | 4/18 |
| 第 2 轮 | **10.404849** | 1352 | 966.3 s | 5/18 |
| 第 3 轮 | **10.369963** | 1292 | 779.5 s | 5/18 |
| 第 4 轮（random） | **10.345139** | 1152 | 216.7 s | 5/18 |
| 第 5 轮（random） | **10.325442** | 1440 | 295.4 s | 5/18 |

- objective 单调下降：10.6538 → 10.4660 → 10.4048 → 10.3700 → 10.3451 → 10.3254；
- 第 3 轮相对第 2 轮下降 **0.335%**，第 2 轮相对第 1 轮下降 **0.584%**，都已经小于 1%；所以继续迭代的理由**不是 objective 没过 1% 阈值**，而是参数还没稳定；
- 原生无参数路径的 `random` 把第 4/5 轮从 13–20 min 降到 **3.6–4.9 min**，是实际可用的加速；
- 第 5 轮后 s* 仍在漂移（rect_fem_r 1.0613，vasti_r 1.9453，gastroc_l 0.5991），说明外层循环还没收敛；
- 但这证明了 **固定 s + 外层闭式更新** 是可行的、单调下降的路线，且比 free-s joint 的 40 min 不收敛更可控。

补充外层验收指标：

| 轮次 | J_total | r_J | R_s | G_s | effort | tracking | 饱和 |
|---|---:|---:|---:|---:|---:|---:|---:|
| s0 | 10.653837 | — | 18.92% | 0.1401 | 7.61172 | 3.04212 | 3/18 |
| s1 | 10.465979 | 1.763% | 8.77% | 0.04189 | 7.52526 | 2.94072 | 4/18 |
| s2 | 10.404849 | 0.584% | 5.81% | 0.02184 | 7.47728 | 2.92757 | 5/18 |
| s3 | 10.369963 | 0.335% | 4.99% | 0.01451 | 7.44613 | 2.92383 | 5/18 |
| s4 | 10.345139 | 0.574% | 4.37% | 0.01091 | 7.42343 | 2.92171 | 5/18 |
| s5 | 10.325442 | 0.190% | 4.00% | 0.00908 | 7.40492 | 2.92052 | 5/18 |

其中 \(R_s=\max_i|\hat s_i-s_i|/s_i\)，\(G_s=J(z,s)-J(z,\hat s)\) 是固定轨迹、只更新参数还能降多少。
**r_J 已经 <1%，但 R_s 仍约 4%、G_s 仍有 0.009，所以不能宣称参数已收敛；** 饱和从 3/18 增到 5/18，目标下降不能直接解释成生理合理性提高。

### 4.3 A/B 形状对比（0.3 s，以 cosθ 为主）

| muscle | cosB | cosA | 说明 |
|---|---:|---:|---|
| rect_fem_r | 0.948 | 0.952 | 基本持平 |
| vasti_r | **0.947** | 0.899 | B 更好 |
| hamstrings_r | **0.755** | 0.810 | B 更差 |
| tib_ant_r | **0.963** | 0.938 | B 更好 |
| gastroc_r | 0.992 | 0.995 | 持平 |
| rect_fem_l | 0.986 | 0.986 | 持平 |
| vasti_l | **0.833** | 0.804 | B 更好 |
| hamstrings_l | 0.961 | 0.961 | 持平 |
| tib_ant_l | 0.997 | 1.000 | 持平 |
| gastroc_l | **0.866** | 0.925 | B 更差 |

- cosθ 不受尺度影响，应作为主要指标；
- RMS 里 B 用 \(e/s\)、A 用 \(e\)，在 s=2 的通道上 B 的 RMS 会被尺度效应压低，不能直接当拟合改进；
- B 明显更好的是 vasti_r、tib_ant_r、vasti_l；明显更差的是 hamstrings_r、gastroc_l。

### 4.4 当前各个 s 一览

bounds 均为 \([0.5,2]\)。

| muscle | 0.3 s joint | 1.5 s 初始 s | 1.5 s s5 实际使用 | s5 后 s_next |
|---|---:|---:|---:|---:|
| rect_fem_r | 0.760193 | 0.760193 | **1.020472** | 1.061272 |
| vasti_r | 2.000000 | 2.000000 | 1.981360 | 1.945299 |
| hamstrings_r | 0.536452 | 0.536452 | 0.661691 | 0.673698 |
| tib_ant_r | 1.482010 | 1.482010 | 1.179160 | 1.158688 |
| gastroc_r | 0.555936 | 0.555936 | 0.701319 | 0.722808 |
| rect_fem_l | 2.000000 | 2.000000 | 1.554307 | 1.496362 |
| vasti_l | 2.000000 | 2.000000 | 1.398266 | 1.403025 |
| hamstrings_l | 2.000000 | 2.000000 | 1.691460 | 1.666919 |
| tib_ant_l | 1.205685 | 1.205685 | 1.341298 | 1.363331 |
| gastroc_l | 0.500000 | 0.500000 | 0.595813 | 0.599080 |

- **s5 实际使用**对应 J=10.325442，是当前 1.5 s 分块解；
- **s5 后 s_next** 只是固定 s5 轨迹算出的下一轮候选，**还没有实际求解**；
- 到 s5 时所有通道都已离开 0.5/2 的硬边界（vasti_r 1.981 接近上界），这是外层更新把贴界通道拉回内部的结果；
- 但 \(R_s\approx4.00\%\)、\(G_s\approx0.00908\)，说明 s 仍未稳定，不能把 s5 当成“辨识完成”的参数。

---

## 5. 现在卡在哪里

### 5.1 卡点

**1.5 s 的 free-s joint B 形式仍没有收敛。**

- fixed-s 1.5 s 能收敛（J=10.65）；
- 一旦把 10 个 s 放开，求解就在 16–20 之间振荡，40 min / 2586 iter 仍没有 EXIT；
- 这不是“维数高导致计算量大”（10 个参数只多 10 列），更可能是 **s 方向平坦、s 与 u 近似可补偿**，联合问题弱辨识/病态；
- 但“Hessian 近奇异”只是提示，**还没有证明**；`inf_pr` 是约束违背，`inf_du` 是对偶不可行性，不能单独证明 Hessian 病态；
- 另一个独立问题是 **初始 guess 的 inf_pr=4.37**：如果 A/B 只差目标函数，A 的可行轨迹理论上也应是 B 的可行轨迹。需要先查时间网格、插值、states/controls/derivatives 是否完整转移，而不是直接把 1.5 s 失败全部归因于平坦方向。

### 5.2 还不能下的结论

- 不能说“B 形式比 A 差”：0.3 s 只有部分通道更好，1.5 s 还没有可比的 joint 结果；
- 不能说“B 形式算不动”：fixed-s 1.5 s 15 min 收敛，慢的是 free-s joint；
- 不能把 0.3 s 的 s 当成 1.5 s 的 s：s* 漂移很大（vasti_l −0.38、rect_fem_l −0.15、tib_ant_r −0.10）；
- **s\* 是在 fixed-s 解的 u 不变时算出的条件最优值，不是完整 joint 解的 s**；它说明当前 s 不是这条轨迹下的驻点，但不能直接解释成“该方向可辨识”；
- A/B 的 objective 不是同一个目标的简单倒数参数化：\((e-su)^2=s^2(u-ge)^2\)，B 改变了误差的尺度，不能直接比 objective，只能比 cosθ、形状和同空间指标。

### 5.3 下一步建议（按专家意见调整顺序）

1. **参数绑定/导数验证**：已用 fixed-s 在 10 个非 1 的 s 值上做间接验证——插件与原生写法 1.5 s 差 0.0031%，s* 逐位一致；还可以补一个单通道 \(s=0.8/1.6\) 的显式 XML 测试，作为更直接的证据；
2. **原生 fixed-s 交叉验证**：已经完成。0.3 s 差 0.044%，1.5 s 差 0.0031%；原生写法是 reference \(e/s\) + per-control weight \(w s^2\)，不调用 `addScaleFactor()`；
3. **外层闭式 s 更新**：已经跑了 5 轮，J 从 10.6538 降到 10.3254；但第 5 轮 \(R_s\) 仍约 4.00%、\(G_s\) 仍有 0.00908，不能按“轮数”结束；建议继续 2–3 轮，直到 \(R_s<1\%\)、连续两轮 \(r_J<0.1\%\)、且 \(G_s\) 接近重复求解噪声。原生无参数路径的 `random` 已把每轮从 13–20 min 降到 3.6–4.9 min，可以继续用；
4. **不改变 s 的完整解重启对照**：已经做了，并找到根因：`--emg-initial-guess` 原来只在 rev 分支里 `setGuess`，原生路径根本没设。修复后 s5 自重启 iter 0 为 10.23/7.11，最终 **228 次迭代 / 53.9 s**，J=10.325376；primal warm start 已有效，dual 仍可选。后续外层可以放心用上一轮解作 guess；
5. 外层终点可以直接作为最终分块解，**不要求必须再做一次 joint 成功**；joint 只作为交叉检查，前提是轨迹块和参数块都满足一阶条件；
6. 若仍太平，再加弱先验 \(w_s(s-1)^2\)，并报告无先验/不同先验强度、贴界通道和 bounds 敏感性；
7. 比较 A/B 时以 cosθ 为主，RMS 要放到同一空间并标注贴界通道；
8. `read_moco_results` 绘图不必等 joint 成功：同时画 EMG 空间 \(e\) 与 \(su\)、control 空间 \(u\) 与 \(e/s\)，并标明 s 是 fixed、外层更新还是 joint。

---

## 6. 代码与结果文件

### 插件

```text
experiments/moco/revgoal/
  MocoReverseScaledControlTrackingGoal.h
  MocoReverseScaledControlTrackingGoal.cpp
  RegisterTypes_revgoal.cpp
  CMakeLists.txt
  build_vs/Release/osimRevGoal.dll
```

### Python

```text
experiments/moco/moco_inverse.py
```

关键参数：

```text
--emg-rev-goal
--emg-rev-plugin
--emg-initial-guess
--emg-s-values
--emg-sparsity-none
--emg-params-initsystem
```

### 主要输出

`experiments/moco/outputs/` 已清理为 38 个文件，保留以下代表结果、输入和说明：

```text
emg_excitation_199p5_204p0.sto        # sEMG reference input
grf_window_200_201p5.mot/.xml        # 1.5 s GRF input
ik_radians.mot                        # processed kinematics
MocoInverse结果文件说明.md            # 输出文件逐项说明
README_输出说明.md                    # 快速索引

moco_inverse_solution_emg200_w20_nogain.sto
moco_inverse_summary_emg200_w20_nogain.txt
moco_inverse_activation_report_emg200_w20_nogain.csv

moco_inverse_solution_rev_joint_tight1e5.sto
moco_inverse_summary_rev_joint_tight1e5.txt
moco_inverse_gains_rev_joint_tight1e5.csv

moco_inverse_solution_A_e2_compare.sto
moco_inverse_summary_A_e2_compare.txt
moco_inverse_solution_A_none_1p5.sto
moco_inverse_summary_A_none_1p5.txt

rev_fixs15_none.sto
moco_inverse_summary_rev_fixs15_none.txt
moco_inverse_gains_rev_fixs15_none.csv

native_b_15.sto
moco_inverse_summary_native_b_15.txt
native_b_15_s13.sto
moco_inverse_solution_native_b_15_s13.sto
moco_inverse_summary_native_b_15_s13.txt
moco_inverse_s_used_native_b_15.csv
moco_inverse_s_used_native_b_15_s13.csv
moco_inverse_activation_report_native_b_15_s13.csv
moco_inverse_residuals_native_b_15_s13.csv
emg_native_b_scaled_native_b_15_s13.sto

moco_inverse_solution_native_b_15_guess_fixed.sto
moco_inverse_summary_native_b_15_guess_fixed.txt
moco_inverse_summary_final_mesh001.txt
moco_inverse_summary_alt_s1_0.txt
moco_inverse_summary_alt_s1_1.txt
moco_inverse_summary_alt_s1_2.txt
final_1p5_tracking_report.csv
```

对应的 activation report 保留；中间 segments、失败/未收敛结果、旧版 random/initial-guess 结果、重复 scaled reference 已删除。

### 关键工程配置

- OpenSim/Moco：`4.6-2026-03-31-d6972f9e6`
- Python：`3.12.13`，env `opensim_scripting`
- 每次 B 形式求解：
  ```python
  solver.resetProblem(problem)
  solver.set_parameters_require_initsystem(False)
  ```
- 可信的 joint 求解：`optim_sparsity_detection="none"`
- `parameters_require_initsystem=True` 会 crash，不要用
- `initial-guess` sparsity 对这个 goal 不可靠

---

## 8. 正式收尾：最终 1.5 s 解与生物力学解读

### 8.1 正式技术路线

从现在起，MocoInverse 中 free-gain sEMG tracking 的正式路线是：

> **原生 fixed-s B 内层 + 外层闭式 s 更新 + `random` + 上一轮完整解作 guess。**

- 内层：`MocoControlTrackingGoal`，reference \(e/s\)，per-control weight \(w s^2\)，无 `MocoScaleFactor`；
- 外层：\(s_i=\mathrm{clip}(\int e_i u_i/\int u_i^2,0.5,2)\)；
- 每轮保存 `s_used` 和 `s_next`，不混用；
- joint、no-gain 不再作为主路线继续研究。

### 8.2 最终迭代

从 s5 继续，修复 guess 后每轮约 1–3 min。最终 13 轮：

| 轮次 | J_total | r_J | R_s | G_s |
|---|---:|---:|---:|---:|
| s5 | 10.325442 | — | 4.00% | 0.00908 |
| s6 | 10.308843 | 0.161% | 3.81% | 0.00806 |
| s7 | 10.294166 | 0.142% | 3.64% | 0.00712 |
| s8 | 10.281040 | 0.128% | 3.52% | 0.00634 |
| s9 | 10.269050 | 0.117% | 3.49% | 0.00607 |
| s10 | 10.257806 | 0.110% | 3.41% | 0.00559 |
| s11 | 10.247550 | **0.100%** | 3.31% | 0.00487 |
| s12 | 10.238180 | **0.091%** | 3.26% | — |
| s13 | **10.229322** | **0.087%** | **3.21%** | **0.00432** |

- 最后两轮 \(r_J<0.1\%\)，目标已稳定；
- 但 \(R_s\approx3.21\%\)，参数仍在漂移；\(G_s\approx0.00432\) 已很小（约总目标的 0.042%）；
- 因此 s13 是**轨迹/目标基本稳定的实用解**，不是“s 已唯一辨识”的证明。

最终 s13 实际使用：

```text
rect_fem_r   1.349567
vasti_r      1.623255
hamstrings_r 0.744262
tib_ant_r    1.107898
gastroc_r    0.862097
rect_fem_l   1.204092
vasti_l      1.475922
hamstrings_l 1.568192
tib_ant_l    1.463234
gastroc_l    0.604136
```

s13 后 `s_next` 仍会继续变化（例如 rect_fem_r 1.3929、vasti_r 1.5878、rect_fem_l 1.1753），
所以不能把 s13 当成跨窗口可迁移的生理增益。

### 8.3 1.5 s tracking 质量

同一 1.5 s / 20 ms / w=20 配置，比较：

- A：A 形式 free-gain + none；
- B0：B 形式 fixed-s 初始 s；
- B13：B 形式外层最终 s13。

均值：

| 指标 | A | B0 | B13 |
|---|---:|---:|---:|
| mean cosθ | 0.816 | 0.903 | **0.900** |
| mean Pearson r | 0.626 | 0.766 | **0.764** |
| mean control-space nRMSE | 0.734 | 0.400 | **0.380** |

B13 逐通道：

| muscle | cosθ | Pearson r | nRMSE |
|---|---:|---:|---:|
| rect_fem_r | 0.985 | 0.955 | 0.177 |
| vasti_r | 0.753 | 0.370 | 0.658 |
| hamstrings_r | 0.788 | 0.571 | 0.616 |
| tib_ant_r | 0.969 | 0.862 | 0.249 |
| gastroc_r | 0.931 | 0.819 | 0.365 |
| rect_fem_l | 0.965 | 0.882 | 0.264 |
| vasti_l | 0.746 | 0.467 | 0.666 |
| hamstrings_l | 0.899 | 0.833 | 0.438 |
| tib_ant_l | 0.989 | 0.932 | 0.147 |
| gastroc_l | 0.976 | 0.948 | 0.219 |

完整逐通道 CSV：`experiments/moco/outputs/final_1p5_tracking_report.csv`。

### 8.4 饱和与 reserves

最终 4/18 的 activation 峰值 >0.99：

```text
glut_max_l   1.000   1.3% 时间，最长约 0.02 s
iliopsoas_l  1.000   7.9% 时间，最长约 0.12 s
bifemsh_l    0.998   3.3% 时间，最长约 0.05 s
soleus_r     0.993   0.7% 时间，最长约 0.01 s
```

**这 4 个都不是被 EMG tracking 的 10 块肌肉**；被 tracking 的肌肉里，activation 峰值最高是
gastroc_l 0.978，没有 activation >0.99。gastroc_l 的 excitation 有约 2% 时间到 1.0，最长约 0.03 s。

Reserves 峰值（控制量）：

```text
hip flexion     R 0.85 / L 1.98
hip adduction   R 53.6 / L 66.2
hip rotation    R 25.5 / L 37.6
knee extension  R 14.8 / L 22.8
ankle dorsiflex R 3.9 / L 4.4
pelvis residual ty 167.9 N, tx 103.8 N
```

非矢状 reserves 的权重只有 1e-3，hip adduction/rotation 数值大是预期代价；但 knee 和 pelvis residual
仍然不低，报告结果时必须说明这是“EMG-informed 轨迹”，不是已通过全部肌骨验证的生理真值。

### 8.5 稳健性筛查

1. **网格 0.02 → 0.01**（固定 s13）：
   ```text
   J: 10.229322 → 10.21262（−0.16%）
   mean cos: 0.9001 → 0.8983
   mean nRMSE: 0.3798 → 0.3822
   逐通道最大差异 <0.02
   ```
   说明该 s 下的轨迹和 tracking 指标对网格不敏感；但 s 本身没有在细网格上重新辨识。

2. **不同初值 s=1**（只跑 2 轮）：
   ```text
   J: 10.834 → 10.5346 → 10.3888
   ```
   仍高于 s13 路径，且 s 分布不同（例如 vasti_r 1.34 vs 1.62）。这提示结果对初值/局部解有依赖，
   不能声称全局唯一。

### 8.6 生物力学解读

- **B 形式确实比 A 更贴近 EMG 形状**：mean cosθ 0.900 vs 0.816，Pearson 0.764 vs 0.626；
  改善最大的是 tib_ant、rect_fem、hamstrings、gastroc；vasti 仍然弱（r 0.37/0.47），
  说明这块肌肉的 EMG 形状与当前模型在 1.5 s 窗口内的动力学/努力分配并不完全一致。
- **s 的符号/大小有明确的“幅值校准”含义，但不是生理增益**：
  - s > 1 的通道（rect_fem_r 1.35、vasti 1.62、hamstrings_l 1.57、tib_ant_l 1.46）表示
    模型 excitation 需要比 EMG 低，才能同时满足动力学和 effort；
  - s < 1 的通道（hamstrings_r 0.74、gastroc 0.60–0.86）表示模型 excitation 需要比 EMG 高；
  - 这些是“在这个模型、这个窗口、这个权重下的有效缩放”，不能直接解释为神经增益或跨任务常数。
- **外层后期主要不是在改善 EMG 形状，而是在重新分配 effort**：
  - B0 → B13：总目标 −3.98%，tracking 项 −4.30%，effort 项 −3.86%；
  - 总下降中约 69% 来自 effort/非 tracking，31% 来自 tracking；
  - 逐通道 cos/r 有升有降（rect、hamstrings_r、tib_ant 升；vasti、rect_fem_l 略降），
    所以“J 下降”不能直接说成“每块肌肉都拟合得更好”。
- **饱和不是 EMG tracking 直接造成的**：4 个 activation 饱和的肌肉都不在 tracking 列表里；
  被 tracking 的肌肉没有 activation 饱和，只有 gastroc_l 的 excitation 短暂到 1.0。
- **reserves 仍是主要局限**：hip adduction/rotation、knee 和 pelvis residual 的峰值不小；
  这说明当前解依赖 residual/reserve 来满足动力学，不能据此声称肌肉力分布已完全可信。

### 8.7 结论

当前可以交付的是：

> **一个可复现的 1.5 s 原生 B fixed-s + 外层 s 更新流程，最终候选 J=10.229322，EMG 形状明显优于 A 基线；s 是弱可辨识的有效缩放，不是已辨识的生理增益；轨迹和 tracking 指标通过一次网格筛查，但 reserves 和 vasti 仍是限制。**

不适合现在声称：

- s 是唯一、可跨窗口迁移的生理参数；
- 肌肉力已经过独立验证；
- A/B 的 objective 可以直接比较。

---

### 8.8 绘图：read_moco_results 的 scale 处理

`read_moco_results.m` 现在按 sidecar 文件自动判断 scale 类型：

- 如果存在 `moco_inverse_s_used_<tag>.csv`（B 形式）：
  - `S.scaleKind = "B_s"`；
  - 图 1 的绿色曲线画 **`sEMG / s`**；
  - 每块肌肉标题显示 `s=...`；
  - 理由：\((e-su)^2=0\) 等价于 \(u=e/s\)，模型 excitation 的 scale 不动。
- 否则（A 形式）：
  - 仍然读 `moco_inverse_gains_<tag>.csv`；
  - 画 `sEMG .* g`，标题显示 `g=...`。

最终 B 解已有：

```text
moco_inverse_s_used_native_b_15_s13.csv
moco_inverse_solution_native_b_15_s13.sto
```

MATLAB 验证：`S.scaleKind = "B_s"`，且 `semgScaled - semg./s` 的最大差为 **0**。
这样不会把 `sEMG * s` 和 excitation 比较，也不会改变模型 excitation 的尺度。

---

## 9. 一句话总结

**正式路线已经确定为“原生 fixed-s B 内层 + 外层闭式 s 更新 + `random` + 上一轮解作 guess”。最终 1.5 s 候选 J=10.229322，mean cosθ=0.900、Pearson r=0.764，明显优于 A 基线（0.816/0.626）；s 是弱可辨识的有效缩放（R_s≈3.2%，s13 后仍在漂移），不是生理增益。轨迹和 tracking 指标通过一次网格筛查，但 vasti 拟合仍弱、hip/knee reserves 和 pelvis residual 仍大。该结果可作为当前模型和数据上的 EMG-informed 轨迹交付，不能声称 s 唯一可迁移或肌肉力已独立验证。**

---

## 10. 30 s 右腿结果（210–240 s）

### 10.1 配置

- 只跟踪右腿 5 块：`rect_fem_r, vasti_r, hamstrings_r, tib_ant_r, gastroc_r`；
- 左腿保留为 actuators，但无 EMG tracking；
- 每块右腿肌肉一个 s，整个 30 s 所有段共享；
- 4 chunk × 4 段并行，mesh 0.02；
- A 基线为 A no-gain（scale=1）。

### 10.2 最终候选

```text
moco_inverse_solution_B30_right_final.sto
moco_inverse_summary_B30_right_final.txt
moco_inverse_s_used_B30_right_final.csv
B30_right_tracking_report_final.csv
B30_right_outer_summary.csv
```

最终 s（第 7 轮实际使用值）：

```text
rect_fem_r   1.662868
vasti_r      1.493015
hamstrings_r 0.851994
tib_ant_r    1.177731
gastroc_r    1.068418
```

### 10.3 结果

| 指标 | A no-gain | B final |
|---|---:|---:|
| J_total | 208.303850 | **203.028772** |
| mean cos | 0.8698 | **0.8914** |
| mean Pearson r | 0.7034 | **0.7424** |
| mean nRMSE | 0.4219 | **0.3956** |

外层 **6 次有效更新**（r3 因脚本重启重复 s2，已排除）：

```text
J: 204.830963 -> 204.490759 -> 204.111790
   -> 203.792004 -> 203.515685 -> 203.244337 -> 203.028772
总下降 0.88%；每轮 r_J 0.11–0.19%
R_s: 3.97% -> 3.80% -> 3.67% -> 3.57% -> 3.47% -> 3.37% -> 3.27%
```

J/tracking 稳定，但 s 仍以约 0.10%/轮下降；**不能据此断言一定会撞 bounds**，
只能确认 6 次更新尚未收敛，现象与尺度补偿/慢收敛相符。

原始 u 检查：

```text
r6 -> r7: u RMS 变 0.3–3%，相关系数≈1，s·u 几乎不变
r0 -> r7: u RMS 变 −18% ~ +4%，s·u 只变 1–4%，形状相关 >0.995
```

所以 EMG 拟合稳定不代表原始 u 稳定；累计幅值确实变了，主要是 s 在补偿。
右腿未 tracking 肌肉 RMS 变化 ±3% 以内，reserves 峰值基本不变。

逐肌肉：B 对 vasti_r 改善最大（r 0.089→0.339，cos 0.584→0.733）；rect_fem_r、
tib_ant_r 也改善；hamstrings_r、gastroc_r 的 B 拟合反而比 A 差。

chunk 边界 seam max 0.0157（activation）。activation >0.99 只有 `soleus_r`
（3.3% 时间）和 `glut_max_l`（8.4%），都不是右腿 tracked 肌肉。

reserves 仍大：hip adduction L 62.5 / R 52.7，hip rotation L 33.9 / R 23.8，
knee L 29.1 / R 11.7；左腿虽不跟踪 EMG，仍参与 effort/动力学。

### 10.4 结论

> **右腿 5-s 的 30 s EMG-informed 轨迹；J/tracking 稳定，但 s 未收敛
> （R_s≈3.27%）；原始 u 累计幅值有变化，主要由 s 补偿；尚不能判定内部/边界收敛，
> 也不能声称参数唯一或肌肉力已验证。**

网格 0.01、初值 s=1、完整 30 s serial 交叉验证尚未做。
