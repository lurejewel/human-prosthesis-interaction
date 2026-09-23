# 用 Python 调用 OpenSim Moco / MocoInverse

本目录提供在本工作区用 Python 调用 OpenSim **MocoInverse** 的完整方案：环境自检脚本、
可直接运行的示例脚本，以及踩坑记录。

---

## 一、结论：可行，无需安装任何东西

本工作区**已经具备**运行 MocoInverse 的全部条件：

| 组件 | 状态 |
|------|------|
| conda 环境 | `opensim_scripting`（`D:\Software\miniconda_py312\envs\opensim_scripting`） |
| Python | 3.12.13 |
| OpenSim | 4.6（`4.6-2026-03-31-d6972f9e6`） |
| Moco | 4.6，78 个 `Moco*` 类，含 `MocoInverse` / `MocoCasADiSolver` |
| CasADi + IPOPT | 已随 conda 包一起安装（`Library\bin\casadi_nlpsol_ipopt.dll`） |
| numpy / scipy / matplotlib | 已安装 |

已验证：示例脚本在 SQR 步态数据上**成功求解**（`EXIT: Optimal Solution Found`，
约束违反 9e-4，51 次迭代，23 s）。

> ⚠️ 注意：系统默认的 `C:\Python314\python.exe`（Python 3.14）**不能**用 —— conda 版
> OpenSim 只提供 3.9–3.12 的绑定。请始终使用上面那个 `opensim_scripting` 的
> `python.exe`。

### 如果换一台机器，需要装什么

```bash
conda create -n opensim_scripting python=3.12
conda activate opensim_scripting
conda install -c opensim-org -c conda-forge opensim   # 自带 Moco + CasADi/IPOPT
conda install numpy scipy matplotlib                  # 脚本用到的辅助库
```

`opensim` 这个 conda 包同时包含 OpenSim 核心库、Moco 和 CasADi 求解器，**不需要**
另外 `pip install casadi`（Python 版 `casadi` 模块只有 `plot_casadi_sparsity.py`
才用得到）。

---

## 二、目录内容

| 文件 | 作用 |
|------|------|
| `moco_bootstrap.py` | 定位 CasADi 插件目录（解决 `ipopt` 找不到的问题），被下面两个脚本导入 |
| `moco_env_check.py` | **环境自检**：版本、求解器、依赖、输入文件、模型与 IK 一致性、运动学平滑度、GRF 覆盖、ModelProcessor 流程；`--solve` 会顺带解一个 1 自由度摆 |
| `moco_inverse.py` | **示例**：在 SQR 步态数据上真正跑一遍 MocoInverse |
| `outputs/` | 运行输出（解、激活、残差报告、日志） |

---

## 三、运行方法

```powershell
# 建议先激活环境（激活后 PATH 里就有 CasADi 插件了）
conda activate opensim_scripting

cd "C:\Users\lurej\OneDrive\Jin Wei\PKU\Muscle-Reflex-Based Predictive Simulation"

# 1) 环境自检
python experiments\moco\moco_env_check.py

# 2) 冒烟测试：0.3 s 窗口，几分钟内跑完
python experiments\moco\moco_inverse.py --smoke

# 3) 默认算例：180.0-181.5 s（约 1.4 个步态周期），mesh 2 ms（750 个网格间隔，
#    输出网格 1 ms）
python experiments\moco\moco_inverse.py

# 4) 指定窗口 / 网格（粗网格用于快速试算）
python experiments\moco\moco_inverse.py --start 180 --end 182 --mesh 0.02
```

脚本**不依赖 `conda activate`**：`moco_bootstrap.py` 会自己去解释器旁边找
CasADi 插件目录并设置 `CASADIPATH` / `PATH`。所以用
`D:\Software\miniconda_py312\envs\opensim_scripting\python.exe` 直接调用也可以。

### 常用参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `--start` / `--end` | 180.0 / 181.5 | 分析窗口 [s]，数据本身覆盖 180–240 s |
| `--mesh` | **0.002** | 配点网格间隔 [s]。默认 2 ms（约 5 个配点/激活时间常数，足够精确）；**0.001 在本 Moco 版本下无法收敛**（见下文「网格精度」）；快速试算用 `--mesh 0.02` |
| `--full` | — | 用满 60 s 窗口，**非常慢**，仅作说明 |
| `--smoke` | — | 0.3 s 窗口（默认网格）的极短算例 |
| `--filter-hz` | 15.0 | IK 坐标低通截止频率；**必须与 GRF 的截止频率一致**（见第五节第 6 条）；`0` 表示不滤波 |
| `--residual-rot` | 400 | 骨盆转动残差作动器容量 [N·m] |
| `--residual-trans` | 500 | 骨盆平动残差作动器容量 [N] |
| `--keep-contacts` | — | 保留模型自带的 HuntCrossley 接触力（一般不要） |
| `--no-dgf` | — | 保留原始 Millard2012 肌肉。**MocoInverse 用不了**，见第五节第 8 条 |
| `--tol` / `--max-iter` | 1e-3 / 3000 | 求解器收敛容差与最大迭代 |
| `--allow-large` | — | 覆盖内存护栏（不建议，见下文「内存与窗口长度」） |
| `--tag` | — | 输出文件名后缀 |

求解耗时会在求解结束时**直接打印在终端**（`>>> solve wall time: ...`，并区分
IPOPT 纯求解时间 `solver duration`），同时写入 `moco_inverse_summary*.txt` 的
`wall_time_s=` / `solver_duration_s=` 字段。

#### 内存与窗口长度（重要）

CasADi 的 NLP 规模 = 网格间隔数 × 每个时间点的变量数（本模型约 93 个状态+控制）。
本机实测 **每个网格间隔约需 2.9 MB 内存**：

| 窗口 × 网格 | 间隔数 | 内存估算 | 结论 |
|-------------|-------|---------|------|
| 1.5 s × 2 ms | 750 | ~2.2 GB | ✅ 已收敛（约 2 min） |
| **2 s 片段 × 2 ms（当前默认，实际求解 3 s）** | **1,500** | **~4.2 GB** | ✅ 已收敛（约 3–5 min/段） |
| 20 s × 2 ms（整段求解） | 10,000 | ~29 GB | ❌ `bad allocation` 崩溃 |
| 20 s × 20 ms | 1,000 | ~2.9 GB | ⚠️ 能收敛，但激励/激活有锯齿，**不可用于肌肉结论** |

脚本在启动前**预估内存并拒绝单段 > 3000 间隔的问题**（可用 `--allow-large`
强行尝试）。长窗口的解决办法就是下面的**分段求解**，不要靠放粗网格。

#### 长窗口：分段求解（默认行为）

`--segment 2 --overlap 0.5`（默认）把 `[--start, --end]` 切成若干 2 s 片段，
每段单独求解后缝合。**关键点是段与段之间没有"时间间隙"，也不需要人为造连续性**：

- **没有间隙**：第 k 段实际求解 `[t_k − 0.5, t_{k+1} + 0.5]`，只保留中间的
  `[t_k, t_{k+1}]`。保留的片段首尾严格相接、**精确平铺**整个窗口（
  `test_segmentation.py` 里有断言验证：无缝、无重叠）。
- **为什么要多解前后各 0.5 s**：MocoInverse 求解的变量分两类——
  - **有记忆的状态**：肌肉激活 a，服从 `da/dt = (x − a)/τ`
    （DGF：τ_act = 10 ms，τ_deact = 40 ms）；
  - **无记忆的量**：兴奋 x、残差/储备控制、以及被样条约束的关节运动学。
  分段时每段的**初始激活是自由变量**（MocoInverse 只约束 a(0)=x(0)，不指定数值），
  所以相邻两段会各自挑一个不同的 a(t_k) → 激活曲线在接缝跳变，而激活一跳、
  肌肉力就跳。多解出来的那 0.5 s 就是**预烧（burn-in）**：它含有"错误初始条件"
  的瞬态，直接丢掉。0.5 s = 12.5 τ_deact，瞬态被指数衰减压到可忽略。
- **瞬态到底衰减多快（实测）**：在 182 s 的接缝上对比两段——

  | 距该段自身求解边界 | max｜Δa｜（激活单位，范围 0–1） |
  |---|---|
  | 0.00 s（段起点，a(0)=x(0) 瞬态） | 0.958 |
  | 0.10 s | 0.005 |
  | 0.20 s | 1e-4 |
  | 0.00 s（段**终点**，自由终端横截条件边界层） | 0.307 |
  | ≥0.10 s（两端都远离边界） | 1e-4 – 4e-4 |

  所以 0.5 s 的 burn-in 有 5 倍余量；脚本在 `--overlap < 0.2` 时会警告。
- **不连续有多大：脚本自己量**。相邻两段在接缝附近有一段重叠，脚本在**对称于
  切口**、且两段都距各自求解边界 ≥ `SEAM_SAFE_MARGIN = 0.4 s` 的窗口内比较它们，
  输出 `moco_inverse_seams.csv` 并在终端打印：激活差的 RMS/最大值、兴奋差的
  RMS/最大值，以及**缝合后序列在切口处的实际跳变量**。判据：`act max` 应与激活
  本身的量级（0–1）相比可忽略（本数据实测 ~4e-4）。若它变大，加大 `--overlap`。
- **没采用的替代方案**：把上一段末尾的激活固定为下一段的初值（链式初值）。它需要
  走 `MocoInverse.initialize()` 自定义流程（该路径构造的问题物理与 `solve()`
  不一致，风险高），而且会把每段误差**累积**下去；burn-in 方案每段都重新"忘掉"
  初值，误差不累积。
- 每段的区间、保留区间和耗时都写进 `moco_inverse_summary.txt`，每段的原始解也单独
  存到 `outputs/segments/seg_XX.sto`，便于逐段核查。

### 输出文件

| 文件 | 内容 |
|------|------|
| `moco_inverse_solution*.sto` | 缝合后的完整轨迹：**状态 + 控制**（56 状态 + 37 控制） |
| `moco_inverse_activations*.sto` | 只含肌肉**激活**（state），方便和 EMG 对比 |
| `moco_inverse_excitations*.sto` | 只含肌肉**兴奋**（control） |
| `moco_inverse_residuals*.csv` | 残差/储备作动器的峰值控制与出力 |
| `moco_inverse_activation_report*.csv` | 每块肌肉的峰值/平均激活 |
| `moco_inverse_seams*.csv` | **接缝连续性诊断**（见上） |
| `moco_inverse_summary*.txt` | 求解摘要 + 每段区间与耗时 |
| `segments/seg_XX.sto` | 每段的**原始解**。中间产物（10 段约 51 MB），**可随时删除**，下次运行会重新生成；`--no-segment-files` 可完全不写 |
| `grf_window_*.mot/.xml` | 裁到窗口的 GRF |
| `ik_radians.mot` | 滤波+弧度化的运动学（**`read_moco_results.m` 读它来画图 2/3 的 IK 对比，别删**） |

#### 运行时产生的临时文件（可随时删）

- `segments/`：每段的原始解，中间产物，下次运行会重新生成；`--no-segment-files` 可完全不写。
- `grf_window_*.mot/.xml`：每次运行的窗口化 GRF，下一次运行会覆盖/新建。
- `delete_this_to_stop_optimization__<时间戳>.txt`：**MocoCasADiSolver 每次求解都会在工作目录
  创建**的哨兵文件（在仓库根目录，因为那里是运行时的工作目录）。它不参与计算，可以随时删除；
  反过来，**在求解过程中删掉它就是让 Moco 优雅停止的官方方法**。已加入 `.gitignore`。
- `<仓库根>/opensim.log`：OpenSim 运行日志（`*.log` 已被 gitignore）。注意它是**追加**的，
  你 MATLAB 的 forward-simulation 也往同一个文件写，删之前留意。
- 注意 `moco_inverse_activations.sto` / `moco_inverse_excitations.sto` 是**只含部分列的
  表格**，头部 `num_states`/`num_controls` 必须与实际列数一致（Moco 要求三者之和 = 列数），
  否则 `osim.MocoTrajectory` 会以 `bad allocation` 报错——脚本里已加断言防复发。

#### 兴奋（excitation）和激活（activation）的关系

MocoInverse **同时解这两者**，它们不是一回事：

- **兴奋 x 是控制变量**（完整解文件里裸的肌肉名列，如 `/forceset/soleus_r`）；
- **激活 a 是状态变量**（`/forceset/soleus_r/activation`），由一阶激活动力学
  `da/dt = (x − a)/τ_a` 决定，且约束 `a(0) = x(0)`。

所以 `moco_inverse_solution.sto` 里两列都有（兴奋在后半段的控制区）；之前生成的
`moco_inverse_activations.sto` 是**故意只导出激活**的辅助文件（为方便与 EMG 对比），
另外还导出了 `moco_inverse_excitations.sto`（只含兴奋）。`_tag` 后缀（如 `_final`）
来自脚本的 `--tag` 参数，只是同一次/另一次运行的标记，内容结构完全一样——建议只保留
一组无后缀的结果，避免混淆。

### 在 MATLAB 里读取结果

仓库自带的 `functions/utils/read_sto_file.m` 依赖 OpenSim 的 MATLAB 接口，而 Moco 的
输出不需要它。本目录的 `read_moco_results.m` 是**纯 MATLAB 读取器**（无需配置任何
OpenSim-MATLAB 接口），读取上面所有文件并绘图：

```matlab
cd('experiments/moco')
S = read_moco_results();                                    % 读取 + 绘图
S = read_moco_results(struct('Plot', false));               % 只读数据
S = read_moco_results(struct('Tag', 'millard'));            % 读带后缀的一组
S = read_moco_results(struct('SaveFigures', true));         % 把图存成 PNG
```

返回的 `S` 里可以直接用：`S.time`、`S.muscles`、`S.activation`（nT×18）、
`S.excitation`（nT×18，与激活一一对应）、`S.coordNames` / `S.coordValue` /
`S.coordSpeed`、`S.actuatorNames` / `S.actuatorControl`、`S.residuals`、
`S.activationReport`、`S.summary`、`S.ikValue` / `S.ikRawValue`（IK 参考，插值到
`S.time`）。四张图分别是：① 激活 vs 兴奋（18 块肌肉）② 坐标值 ③ 坐标速度
④ 残差/储备作动器报告。

**图 2/图 3 会把 Moco 与 IK 画在一起**：实线 = Moco 解，虚线 = 喂给 Moco 的滤波后
IK（`outputs/ik_radians.mot`），点线 = 原始未滤波 IK；每个子图标题给出 Moco 与滤波
IK 的 RMS 差异，控制台也会打印逐坐标 RMS 表。注意：MocoInverse 把坐标**约束成样条**
（滤波 IK 的插值），所以 Moco 与滤波 IK 几乎完全重合（RMS ≈ 1e-14），肉眼可见的差异
来自与**原始未滤波 IK** 的对比，即 15 Hz 低通滤波去掉的高频成分。

---

## 四、这个流程做了什么

MocoInverse 解决的是**肌肉冗余问题**：给定实测运动学（以及实测地面反力），求能产生
该运动的肌肉激活/兴奋轨迹。和 `experiments/scripts/run_ik_id.py` 里的逆动力学互补 ——
逆动力学得到关节力矩，MocoInverse 得到肌肉层面结果。

脚本的处理链：

1. **模型处理**（`ModelProcessor`）
   - 载入 `SQR_simbody.osim`
   - `ModOpAddExternalLoads` 挂上实测 GRF
   - `ModOpIgnoreTendonCompliance` + `ModOpReplaceMusclesWithDeGrooteFregly2016`
     + `ModOpIgnorePassiveFiberForcesDGF` + `ModOpScaleActiveFiberForceCurveWidthDGF(1.5)`
     —— DeGrooteFregly2016 是 Moco 推荐、数值条件更好的肌肉模型
   - `ModOpAddResiduals` 给 6 个骨盆自由度加残差作动器
   - `ModOpAddReserves(1.0)` 给其余 13 个自由度（含无肌肉跨越的 3 个腰椎自由度）加储备作动器
2. **关闭模型自带的 HuntCrossley 足部接触力** —— 已经用了实测 GRF，接触模型再算一遍
   就是重复计力
3. **预处理运动学**：15 Hz 零相位低通滤波 + **角度由度转弧度**
4. **预处理 GRF**：裁到分析窗口
5. `osim.MocoInverse` 求解，写出解、激活、残差报告

---

## 五、踩坑记录（都已在脚本里处理）

这几条是本工作区实测出来的，按"踩到的顺序"排列。如果你自己写脚本，很可能撞上同样的坑。

### 1. `MocoCasADiSolver.isAvailable()` 返回 True，但一求解就报 `Plugin 'ipopt' is not found`

`isAvailable()` 是**编译期**标志，不代表运行期能加载到插件。conda 把 CasADi 的
插件 DLL 放在 `<env>\Library\bin`，只有 `conda activate` 之后才会进 `PATH`。直接调用
`python.exe`（VS Code 运行按钮、任务计划、`subprocess`……）就会漏掉。

→ `moco_bootstrap.py` 自动定位并设置 `CASADIPATH` / `PATH` / `os.add_dll_directory`。
这是 conda 打包的已知问题，见
[opensim-org/conda-opensim#44](https://github.com/opensim-org/conda-opensim/issues/44)。

### 2. IK 文件是「度」，Moco 要的是「弧度」—— 差 57.3 倍

`level_walking_ik_filtered.mot` 头部写着 `inDegrees=yes`，`knee_extension_r` 的值是 **-67.0（度）**。
但是 `osim.TimeSeriesTable` **不会**根据 `inDegrees` 做转换，它原样返回数值；而 Moco 把
坐标值当 **SimTK 内部单位（弧度）** 用。直接把文件喂进去，等于要求膝盖转到 -67 rad
（≈ -3840°）。

已实测确认：`TimeSeriesTable` 的读数与文本原值逐位相同，没有转换。

→ `prepare_kinematics_radians()` 用**模型**的 `Coordinate.getMotionType()` 判断哪些是
转动坐标（而不是靠名字猜），只对这些列做 `deg2rad`，并把头部改成 `inDegrees=no`。
`pelvis_tx/ty/tz` 是平动，保持不变。

### 3. 骨盆残差作动器容量太小 → 问题**不可行**

这是最难查的一个：现象是 IPOPT 迭代 3000 次也不收敛，约束违反量在 1e2–1e8 之间乱跳，
解里 `residual_..._pelvis_tx/ty/tz` 的控制量**恰好等于 1.0000**（顶到边界）。

根因：`ModOpAddResiduals(rot, trans, bound)` 里，`bound=1.0` 把控制限在 [-1, 1]，
所以平动残差的最大出力就是 `trans × 1.0`。而用本工作区已有的逆动力学结果
（`ResultsInverseDynamics/level_walking_id.sto`）一查就知道，这位受试者在 180–181.5 s
内需要：

| 骨盆自由度 | 实测需求 | 原设定 `trans=50` | 现设定 `trans=500` |
|-----------|---------|------------------|-------------------|
| `pelvis_tx_force` | **131.9 N** | 50 N ❌ 不可行 | 500 N ✓ |
| `pelvis_ty_force` | **98.1 N** | 50 N ❌ 不可行 | 500 N ✓ |
| `pelvis_tz_force` | 52.5 N | 50 N ❌ 不可行 | 500 N ✓ |
| `pelvis_tilt_moment` | 78.0 N·m | 250 N·m ✓ | 400 N·m ✓ |

OpenSim 官方示例用 `ModOpAddResiduals(250.0, 50.0, 1.0)`，但对这份数据 50 N 明显不够。
改成 `(400, 500, 1.0)` 之后立刻收敛。

改容差**没用**——问题不是收敛慢，而是本来就没有可行解。**判断技巧：看解里残差/储备
作动器的控制量是否顶在边界上；顶住了就说明容量不够。**

收敛后的残差占用率（默认算例，1.5 s 窗口）：

| 作动器 | 峰值控制 | 占容量 |
|--------|---------|--------|
| `residual_..._pelvis_tilt` | 0.09 | 9% |
| `residual_..._pelvis_tx` | 0.25 | 25% |
| `residual_..._pelvis_ty` | 0.19 | 19% |
| `residual_..._pelvis_tz` | 0.12 | 12% |

### 4. GRF 文件太大 → 求解前先花几分钟读文件

`level_walking_grf_filtered.mot` 有 249,600 行（0–249.6 s，约 54 MB 文本）。Moco 在构建
优化问题的过程中会**反复初始化模型**，每次初始化 `ExternalForce` 都要把这个 Storage
完整解析一遍。实测一次求解前发生了 **114 次全量读取 ≈ 6 GB 文本解析**。

→ `prepare_grf_window()` 把 GRF 裁到分析窗口（两侧各留 1 s），行数从 249,600 降到
2,301，**解完全不变**。

### 5. 模型自带的接触力会和实测 GRF 重复计力

`SQR_simbody.osim` 里有 `HuntCrossleyForce`（`foot_r` / `foot_l`，stiffness 5e6）。
既然已经把实测 GRF 作为外力加上去了，接触模型就不能再算一遍。
→ 脚本默认用 `set_appliesForce(False)` 关掉它们（保留几何，便于随时恢复）。

### 6. IK 未滤波 → 二阶导爆炸

MocoInverse 对坐标拟合样条并**求两次导**，所以 IK 里的噪声会被放大成巨大的关节角加速度。
本数据 180–181.5 s 窗口内的峰值：

| 坐标 | 原始 | 6 Hz 滤波后 | 15 Hz 滤波后（当前） |
|------|------|------------|--------------------|
| `ankle_dorsiflexion_l` | 24,488 °/s² | 6,468 °/s² | 17,116 °/s² |
| `ankle_dorsiflexion_r` | 19,881 °/s² | 7,050 °/s² | 16,473 °/s² |
| `knee_extension_l` | 10,606 °/s² | 6,665 °/s² | 10,270 °/s² |

→ 默认按 **15 Hz** 零相位 Butterworth 滤波，和 `Setup_InverseDynamics.xml` 里逆动力学用的
截止频率、以及 `grf_filter.py` / `ik_filter.py` 里 GRF 与运动学的截止频率**三者必须一致**：
逆动力学与 MocoInverse
都要把「二次求导后的坐标」和「外力」相加，两者带限不同则残差不能相消。
代价是角加速度比 6 Hz 时大约 2.5 倍（上表），这是为了让动力学自洽而接受的。
`--filter-hz 0` 可关闭。

### 7. 求解失败时解是「封存」的，写盘会抛异常

求解不收敛时 `MocoTrajectory.write()` 会报 *"This trajectory is sealed"*。脚本会在
失败时自动 `unseal()`，这样最后一步迭代结果仍可写出来检查。

### 8. MocoInverse 实际上必须用 DeGrooteFregly2016 肌肉

试过 `--no-dgf`（保留模型原始的 Millard2012EquilibriumMuscle），会在构建问题时直接失败：

```
RuntimeError: ... No info available for state '/forceset/hamstrings_r/fiber_length'.
```

Millard2012 不提供 MocoInverse 初始化所需的肌纤维长度状态信息。所以
`ModOpReplaceMusclesWithDeGrooteFregly2016()` 不是可选项，而是前提。
`--no-dgf` 仅保留用于实验和演示这个限制。

---

## 六、结果解读 / 怎么判断结果可信

### 残差/储备作动器（`outputs/moco_inverse_residuals*.csv`）

MocoInverse 是"给定运动学反推肌肉"，模型和实测数据之间必然有出入，这些出入由残差/
储备作动器吸收。判断标准：

- **残差控制量远小于 1**（比如 < 0.3）→ 模型能很好地复现该运动
- **残差控制量接近或等于 1** → 容量不够或模型不匹配；先加大 `--residual-rot/-trans`
  看是否缓解；如果加大后残差**仍然很大**，说明模型/GRF/运动学三者不自洽
- **储备作动器控制量很大**（> 10）→ 对应关节的肌肉无法产生所需力矩

默认算例（180–181.5 s）的结果，和 `ResultsInverseDynamics/level_walking_id.sto`
里独立算出来的逆动力学结果对比 —— 两者吻合得很好，说明整个设置是对的：

| 骨盆自由度 | MocoInverse 残差出力 | 逆动力学需求 |
|-----------|---------------------|-------------|
| `pelvis_tx` | 123.0 N | 131.9 N |
| `pelvis_ty` | 94.8 N | 98.1 N |
| `pelvis_tz` | 61.3 N | 52.5 N |
| `pelvis_list` | 44.3 N·m | 43.5 N·m |
| `pelvis_rotation` | 21.3 N·m | 19.4 N·m |

残差占容量比例：tx 25%、ty 19%、tz 12%、tilt 11%、list 11%、rotation 5% —— 均远低于
饱和，属于健康范围。

### 肌肉激活（`outputs/moco_inverse_activation_report*.csv`）

**这是本数据最需要注意的一点。** 默认算例里 **18 块肌肉中有 12 块的峰值激活达到 1.0
（饱和）**，因此储备作动器承担了相当大的关节力矩（髋内收 67 N·m、髋屈曲 49 N·m、
踝背屈 40 N·m、膝 44 N·m）。

也就是说：**`SQR_simbody.osim` 这个 18 肌肉简化模型，对于该受试者这段步态的力矩需求
来说偏弱。** 这不是脚本的问题（残差与逆动力学吻合、骨盆残差远未饱和都证明了设置正确），
而是模型本身的能力上限。脚本会在饱和肌肉超过 30% 时自动打印提示。

如果要做肌肉层面的定量结论，建议先解决这一点，例如：给模型做肌肉强度标定
（`ModOpScaleMaxIsometricForce`）、换用更完整的模型（如
`model/NOT YET TO USE/Rajagopal2016_OpenSim4.5.osim`，但需要配套 IK），或者明确
在论文里说明储配作动器的贡献比例。

### 网格精度（mesh 选择）

实测结果（默认 1.5 s 窗口，本机）：

| mesh | 网格间隔数 | 结果 |
|------|-----------|------|
| 0.02 | 75 | 收敛，38.8 s（68 次迭代） |
| **0.002（默认）** | 750 | **收敛，187.8 s**（93 次迭代，IPOPT 180.9 s），输出网格 1 ms |
| 0.001 | 1500 | **不收敛**：IPOPT 约束违反卡在 ~60 后进入 restoration 模式并发散（约 30 min 仍无解，已终止） |

> ⚠️ **重要更正**：本节早先版本根据"骨盆残差在 0.02 与 0.002 下一致"就断言 2 ms 足够、
> 20 ms 也可用。**这个结论是错的**——残差只反映整体动力学，掩盖了肌肉层面的网格尺度
> 振荡。请以下面的实测为准。

**网格太粗会让激励/激活出现锯齿（chatter）。** 用网格尺度上的局部交替分量
`alt_i = y_i − (y_{i−1}+y_{i+1})/2`（在**配点网格**上计算，跨网格可比）实测：

| 网格 | alt/std 激励 x | alt/std 激活 a | a 的交替 RMS（绝对量） |
|------|------|------|------|
| 0.02 s | 0.412 | **0.237**（18 块肌肉中 13 块判为 CHATTERING） | 最大 **0.110** |
| 0.002 s | 0.121 | **0.011** | 最大 0.004 |

原因：DeGrooteFregly2016 的**激活时间常数 τ_act = 10 ms**（去激活 40 ms），而 h = 20 ms
时 **h/τ = 2**。梯形配点在该比例下存在**近交替的零空间方向**（相邻两点的激励一升一降，
平均后激活几乎不变），而代价函数只惩罚激励的**幅值**（`excitation_effort = Σx²`）、
不惩罚其变化率，于是优化器会激发该模式。网格细化到 2 ms（h/τ = 0.2）后**激活的锯齿
降约 20 倍**，肌肉力才可信。

已排除的其他原因：
- **不是外力造成的**：`level_walking_grf_filtered` 在窗口内是带限的
  （15 Hz 零相位低通），换滤波参数锯齿依旧；
- **`minimize_sum_squared_activations=True` 无效**：实测 alt RMS a 0.0636 → 0.0688；
- **不是"网格太细"**：细化网格是**改善**而非恶化。

**结论**：肌肉层面的结论（激活、兴奋、肌肉力、与 EMG 对比）**必须用 ≤ 2 ms 网格**。
`--mesh 0.02` 只适合冒烟测试和快速探索，不能用于定量结果。20 s 窗口在 2 ms 下需要
10,000 个配点间隔（约 29 GB 内存，本机不可行），所以长窗口要么**分段求解**
（例如每段 2 s ≈ 1000 个间隔 ≈ 2.9 GB、各约 3 min），要么在论文里明确说明只用了粗网格。

**1 ms 不收敛**是另一回事（数值条件限制）：Moco 用前向有限差分求导，网格到 1 ms 时
KKT 系统条件数恶化，IPOPT 无法推进。

**结果质量自检**：`check_solution_quality.py` 会报告每个解的锯齿指标、饱和肌肉数和
残差是否顶到边界，并直接给出可用/不可用的结论：

```powershell
python experiments\moco\check_solution_quality.py                       # 默认解
python experiments\moco\check_solution_quality.py outputs\moco_inverse_solution_xxx.sto
```

---

## 七、和其他脚本的关系

```
experiments/scripts/run_ik_id.py     IK + ID（关节力矩）
        │  level_walking_ik.mot         (原始 IK，未滤波)
        │  level_walking_ik_filtered.mot (15 Hz，由 ik_filter.py 生成)
        │  level_walking_grf_filtered.xml
        ▼
experiments/moco/moco_inverse.py    MocoInverse（肌肉激活/兴奋）
```

如果重新跑了 IK（例如改了 `Setup_IK.xml`），MocoInverse 会自动用新的
`level_walking_ik_filtered.mot`，因为脚本每次都会重新做滤波+单位转换到 `outputs/`。

注意两个 IK 文件的区别：`level_walking_ik.mot` 是 `InverseKinematicsTool` 的**原始**
输出，`level_walking_ik_filtered.mot` 是它经过 `experiments/scripts/ik_filter.py`
（15 Hz、4 阶 Butterworth、零相位）滤波后的结果。逆动力学、RRA、MocoInverse 一律读
**filtered** 那个。

## 八、其他可用的模型

`model/human0918.osim`（18 肌肉，前向仿真在用的模型）同样可以用来做 MocoInverse，
但**需要先有与之匹配的 IK 结果**（即该模型坐标系下的坐标 .mot 文件）。目前仓库里
没有这份数据，所以示例用的是 `SQR_simbody.osim` + `level_walking_ik_filtered.mot`
这一套（二者坐标完全对应，19/19 匹配）。要换模型，改脚本顶部的 `MODEL_FILE` /
`IK_FILE` 即可，其余逻辑不变。
