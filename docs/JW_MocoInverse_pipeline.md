# JW_* 七组数据 → MocoInverse 肌肉力：前置流程梳理

对象：`experiments/data/JW_{walking1p0, fastwalking1p4, run, uphillwalking0p8, uphillwalking1p0, downhillwalking0p8, downhillwalking1p0}`。
参照：`experiments/data/SQR_walking` 已跑通的链路（`experiments/scripts/` + `experiments/moco/`）。
环境：`D:\Software\miniconda_py312\envs\opensim_scripting\python.exe`（Python 3.12，含 OpenSim/Moco/CasADi-IPOPT、numpy/scipy/matplotlib）。系统 `C:\Python314` 不可用。

---

## 0. 数据核对（已实测，非假设）

| 组 | 时长 (s) | GRF fs | IK fs | sEMG fs | 稳定段 (s) | 分析窗 = 稳定段末 30 s | 步频 (步/min) |
|---|---|---|---|---|---|---|---|
| JW_walking1p0 | 185.29 | 2000 | 100 | 2000 | 5.69–184.69 | **154.69–184.69** | 102.4 |
| JW_fastwalking1p4 | 92.94 | 2000 | 100 | 2000 | 0.44–92.44 | **62.44–92.44** | 121.2 |
| JW_run | 102.99 | 2000 | 100 | 2000 | 4.41–95.41 | **65.41–95.41** | 150.0 |
| JW_uphillwalking0p8 | 70.11 | 2000 | 100 | 2000 | 0.31–69.31 | **39.31–69.31** | 97.6 |
| JW_uphillwalking1p0 | 69.93 | 2000 | 100 | 2000 | 0.42–69.42 | **39.42–69.42** | 111.8 |
| JW_downhillwalking0p8 | 74.63 | 2000 | 100 | 2000 | 0.27–74.27 | **44.27–74.27** | 99.5 |
| JW_downhillwalking1p0 | 78.52 | 2000 | 100 | 2000 | 0.34–78.34 | **48.34–78.34** | 111.4 |

- GRF `.mot` 19 列，与 SQR 同布局：`time, ground_force{1,2}_{vx,vy,vz,px,py,pz}, ground_torque{1,2}_{x,y,z}`；**belt1=左、belt2=右**（与各组 `*_grf.xml` 的 `force_identifier` 一致）。单位 N / m / N·m，ground 系。
- IK `.mot` 20 列、`inDegrees=yes`，19 个坐标与 `JW_simbody.osim` 一致。
- 步频×步长反推与跑台速度自洽（1.0 m/s→0.59 m 步长；2.0 m/s 跑→0.80 m），说明着地检测正确。
- 每组 `notes.txt` 已写入稳定段与分析窗（见 §5）。

### 0.1 两个必须先知道的数据问题

**(A) 四个坡度组的 GRF 相对铅垂倾斜 15.05° —— 这是预期的，不需要校正。**（已确认）
世界坐标系原点设在跑台上，跑台真实倾斜而 mocap 不动，因此**标记点与力数据都表达在跑台对齐系**里：数据上看跑台是平的、mocap/力反向倾斜。各组 `JW_simbody.osim` 的**重力已预旋转**匹配（平地 `[0,-9.80665,0]`；上坡 `[-2.5381,-9.4725,0]`；下坡 `[+2.5381,-9.4725,0]`，模长均 9.8066），动力学自洽。
结论：**按部就班执行、不做任何反旋转**；唯一要求是 MocoInverse 取用**本组自己的** `JW_simbody.osim`（`trial_paths.py` 已保证）。`derotate_grf.py` 仅作诊断报告（`grf_frame_report.json`），**禁止**用它改力数据。

**(B) SQR 的 heel-strike 剪切力符号修正不能照搬到 JW。**
`fix_grf_heel_strike.py` 的规则（左 vx≤0/vz≥0，右 vx≤0/vz≤0）编码的是 SQR 朝 **+x** 走的制动剪切方向。JW 跑带朝 **−x**（实测带速 −0.8…−2.0 m/s），早期 stance 的 vx 本应为正。实测照搬会"修正" **611 段**（≈全部着地），即把正确信号反转。
结论：**JW 默认跳过 `fix` 阶段**；15 Hz stance 分段滤波 + 去抖 + taper 已足够处理冲击尖峰。`run_pipeline.py` 默认 stages 不含 `fix`。

---

## 1. 前置工作总览（按执行顺序）

```
raw *_grf.mot ──[1 GRF 滤波]──> *_grf_filtered.mot/.xml ─
raw *_ik.mot  ──[2 IK 滤波]───> *_ik_filtered.mot        ─┼─> [6 MocoInverse]
raw *_sEMG.csv ─[3 sEMG 链]───> *_semg_processed.csv ─[4 重采样]─> *_emg_reference.sto ─
stable window ─[5 稳定区间]───> notes.txt / steady_state_summary.json ─┘
```

1. **GRF 滤波**（`grf_filter.py`）：必做。
2. **IK 滤波**（`ik_filter.py`）：必做，截止频率必须与 GRF 相同。
3. **sEMG 预处理**（`process_semg_jw.py`，新写）：必做。
4. **EMG 参考表**（`make_emg_reference.py`）：必做（若用 EMG tracking）。
5. **稳定区间判定**（`detect_steady_state.py`，新写）：必做，决定 30 s 窗。
6. **MocoInverse**（`moco_inverse.py`）：求解。
- ~~heel-strike 符号修正~~：JW 跳过（§0.1B）。
- ~~RRA / ID~~：MocoInverse 自带 residual actuators，不需要先跑 RRA；ID 仅用于核对 residual 需求（可选）。

---

## 2. sEMG 预处理（具体实现）

### 2.1 JW 导出结构与 SQR 不同
JW `*_sEMG.csv` 是**多设备合并导出**：5 行表头（`Devices` / `2000` / 设备块名 / 列名 / 单位），数据分两块——
- 块 1：370580 行 = 2000 Hz EMG（`Frame` 1..18529、`SubFrame` 0..19），**只取这一块**；
- 块 2：18529 行 = 100 Hz 设备块（力台等），丢弃。
列按**名字**取：`EMG7..EMG16` 在 0 基列 26..35；`EMG1..6` 全零；另有 `IM EMG1..16` 是重复副本（幅值一致），不用。
时间基 `t = ((Frame-1)*20 + SubFrame)/2000`，dt=0.0005 s，与 GRF/IK 对齐。
SQR 的 `process_semg.m` 用固定 12 列 `textscan` + `data(:,3:12)`，**读 JW 会读错信号**，故新写 Python 版 `process_semg_jw.py`。

### 2.2 处理链与参数（与 SQR 完全一致）
1. 裁剪到分析窗 ±1 s（留 solver burn-in）。
2. **去尖峰**：`|x| > 2.5 × P99.9` 判 artifact，掩膜 ±10 样本（5 ms）膨胀，每段用两侧各 50 个干净样本（25 ms）做 **PCHIP** 桥接。
3. **高通** butter(2, 10 Hz) + `filtfilt`（零相位，等效 4 阶）。
4. **整流** `abs`。
5. **低通包络** butter(2, 5 Hz) + `filtfilt`，clamp ≥0。
6. **归一化**：每通道除以**分析窗内包络最大值**（max-of-window，**不是 MVC**）。输出 0–1。
7. 输出保持 2000 Hz；重采样到 100 Hz 在第 4 步做。

### 2.3 通道→肌肉映射（与 SQR 相同，已用包络图验证）
| EMG | 肌肉 | Moco control |
|---|---|---|
| 7 / 12 | 股直肌 R/L | rect_fem_r/l |
| 8 / 13 | 股外侧肌 R/L | vasti_r/l |
| 9 / 14 | 股二头肌 R/L | hamstrings_r/l |
| 10 / 15 | 胫骨前肌 R/L | tib_ant_r/l |
| 11 / 16 | 腓肠肌 R/L | gastroc_r/l |

验证：`*_semg_envelope.png` 中 10 通道均呈清晰步态调制 burst，左右反相（vasti/gastroc 交替），符合生理。其余 8 块模型肌肉无 EMG，tracking 权重 0。
幅值量级与 SQR 同（几十 µV，Delsys 以 V 导出）；max 归一化后绝对幅值不影响结果。

### 2.4 EMG 参考表
`make_emg_reference.py --input *_semg_processed.csv --start (a0-1) --end (a1+1) --dt 0.01 --out *_emg_reference.sto`
→ `TimeSeriesTable`，列名 `/forceset/<muscle>`，值为**激励 excitation**（0–1），100 Hz。覆盖分析窗 ±1 s 以含分段 burn-in。

---

## 3. GRF 预处理（滤波 / 去抖 / 其它）

`grf_filter.py`，逐 belt：
1. **stance 掩膜** `vy > 10 N`；**去抖**：丢弃 <0.10 s 的 stance 段、填补 <0.04 s 的 swing 隙。
2. **力 vx,vy,vz + 力矩 tx,ty,tz**：butter(4, **15 Hz**) `filtfilt`，**按 stance 分段**、每段镜像_pad 0.30 s；swing 段严格置 0。
3. **Tukey 包络**：每 stance 端 25 ms 余弦 ramp，避免分段边界跳变。
4. **tx, tz 置 0**（舍入噪声 ~1e-5 N·m）；**ty（绕铅垂自由力矩）保留**。
5. **CoP px,py,pz**：同样 15 Hz 分段 filtfilt，但 **swing 段不动**（保留台面中心占位值）。
6. 重写 ExternalLoads XML，`<datafile>` 指向滤波后 mot。

要点：
- **滤波是必须的**：原始 2000 Hz 含冲击尖峰与高频噪声；MocoInverse 网格 2 ms 对应 500 Hz，但 IK/GRF 带宽必须一致且远低于网格频率，15 Hz 是 SQR 验证过的选择。
- **去抖是必须的**：否则阈值附近的噪声会被切成大量假 stance 段。
- **不需要**额外"去抖"CoP 或钳位：SQR 校验器 C1–C8 已覆盖；JW 沿用同法。
- **已知衰减**：15 Hz 后 |ty| 峰值仅保留 raw 的 32–45%（SQR 校验阈值 75% 是按 1000 Hz 调的）。JW 2000 Hz 原始 ty 高频成分更多。若后续发现 MocoInverse 平衡/残差异常，优先复查 ty 带宽，而不是盲目改 IK。
- **坡度组不做反旋转**（§0.1A）：数据保持跑台对齐系，模型重力已预旋转匹配；直接滤波即可。

---

## 4. IK 滤波

`ik_filter.py`：butter(4, **15 Hz**) `filtfilt`，全部非时间列，整段记录。输出 `*_ik_filtered.mot`。
**截止频率必须 = GRF 截止频率**，否则外载与运动学带宽不一致，MocoInverse 不自洽。
MocoInverse 内部再把旋转坐标 deg→rad（`prepare_kinematics_radians`）。
实测 walking1p0：去除 0.92% RMS，DC 漂移 1.6e-5（零相位、DC 单位增益，正常）。

---

## 5. 稳定区间判定（简单、可复现）

`detect_steady_state.py`，只用 GRF + IK：
1. 每 belt 取 stance 起始 = **足跟着地事件**（vy>20 N，去抖同 §3）。
2. 左右事件合并排序 → **步间期**序列。
3. **10 s 滑窗**（hop 1 s）内计算：步间期 **CV**、局部平均步间期相对全 trial 中位数的偏差、IK `pelvis_tx` 线性漂移。
4. 稳定判据：`CV ≤ 0.05` 且 `rate_dev ≤ 10%` 且 `|drift| ≤ 0.25 m/s` 且窗内 ≥8 步。
5. 取**最长连续稳定段**；分析窗 = 其**末 30 s**。
6. 输出 `steady_state_{summary.json, windows.csv, steps.csv, png}` 并重写 `notes.txt`。

为什么这样够用：跑台速度恒定 ⇒ 步态稳定 ⇔ 步间期规则 + 身体位置不漂移。CV 抓节奏抖动，rate_dev 抓速度切换，drift 抓人在带上前后移动。三者都通过即稳态。
实测区分力：JW_run 的 0–4.4 s 加速段与 95.4 s 后减速段被正确排除（见 `JW_run/steady_state.png`）；JW 各组稳态 CV 中位 0.02–0.04、drift 0.003–0.014 m/s，确实走得很稳，故稳定段接近全程。

---

## 6. 每组脚本就位与运行

`experiments/scripts/setup_group_pipeline.py` 已把 10 个脚本拷入每组：
`grf_filter.py, fix_grf_heel_strike.py, ik_filter.py, make_emg_reference.py, detect_steady_state.py, process_semg_jw.py, trial_paths.py, run_pipeline.py, moco_inverse.py, moco_bootstrap.py`。

`trial_paths.py` 从**脚本所在目录**自动发现本组文件（`*_grf.mot / *_grf.xml / *_ik.mot / *_simbody.osim`），`run_pipeline.py` 用这些覆盖各模块的路径全局量并把 `OUT_DIR` 设到本组 `moco_outputs/`，且 `os.chdir` 到本组 ⇒ **所有输出（含 opensim.log、solver sentinel、分段文件）都留在本组目录**。原 SQR 脚本不受影响。

运行（在组目录内）：
```
python run_pipeline.py                 # 前置：grf,ik,semg,emg
python run_pipeline.py --solve         # 前置 + MocoInverse（30 s 窗）
python run_pipeline.py --solve -- --mesh 0.002 --emg-weight 0.1
```
MocoInverse 关键参数（沿用 SQR 验证值）：`--mesh 0.002`（2 ms；1 ms 不收敛、20 ms  chatter）、`--segment 2.0 --overlap 0.5`（30 s ⇒ ~15 段，每段 ~1500 区间 ≈4 GB、3–5 min；整组 45–75 min）、`--emg-tracking *_emg_reference.sto --emg-weight 0.1`、residual 容量 400 N·m / 500 N。
**实测性能警告（本机）**：JW_walking1p0 在 2 ms 网格下单段（3 s、1500 区间）≈**29 min**（465 迭代），故 30 s/15 段 ≈**7 h**，比 README 估计慢 ~6–10 倍。后续 7 组生产求解的算力规划须按实测而非 README。慢因待查（RAM/swap、MUMPS 线程、EMG-gain 开销）。
QC：`check_solution_quality.py`（chatter alt/std >0.15 判不可用）+ `read_moco_results.m`。

**链路验证（已完成）**：JW_walking1p0 冒烟求解（0.3 s @ mesh 0.02）`solved=True, iterations=26, objective=172.26`，模型/IK 弧度化/GRF 窗/ExternalLoads/EMG 参考全部从组目录加载成功；residual 峰值 83 N / 46 N·m（远低于 500/400 容量），2/18 肌肉激活饱和。`.vtp` 几何缺失 warning 仅影响渲染，不影响求解。完整 30 s 求解（mesh 0.002、分段 2 s+0.5 overlap）随后执行。
**输入跨度审计（7/7 通过）**：各组 IK/GRF/EMG 参考均覆盖 `[a0-0.5, min(a1+0.5, 记录末端)]`（burn-in 钳位到数据末端）；EMG 参考左端留 1 s 余量。
**30 s 生产配置验证（JW_walking1p0，已完成）**：15/15 段 `solved=True`，window 154.689–184.689、mesh 0.002，总耗时 383.8 min，objective 38498.6，5219 迭代。`check_solution_quality.py` 判定 **OK：激活平滑、肌肉力可信**（worst activation chatter 0.016 ≪ 0.15；残差近界 0 个）。EMG tracking 生效（10 个 gain 解出于 [0.5,2]）。
**质量注意事项**：13/18 肌肉激活饱和（模型强度限制，SQR 亦 12/18）；14 个拼接刀口中 3 个（172.689 / 174.689 / 180.689 s）activation 跳变 0.14–0.28，属分段求解固有 seam 不连续，用肌肉力时避开这些刀口邻域或加大 overlap 重解。

---

## 7. 待确认 / 风险清单

1. ~~坡度组 GRF 坐标系~~：**已确认无需校正**（§0.1A），模型重力预旋转、取本组模型即可。
2. ~~sEMG 映射~~：**已确认与 SQR 一致**。
3. **ty 自由力矩衰减**（§3）：32–45% 峰值保留，留意平衡残差。
4. **residual 容量**：400/500 是按 SQR ID 需求调的；JW 跑步/坡度需求更大，若 residual 控制饱和（ pinned ±1）需上调。
5. **跑步 flight phase**：双足腾空段 GRF≈0，MocoInverse 靠 residual 维持平衡；检查该段 residual 是否异常大。
6. **模型强度**：SQR 模型 12/18 肌肉饱和；JW 同模型，跑步/坡度预期更差，必要时复核 Fmax。

---

## 8. 本次新增/修改文件

新增（`experiments/scripts/`）：
- `detect_steady_state.py` — 稳定区间判定 + notes.txt
- `process_semg_jw.py` — JW 多设备 sEMG 预处理
- `derotate_grf.py` — 坡度组 GRF 台面系倾斜的**纯诊断**报告（禁止改力数据）
- `trial_paths.py` — 组内路径自动发现
- `run_pipeline.py` — 组内流水线编排
- `setup_group_pipeline.py` — 拷贝脚本到各组

各组目录新增产物：`steady_state_*`、`*_grf_filtered.mot/.xml`、`*_ik_filtered.mot`、`*_semg_processed.csv`、`*_semg_artifacts.csv`、`*_semg_envelope.png`、`*_emg_reference.sto/.csv`、`grf_frame_report.json`（坡度组）、`notes.txt`（已填）、12 个脚本副本（含 `check_solution_quality.py`）。