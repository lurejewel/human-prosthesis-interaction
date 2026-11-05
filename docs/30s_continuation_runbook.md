# 30 s 右腿 B 形式续跑手册（移植到另一台机器）

本手册对应：

```text
run_B30_continue.py
experiments/moco/moco_inverse.py
experiments/moco/outputs/
```

目标不是“再跑固定 13/40 轮”，而是：

> 在**目标、权重、bounds、求解精度不变**的前提下，检验 s 与原始 excitation/activation
> 是否达到数值稳定；并明确区分内部收敛、边界收敛、预算耗尽。

---

## 0. 先确认一个关键口径问题

当前 v1 的 6 次有效更新里：

- 内层优化的是 **3 s solve-window**（含两侧 0.5 s overlap/burn-in）的目标之和；
- 外层 s 却是用 **30 s retained（拼接后 2 s 段）** 的 controls 更新。

**两者不是同一个积分域。** 所以 v1 不能解释为严格的交替最小化；它只能算一个
“retained-domain heuristic outer update”。

`run_B30_continue.py` 因此提供两种模式：

```text
--s-domain solve     推荐；s 用 3 s solve-window 更新，与内层目标一致
--s-domain retained  旧版；s 用 30 s retained 更新，只作对照
```

**不要把两种模式的 s 序列混在一起。** 新长跑默认用 `solve`，并使用新的
`--tag-prefix` 和新的 checkpoint。

注意：r7 的 per-segment 文件已经被清理，所以第 1 次 `solve` 更新会先用 retained
估计作为 bootstrap（程序会打印 `s_hat bootstrap=retained`）；从第 2 次开始才完全用
3 s solve-window。这个 bootstrap 是显式标注的，不是偷偷混口径。

---

## 1. 需要带到新机器的文件

### 1.1 代码

```text
run_B30_continue.py
experiments/moco/moco_inverse.py
experiments/moco/moco_bootstrap.py
experiments/moco/read_moco_results.m          # 可选，画图用
```

### 1.2 数据

整个目录：

```text
experiments/data/SQR_walking/
```

至少包含：

```text
SQR_simbody.osim
level_walking_ik_filtered.mot
level_walking_grf_fixed.xml  或  level_walking_grf_filtered.xml
对应的 GRF .mot
```

### 1.3 续跑起点

```text
experiments/moco/outputs/emg_excitation_210_240_right.sto
experiments/moco/outputs/moco_inverse_solution_B30_right_final.sto
```

建议同时带上参考文件：

```text
moco_inverse_summary_B30_right_final.txt
moco_inverse_s_used_B30_right_final.csv
B30_right_outer_summary.csv
B30_right_goal_components.csv
B30_right_u_stability_check.csv
```

`moco_inverse_solution_B30_right_final.sto` 是 r7 的 30 s 轨迹；`s_used` 是：

```text
1.662868, 1.493015, 0.851994, 1.177731, 1.068418
```

### 1.4 Python 环境

```text
Python 3.12.x
OpenSim 4.6（与当前一致）
Moco 4.6
numpy / scipy
CasADi IPOPT（OpenSim conda 包自带）
```

检查：

```powershell
python -c "import opensim as osim; print(osim.__version__)"
python -c "import opensim as osim; print(osim.MocoCasADiSolver.isAvailable())"
```

第二个必须输出 `True`。

---

## 2. 第一次在新机器上运行

### 2.1 先做静态检查

```powershell
python -m py_compile run_B30_continue.py
python -m py_compile experiments\moco\moco_inverse.py
```

### 2.2 dry-run：只打印命令，不求解

```powershell
python .\run_B30_continue.py `
  --dry-run `
  --s-domain solve `
  --max-updates 1
```

应看到：

```text
B30 continuation | domain=solve | workers=4
start s=[...]
s_hat bootstrap=retained (first solve-domain update)
DRY: ... moco_inverse.py --start 210.0 ...
```

如果 dry-run 都失败，不要开始长跑，先按第 6 节查路径/环境。

### 2.3 第一次真跑：只跑 1 个 update

先验证 4 个 chunk 都能正常结束、能写 checkpoint：

```powershell
python .\run_B30_continue.py `
  --s-domain solve `
  --max-updates 1 `
  --max-total-updates 40 `
  --workers 4 `
  --init-s "1.662868,1.493015,0.851994,1.177731,1.068418" `
  --init-guess "moco_inverse_solution_B30_right_final.sto" `
  --tag-prefix B30_cont
```

成功后应生成：

```text
experiments/moco/outputs/B30_cont_checkpoint.json
outputs/B30_cont_progress.csv
outputs/B30_cont_u001_c0..c3.sto
outputs/B30_cont_u001_merged.sto
outputs/continuous_logs/B30_cont_u001_c*.log
outputs/segments_B30_cont_u001_c*/seg_XX.sto
```

检查 `B30_cont_progress.csv` 有 1 行，且 `solver_ok=True`。

### 2.4 正式长跑：从 checkpoint 续跑

```powershell
python .\run_B30_continue.py `
  --resume "experiments\moco\outputs\B30_cont_checkpoint.json" `
  --s-domain solve `
  --max-updates 10 `
  --max-total-updates 40 `
  --workers 4
```

- `--max-updates 10`：本次最多 10 个 update；
- `--max-total-updates 40`：包含之前的累计最多 40 个；
- 到 3 次连续稳定会自动停；
- 没到稳定就到预算，程序会正常结束，checkpoint 保留。

之后还可以继续：

```powershell
python .\run_B30_continue.py `
  --resume "experiments\moco\outputs\B30_cont_checkpoint.json" `
  --s-domain solve `
  --max-updates 10 `
  --max-total-updates 40 `
  --workers 4
```

**不要**用 `--init-s`/`--init-guess` 重新开始同一序列；resume 才用 checkpoint。

---

## 3. PyCharm 运行方式

1. 打开整个仓库目录；
2. `File > Settings > Project > Python Interpreter` 选择新机器的 `opensim_scripting` 环境；
3. `Run > Edit Configurations > + > Python`：
   - Script path：`run_B30_continue.py`
   - Parameters：
     ```text
     --s-domain solve --max-updates 1 --max-total-updates 40 --workers 4
     ```
   - Working directory：仓库根目录；
4. 先跑 dry-run，再跑 1 update，再改成 `--max-updates 10` 长跑。

PyCharm 控制台会显示 wrapper 的 `update`、`s`、`J`、`R_s`、`u_rel` 等行。  
子进程的 IPOPT 迭代输出在：

```text
outputs/continuous_logs/B30_cont_uXXX_cY.log
```

不要在 PyCharm 里同时开多个相同 tag 的运行；会抢同一个 `outputs` 和 `segments_*`。

---

## 4. PowerShell 运行方式

### 4.1 前台 + 同时看日志

```powershell
$env:OMP_NUM_THREADS = "2"
$env:MKL_NUM_THREADS = "2"

python .\run_B30_continue.py `
  --s-domain solve `
  --max-updates 10 `
  --max-total-updates 40 `
  --workers 4 `
  --resume "experiments\moco\outputs\B30_cont_checkpoint.json" `
  2>&1 | Tee-Object -FilePath ".\b30_cont_console.log"
```

如果 PowerShell 管道缓冲导致输出不及时，改用重定向：

```powershell
python .\run_B30_continue.py `
  --s-domain solve --max-updates 10 --max-total-updates 40 --workers 4 `
  --resume "experiments\moco\outputs\B30_cont_checkpoint.json" `
  *> ".\b30_cont_console.log"
```

另开一个窗口实时看：

```powershell
Get-Content .\b30_cont_console.log -Wait -Tail 40
```

### 4.2 后台运行

```powershell
$args = @(
  ".\run_B30_continue.py",
  "--s-domain","solve",
  "--max-updates","10",
  "--max-total-updates","40",
  "--workers","4",
  "--resume","experiments\moco\outputs\B30_cont_checkpoint.json"
)
Start-Process -FilePath "python" `
  -ArgumentList $args `
  -WorkingDirectory (Get-Location) `
  -RedirectStandardOutput ".\b30_cont_console.log" `
  -RedirectStandardError ".\b30_cont_err.log"
```

之后看 `b30_cont_console.log`、`b30_cont_err.log`。

### 4.3 并发/内存

- `--workers 4` 是默认；如果内存不足改为 2 或 1；
- 每个 worker 内不要再开大量 BLAS 线程，建议 `OMP_NUM_THREADS=2`；
- 不要同时跑多个 `run_B30_continue.py` 用同一个 `--tag-prefix`。

---

## 5. 输出怎么看

### 5.1 控制台每轮关键行

```text
--- update 3 / total 3 ---
s_used = [...]
s_hat  = [...]
J=... effort=... tracking=... blocks=16
mean cos=... r=... nRMSE=...
R_s=...% rJ=...% G_s=... u_rel=...% act_rel=...%
stable=False (0/3) bound=[] elapsed=...
```

解释：

| 字段 | 含义 |
|---|---|
| `update_id` | generation/有效更新编号 |
| `s_used` | 本轮内层实际用的 5 个 s |
| `s_hat` / `s_hat_new` | 本轮轨迹上重新估计的 5 个 s |
| `J` | 4 个 chunk 的 solve-window objective 之和 |
| `effort` / `tracking` | 日志分项和；3 s solve-window，不是 30 s retained |
| `R_s` | `max_i |s_hat_i - s_used_i| / s_used_i` |
| `rJ` | 与上一有效轮的 `J` 相对变化 |
| `G_s` | 固定当前轨迹，仅把 `s_used` 换到 `s_hat_new` 的 tracking 下降 |
| `u_rel` | 5 块 tracked 肌肉原始 excitation 的最大相对 L2 变化 |
| `act_rel` | activation 的最大相对 L2 变化 |
| `stable_count` | 连续满足停止条件的次数，满 3 次停止 |
| `bound` | 接近 0.5/2.0 的 s 名称 |
| `blocks` | 成功解析到目标分项的 segment 数；应为 16 |

`G_s` 是当前轨迹上的半步收益，不是全局最优距离上界；不要拿 `G_s/203` 判断“已经很小”。

### 5.2 `B30_cont_progress.csv`

每行一个 update，列为：

```text
update_id,timestamp,J,effort,tracking,mean_cos,mean_r,mean_nRMSE,
s0..s4,s_hat0..s_hat4,R_s,rJ,G_s,u_rel,act_rel,
stable_count,bound,solver_ok,elapsed_s
```

用 Excel/Python 画：

- `R_s` vs update；
- `u_rel`、`act_rel` vs update；
- `J`、`tracking`、`effort` vs update；
- `s_i` vs update。

### 5.3 `B30_cont_checkpoint.json`

最重要的恢复文件，包含：

```text
domain, update_id, total_updates,
s_current, s_hat, solution, guess,
J, effort, tracking, R_s, rJ, G_s, u_rel, act_rel,
stable_count, bound, saturation, timestamp
```

- `s_current` 是最后用于求解的 s；
- `s_hat` 是下一次候选 s；
- `solution` 是当前 30 s 轨迹；
- 恢复时用 `--resume`。

### 5.4 `continuous_logs/*.log`

每个 chunk 一个日志，包含：

- `fixed s: ...`
- 每个 segment 的 `solved=True/False`、`objective=`、`wall=`
- 目标分项：
  ```text
  excitation_effort: ...
  emg_tracking: ...
  ```
- IPOPT 迭代。

如果 chunk 失败，先看这个日志的最后 100 行。

### 5.5 结束时看哪里

程序会在以下情况结束：

1. `stable_count >= 3`：连续 3 次满足：
   - `R_s < 1%`
   - `u_rel < 1%`
   - `act_rel < 2%`
   - `rJ < 0.1%`
2. 达到 `--max-total-updates`；
3. 某个 chunk 返回非 0，程序停止且**不推进 checkpoint**。

结束后：

```text
outputs/B30_cont_final.sto
outputs/B30_cont_final_s_used.csv
experiments/moco/outputs/B30_cont_checkpoint.json
outputs/B30_cont_progress.csv
```

`final.sto` 是最后成功 merge 的轨迹；`final_s_used.csv` 是求解该轨迹实际使用的 s。
如果 checkpoint 的 `stable_count<3`，应标注为“预算内未达到 3 次稳定”，不要标成收敛。

---

## 6. 报错自检表

| 现象 | 先检查 | 处理 |
|---|---|---|
| `ModuleNotFoundError: opensim` | 当前 Python 是否是 `opensim_scripting` | 在 PyCharm 选对 interpreter；PowerShell 用完整路径 |
| `MocoCasADiSolver.isAvailable()` 为 False | CasADi/IPOPT 插件路径 | 让 `moco_bootstrap.ensure_casadi_plugins()` 执行；检查 `Library\bin\casadi_nlpsol_ipopt.dll` |
| `FileNotFoundError` 模型/IK/GRF | `experiments/data/SQR_walking/` 是否完整 | 补齐数据；确认 XML 名字是 `fixed` 还是 `filtered` |
| `FileNotFoundError` reference | `emg_excitation_210_240_right.sto` 是否在 outputs | 从原机器复制；确认 5 列、时间覆盖 209.5–240 |
| `Bad allocation` / 内存不足 | `--workers` 是否太大 | 改 2 或 1；关闭其他程序；不要提高 mesh |
| 某个 chunk `rc=1` | `continuous_logs/<tag>.log` 最后 100 行 | 先单独重跑该 chunk；检查该段 GRF/IK/时间范围 |
| `IPOPT ... restoration failed` | 该段数据/初值/迭代上限 | 单独重跑；可提高 `--max-iter`，不要改目标 |
| `seg_XX.sto` 找不到 | `--seg-dir` 是否唯一、目录是否可写 | 检查 `outputs/segments_<tag>/`；不要和别的 run 共用 |
| `s_hat` 是 NaN | reference 列名/时间、controls 是否为 5 块 | 检查 `emg_excitation` 列名 `/forceset/...`；检查 `--emg-leg right` |
| `resume` 后像重复了一轮 | checkpoint 的 `update_id`/`s_hat` | 只用最新 checkpoint；不要手工改 `s_current` |
| 跑到一半没输出 | 控制台缓冲 / 子进程日志 | 看 `continuous_logs`；用 `Get-Content -Wait`；不要重复启动 |
| `R_s` 不降 | 不是程序错误 | 看 `u_rel/act_rel/bound`；按第 7 节判断 |

---

## 7. 什么时候算结束/怎么解释

### 情况 A：`stable_count=3`

说明在当前目标、权重、bounds、求解精度下达到数值稳定。  
可以写：

> 当前设置下获得数值收敛解。

但不能自动写“s 是唯一生理增益”；需要其他初始化/profile 才能谈可辨识性。

### 情况 B：部分 s 在 `bound`

`bound` 列非空，且输出稳定：

- 可能是当前约束下的边界最优；
- 报告中明确哪些 s 活跃在边界；
- 不能把边界 s 解释为数据独立确定的生理参数；
- 建议对关键输出做 bounds 敏感性。

### 情况 C：预算用完仍未稳定

`stable_count<3` 且 `max_total_updates` 到：

> 当前预算/求解设置下尺度方向仍未稳定；需要进一步处理慢收敛或幅值约束。

这时再考虑独立标定、先验、敏感性或算法改进；不要现在混入新约束。

---

## 8. 需要发回给我哪些结果

如果报错或跑完，请发：

1. 完整命令（包括 `--s-domain`、`--resume`、`--workers`、`--tag-prefix`）；
2. Python/OpenSim/Moco 版本；
3. `B30_cont_progress.csv`；
4. `B30_cont_checkpoint.json`；
5. `continuous_logs/` 里对应 update 的 `.log`（失败时至少失败 chunk 的 log）；
6. `B30_cont_final_s_used.csv`；
7. `B30_cont_final.sto` 的 summary（或 `moco_inverse_summary_*`）；
8. 报错窗口的完整 traceback。

不要只发一行 `objective`；没有 `s_used`、`s_hat`、`R_s`、`u_rel`、`act_rel` 和 log，无法判断是收敛、边界还是求解失败。

---

## 9. 当前已确认的检查结果

- r3 是重复 s2（脚本重启 bug），已从有效更新中剔除；
- 有效更新 6 次，最终 r7 `J=203.028772`，`R_s=3.27%`；
- 目标分项（3 s solve-window 和）：
  - effort 169.618→169.517（−0.06%）
  - tracking 35.213→33.512（−4.8%）
- 原始 u：
  - r6→r7：RMS 变 0.3–3%，corr≈1；
  - r0→r7：rect_fem_r −18%、gastroc_r −17%、hamstrings_r −10%、tib_ant_r −4%、vasti_r +4%；
  - `s·u` 只变 1–4%，形状相关 >0.995；
- 右腿未 tracking 肌肉 RMS 变化 ±3% 以内，reserves 峰值基本不变；
- 所以：**EMG 拟合稳定 ≠ 原始 u 稳定**；累计幅值变化主要由 s 补偿。

这些是续跑前的基线，不是最终收敛结论。

