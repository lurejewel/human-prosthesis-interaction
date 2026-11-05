# MocoInverse 输出文件说明（sEMG tracking / B 形式）

本目录只保留当前有效的代表结果、最终解和后续对比需要的文件。旧的中间 `segments/`、
失败的 joint、未收敛的 random/initial-guess 结果和重复 scaled reference 已删除。

## 1. 先看结论

最终 1.5 s free-s 解：

```text
Tag: native_b_15_s13
J_total = 10.229322
mean cosθ = 0.900
mean Pearson r = 0.764
mean control-space nRMSE = 0.380
```

最终 s 和逐通道 tracking 指标：

```text
moco_inverse_s_used_native_b_15_s13.csv
final_1p5_tracking_report.csv
```

对应技术路线：

> 原生 fixed-s B 内层（reference = e/s，per-control weight = w s²）
> + 外层闭式 s 更新
> + `random`
> + 上一轮完整解作 guess。

A 形式只作为对比基线；joint/no-gain 不再作为主路线。

## 2. 文件命名规则

所有文件都有公共前缀：

```text
moco_inverse_<类型>_<tag>.<扩展名>
```

其中 `<tag>` 对应一次运行。例如：

```text
moco_inverse_solution_native_b_15_s13.sto
moco_inverse_summary_native_b_15_s13.txt
moco_inverse_s_used_native_b_15_s13.csv
```

没有 `_<tag>` 的文件是默认/旧结果；当前分析优先看带 `native_b_15_s13` 的文件。

## 3. 输入文件

| 文件 | 含义 |
|---|---|
| `emg_excitation_199p5_204p0.sto` | 0.3 s 验证窗使用的 sEMG reference（列名 `/forceset/<muscle>`） |
| `grf_window_200_201p5.mot` | 1.5 s 窗口的 GRF 输入 |
| `grf_window_200_201p5.xml` | 1.5 s 窗口的 ExternalLoads 配置 |
| `ik_radians.mot` | 处理后的 IK 运动学（弧度） |
| `emg_native_b_scaled_native_b_15_s13.sto` | B 形式的 scaled reference：`e_i / s_i`（最终 s13） |

`sEMG reference` 的时间范围必须覆盖 MocoInverse 窗口；`grf`/`ik` 也必须覆盖同一窗口。

## 4. 解文件 `.sto`

### `moco_inverse_solution_<tag>.sto`

完整的 Moco 轨迹，列包括：

- 第 1 列：`time`；
- 18 个状态：`/forceset/<muscle>/activation` + 关节坐标/速度；
- 37 个控制：18 个肌肉 excitation + 6 个骨盆 residual + 13 个 reserve。

这是最完整的轨迹，`read_moco_results` 默认读它。

### `moco_inverse_activations_<tag>.sto`

只含 18 块肌肉的 activation（state），方便和 EMG 对比。

### `moco_inverse_excitations_<tag>.sto`

只含 18 块肌肉的 excitation（control），是 `MocoControlTrackingGoal` 实际跟踪的量。

### 当前保留的 `solution` 文件

| 文件 | 用途 |
|---|---|
| `moco_inverse_solution_native_b_15_s13.sto` | **最终 B 解**（1.5 s，s13） |
| `native_b_15_s13.sto` | 最终解的分段轨迹副本（65 列，含 s 参数时用于诊断） |
| `native_b_15.sto` | B0：初始 fixed-s（s0） |
| `rev_fixs15_none.sto` | 插件 fixed-s 1.5 s（用于交叉验证） |
| `moco_inverse_solution_native_b_15_guess_fixed.sto` | P0 修复 guess 后的自重启结果 |
| `moco_inverse_solution_A_none_1p5.sto` | A 形式 free-gain 1.5 s 基线 |
| `moco_inverse_solution_A_e2_compare.sto` | A 形式 no-gain（scale=1）0.3 s 对照 |
| `moco_inverse_solution_rev_joint_tight1e5.sto` | 0.3 s B joint（历史对照，不作为主路线） |
| `moco_inverse_solution_emg200_w20_nogain.sto` | A no-gain 200–203.5 s 主结果（不同窗口，不能直接和 1.5 s B 比） |

## 5. summary `.txt`

`moco_inverse_summary_<tag>.txt` 包含：

- `window`、`mesh`、`segments`；
- `solved=True/False`；
- `objective`；
- `solver_iterations`；
- `solver_duration_s` / `wall_time_s`；
- `activation_saturated_muscles`；
- 分段 wall time。

当前关键 summary：

```text
moco_inverse_summary_native_b_15_s13.txt   # 最终 B，J=10.229322
moco_inverse_summary_native_b_15.txt       # B0，J=10.653837
moco_inverse_summary_native_b_15_guess_fixed.txt  # 同一 s 自重启，228 iter/53.9 s
moco_inverse_summary_final_mesh001.txt     # 细网格筛查，J=10.21262
moco_inverse_summary_alt_s1_0/1/2.txt      # 从 s=1 出发的初值稳健性检查
```

## 6. `moco_inverse_gains_<tag>.csv`

A 形式的每块肌肉增益 `g`：

```text
segment,muscle,gain
```

用于 A 形式绘图：`sEMG .* g`。

B 形式的 s **不放在这里**，而是放在：

```text
moco_inverse_s_used_<tag>.csv
```

格式：

```text
segment,muscle,s
```

B 形式绘图用 `sEMG ./ s`，不是 `sEMG .* s`。

当前 B 文件：

```text
moco_inverse_s_used_native_b_15.csv      # s0
moco_inverse_s_used_native_b_15_s13.csv  # 最终 s13
```

## 7. `moco_inverse_activation_report_<tag>.csv`

逐肌肉激活峰值/均值：

```text
muscle,peak_activation,mean_activation
```

用于快速检查哪块肌肉饱和；但它只是 `activation` state，不是 excitation control。
报告中要区分 activation 与 excitation。

## 8. `moco_inverse_residuals_<tag>.csv`

残差/储备作动器的峰值报告：

```text
actuator,peak_control,optimal_force,peak_force
```

`peak_control` 是控制峰值，`peak_force = peak_control * optimal_force`。
当前最终解的 reserves 重点看：

```text
moco_inverse_residuals_native_b_15_s13.csv
```

注意：非矢状 reserves 的权重只有 1e-3，hip adduction/rotation 数值偏大是权重设置的结果；
不能只用 EMG 拟合好坏判断这些 reserve 是否合理。

## 9. `final_1p5_tracking_report.csv`

最终 1.5 s 逐通道 A/B 对比：

```text
muscle,A_cos,A_pearson,A_nRMSE,
B0_cos,B0_pearson,B0_nRMSE,
B13_cos,B13_pearson,B13_nRMSE,s13
```

定义：

- `cos`：模型 control 与目标 `e/s`（B）或 `g e`（A）的 cosine；
- `pearson`：同一对的 Pearson r；
- `nRMSE`：control-space 归一化 RMSE。

A 是 A 形式 free-gain，不是 scale=1；scale=1 的 A 是 `A_e2_compare`（0.3 s，窗口不同）。

## 10. 如何用 `read_moco_results` 绘图

MATLAB：

```matlab
S = read_moco_results(struct('Tag', 'native_b_15_s13'));
```

它会：

- 读 `moco_inverse_solution_native_b_15_s13.sto`；
- 读 `moco_inverse_s_used_native_b_15_s13.csv`；
- 设置 `S.scaleKind = "B_s"`；
- 图 1 绿线画 `sEMG ./ s`，标题标 `s=...`。

A 形式没有 `s_used` 文件时，退回旧逻辑：

```text
S.scaleKind = "A_gain"
绿线 = sEMG .* g
标题 = g=...
```

## 11. 结果解读时必须带上的限制

- B 的 s 是**这个模型、这个窗口、这个权重下的有效幅值校准**，不是生理增益；
- s 弱可辨识：最终 \(R_s\approx3.2\%\)，s13 后 `s_next` 仍在漂移；
- B13 的 mean cosθ 0.900、Pearson 0.764，明显优于 A 的 0.816/0.626；
- 但 vasti 仍弱（r 0.37/0.47），reserves（hip adduction/rotation、knee、pelvis residual）仍不小；
- 网格 0.02→0.01 只验证该 s 下的轨迹，不验证 s 本身跨网格稳定；
- A/B 的 objective 不能直接比较：\((e-su)^2=s^2(u-ge)^2\)。

## 12. 当前不保留的内容

- 旧 `segments/`；
- 未收敛的 free-s joint；
- random/initial-guess 的旧 5.0/5.8/6.3 结果；
- 重复的 `emg_native_b_scaled_*`；
- 旧 plots/logs。

如果需要复现某个历史轮次，按 `docs/sEMG_tracking_工作总结.md` 的配置重新跑；不要从已删除的中间文件恢复。

## 13. 30 s 右腿结果（210–240 s，per-muscle s，跨段共享）

最终候选：

```text
moco_inverse_solution_B30_right_final.sto
moco_inverse_summary_B30_right_final.txt
moco_inverse_s_used_B30_right_final.csv
B30_right_tracking_report_final.csv
B30_right_outer_summary.csv
```

A no-gain 基线：

```text
moco_inverse_solution_A30_right_merged.sto
moco_inverse_summary_A30_right_merged.txt
```

B 稳定性检查与最终解：

```text
moco_inverse_solution_B30_right_final.sto       # r7 / s6 final
moco_inverse_solution_B30_right_r0_recheck_merged.sto
moco_inverse_solution_B30_right_r6_recheck_merged.sto
B30_right_goal_components.csv
B30_right_u_stability_check.csv
B30_right_outer_summary.csv
```

- `B30` 只跟踪右腿 5 块；左腿保留但无 EMG tracking；
- s 是 5 个 per-muscle 值，30 s 所有段共享；
- `moco_inverse_s_used_B30_right_final.csv` 有 5 行；
- 最终 \(R_s\approx3.27\%\)，s 未收敛；r3 因脚本重启重复 s2，已排除；
- `B30_right_goal_components.csv` 的 effort/tracking 是 3 s solve 窗口和，
  **不是** 30 s retained-window 目标；
- `B30_right_u_stability_check.csv` 是 r0/r6/r7 原始 u 的 RMS、相关和 `s·u` 对比。
- 读图时用：
  ```matlab
  S = read_moco_results(struct('Tag','B30_right_final'));
  ```
  默认 `RightOnly=true`，肌肉图只画右腿 9 块（含 4 块无 sEMG 的
  bifemsh_r/glut_max_r/iliopsoas_r/soleus_r），GRF 只画右腿。
  这些无 sEMG 的肌肉标题显示 `s=n/a`、没有绿色 sEMG 线。
  当前 `--sagittal-cost` 解没有 Moco 坐标状态；坐标图会用 IK reference 占位，
  并标注 `IK reference only`，不会再出现空白 figure。
