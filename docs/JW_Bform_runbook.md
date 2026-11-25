# JW 各组 B-form 内外交替 MocoInverse 运行手册

适用脚本：各组目录下的 `run_bform_outer.py` 与 `moco_inverse.py`。**只支持 B-form；A-form 已移除。**

## 1. 求解流程

对每块有 sEMG 的肌肉 \(m\)（左右共 10 块）：

- 内层固定 \(s_m\)，MocoInverse 使用参考 \(e_m/s_m\) 和控制权重 \(w s_m^2\)，即最小化 \(\sum w(e_m-s_m u_m)^2\)。
- 内层收敛后，外层由轨迹闭式更新
  \[
  s_m \leftarrow \operatorname{clip}\left(\frac{\sum w e_m u_m}{\sum w u_m^2},0.5,2\right).
  \]
- 默认 `--s-domain solve`：在包含 burn-in 的内层求解区间上更新，与内层目标积分域一致。
- 每一代将 30 s 分成 4 个 7.5 s chunk，以独立子进程并行；每 chunk 内仍按 2 s retained + 两侧 0.5 s burn-in 分段。
- 每代合并 4 个 chunk，并作为下一代的 initial guess。

### 为什么需要第 0 代

是的，**先跑一代形成初值**。第 0 代默认 10 块肌肉均为 `s=1`，不使用历史轨迹；成功后得到首条 30 s 拼接轨迹。第 1 代起才使用“上一代轨迹作初值 + 上一代闭式更新的 s”。

**不要从之前的 A-form 结果开始；请重新跑 B-form 第 0 代。**

## 2. 拷贝到另一台电脑

直接拷贝整个 `JW_<组名>` 文件夹即可。必须包含：

- `JW_simbody.osim`（各组自己的模型，坡度组重力方向不同）；
- `*_ik_filtered.mot`、`*_grf_filtered.mot/.xml`、`*_emg_reference.sto`；
- `steady_state_summary.json`（提供 30 s 窗口）；
- `run_bform_outer.py`、`moco_inverse.py`、`moco_bootstrap.py`。

预处理已经完成，不需要重跑。若原始数据改变，先在该组运行 `python run_pipeline.py`。

## 3. 环境

使用 conda 的 OpenSim Python 3.12，不要用系统 Python 3.14：

```
conda create -n opensim_scripting python=3.12
conda activate opensim_scripting
conda install -c opensim-org -c conda-forge opensim
conda install numpy scipy matplotlib
```

`moco_bootstrap.py` 会从当前解释器旁查找 CasADi/IPOPT 插件；不需要单独 `pip install casadi`。

资源建议：

- 4 路并行：建议 16 个逻辑核、至少 32 GB RAM；
- 内存较小或出现 swap：`--workers 2`；
- 不要同时以 4 路并行跑多个组，避免 CPU/内存超订。

## 4. 推荐运行方式：先 bootstrap，再续跑

进入某一组目录并激活 conda 环境。

### 4.1 只打印命令，不求解

```
python run_bform_outer.py --workers 4 --dry-run
```

应看到 4 个 chunk，均含 `--emg-native-b-fixed`、10 个 `s=1`、独立 `--out-dir` 和 `--seg-dir`。

### 4.2 第 0 代：形成初值

```
python run_bform_outer.py --workers 4 --bootstrap-only
```

成功后关键文件：

- `bform_outputs/BFORM_u000_merged.sto`：首条 30 s 轨迹；
- `bform_outputs/checkpoint.json`：记录下一代应使用的 `s_next`；
- `bform_outputs/progress.csv`：第 0 代指标；
- `bform_outputs/logs/BFORM_u000_c*.log`：各 chunk 日志。

### 4.3 从 checkpoint 做外层迭代

```
python run_bform_outer.py --resume bform_outputs/checkpoint.json --workers 4 --max-updates 40
```

默认连续 3 代满足以下条件即停止：

- `R_s < 1%`；
- `u_rel < 1%`；
- `act_rel < 2%`；
- `rJ < 0.1%`。

中途停止或断电后，重复同一条 `--resume` 命令即可；每代成功后才推进 checkpoint。

## 5. 也可以一次跑完

```
python run_bform_outer.py --workers 4 --max-updates 40
```

该命令自动完成第 0 代和后续外层迭代，不需要手动提供初值。

常用选项：

- `--workers 2`：降低并行度；
- `--max-updates 10`：限制本次外层更新数；
- `--init-s "1,1,1,1,1,1,1,1,1,1"`：修改第 0 代的 10 个 s；
- `--init-guess <B-form轨迹.sto>`：仅给第 0 代提供自定义初值；
- `--mesh 0.002 --tol 1e-4 --max-iter 5000`：默认生产数值设置；
- `--emg-weight 20`：默认与 SQR B-form continuation 一致。

## 6. 最终结果与检查

正常结束：

- `bform_outputs/BFORM_final.sto`：最终 30 s 轨迹；
- `bform_outputs/BFORM_final_s_used.csv`：产生该轨迹时实际固定的 10 个 s；
- `bform_outputs/progress.csv`：每代的 s、R_s、rJ、轨迹变化及收敛状态；
- `bform_outputs/chunks/`：每代、每 chunk 独立输出；
- `bform_outputs/logs/`：运行日志。

可对最终结果运行（时间换成该组 `notes.txt` 中的分析窗）：

```
python check_solution_quality.py bform_outputs/BFORM_final.sto --window 154.689 184.689
```

生产建议仍使用 2 ms 网格；20 ms 只用于快速测试，不能作为肌肉力结果。并行减少墙钟时间，但不减少总迭代计算量，也不会改变 B-form 的数学定义。