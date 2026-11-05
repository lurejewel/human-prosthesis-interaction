# outputs 快速索引

完整说明见 [MocoInverse结果文件说明.md](MocoInverse结果文件说明.md)。

## 最终结果

```text
moco_inverse_solution_native_b_15_s13.sto
moco_inverse_summary_native_b_15_s13.txt
moco_inverse_s_used_native_b_15_s13.csv
final_1p5_tracking_report.csv
```

## 基线 / 对照

```text
A 形式 free-gain: moco_inverse_solution_A_none_1p5.sto
A 形式 no-gain 0.3 s: moco_inverse_solution_A_e2_compare.sto
B 形式初始 fixed-s: native_b_15.sto
插件 fixed-s 交叉验证: rev_fixs15_none.sto
0.3 s B joint（历史）: moco_inverse_solution_rev_joint_tight1e5.sto
```

## 绘图

```matlab
S = read_moco_results(struct('Tag', 'native_b_15_s13'));
```

B 形式画 `sEMG ./ s`；A 形式画 `sEMG .* g`。
