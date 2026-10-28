# Ni–14 at.%Al、1000 K：二维熟化结果

计算完成：模型时间 `t*=0–60000`，1147 个接受时间步，231 个保存帧。采用 Thermo-Calc Educational 2026b / NIDEMO 2.0 的局部二次自由能及参考界面能、扩散系数；时间轴尚未实验标定。

## 组织演化

![熟化动画](coarsening_preview.gif)

![相场及 Al 浓度](evolution.png)

## 颗粒等效直径

![颗粒直径与求解时间步](particle_diameter_vs_step.png)

| 颗粒 | 初始直径 / nm | 最终直径 / nm | 结果 |
| --- | ---: | ---: | --- |
| 1，左下 | 18.47 | ≈0 | 溶解 |
| 2，右下 | 19.47 | 18.46 | 净收缩 |
| 3，左上 | 20.66 | 23.85 | 长大 |
| 4，右上 | 21.66 | 25.69 | 长大 |

等效直径 `d=2√(∫h(η)dA/π)`，其中 `h(η)=η³(6η²−15η+10)`。四个颗粒保持分离，因此按初始象限跟踪；此定义不适用于颗粒跨区移动或并合后的普遍情况。平滑界面使积分直径与初始几何直径略有差别。

[particle_diameter_vs_step.csv](particle_diameter_vs_step.csv) 将每个保存帧的模型时间与 [diagnostics.csv](diagnostics.csv) 的接受步逐一匹配，最大时间匹配差约 `5×10⁻⁸`。求解步编号以初始状态为 0；保存帧数不是时间步数。记录中的半径乘 2 得到直径，未对求解结果进行拟合或补造。

颗粒 1 首次在第 1040 步的保存帧（`t*≈51141.63`）中不再含 `η≥0.5` 的区域，最终连通颗粒数从 4 降到 3。此阈值消失判据和面积积分直径为零是两个不同的量；最终直径约 `5.7×10⁻⁹ nm`，表中记为近似零。

各象限 Al 积分显示小颗粒区域失去溶质、较大颗粒区域获得溶质，空间证据见 [diffusion_evidence.png](diffusion_evidence.png)。

## 校验与文件

- 总体 Al 摩尔分数维持 0.14，节点积分最大偏差约 `1.47×10⁻¹⁴`。
- 总约化自由能从 `35.91768` 降至 `30.24406`，接受步中未出现上升。
- 原计算的 Exodus 全场回读通过，时间严格递增，场值有限且在物理范围内，见 [validation.json](validation.json)。
- [particle_history.csv](particle_history.csv) 保存半径、面积、颗粒数、区域 Al 积分；[diameter_summary.json](diameter_summary.json) 保存直径统计和源 CSV 哈希。
- [diagnostics.png](diagnostics.png) 展示原半径曲线、守恒和能量；[model_parameters.json](model_parameters.json) 保存模型参数。
- 主输入见 [coarsening.i](../../examples/ni_al_coarsening/coarsening.i)，方程和适用边界见 [MODEL.md](../../examples/ni_al_coarsening/MODEL.md)。

`coarsening.csv` 是原始追加式记录；生成 `diagnostics.csv` 时排除了 92 条被检查点恢复替代的记录。原 Exodus 的一个超过最终时刻的过期尾部记录被剔除，保留的场值未修改，操作记录见 [restart_tail_cleanup.json](restart_tail_cleanup.json)。本仓库包含清理前的 CSV 和清理后的诊断表；完整 Exodus、旧版本及检查点保留在原计算机器，不随仓库上传。

本地最终场文件 `coarsening.e` 大小为 139,242,316 字节，哈希见 [exodus_sha256.json](exodus_sha256.json)。克隆仓库不会取得该场文件；按照根目录 README 运行主输入可重新生成场结果。原运行含检查点恢复和步长控制调整，干净重算不保证每个时间步完全一致。

这些检查验证本次演示的数值一致性和熟化现象，不能替代系统收敛研究、三维建模或实验标定。
