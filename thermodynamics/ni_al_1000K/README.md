# Ni–Al 在 1000 K 的热力学数据

本机实际计算：**Thermo-Calc Educational 2026b，NIDEMO 2.0**。导出日期：2026-10-19。
条件为 `T=1000 K`、`P=101325 Pa`、`N=1 mol`；目标总体成分 `x_Al=0.14`、`x_Ni=0.86`。
这里的成分均为原子/摩尔分数，不是质量分数。本结果描述体相化学平衡，未计入相干弹性能和曲率效应。

## 目标成分下的平衡

| 量 | γ | γ′ |
| --- | ---: | ---: |
| 平衡 Al 含量 / at.% | 12.4644 | 23.7162 |
| 相摩尔分数 | 0.863523 | 0.136477 |
| 相体积分数 | 0.861391 | 0.138609 |
| 摩尔体积 / m³ mol⁻¹ | 6.92577 × 10⁻⁶ | 7.05137 × 10⁻⁶ |
| 平衡计算中的相名 | FCC_L12#1，DISORD | FCC_L12#2，ORD |

γ′约占 **13.86% 体积**。这是无曲率、无相干弹性约束的体相平衡结果，不应作为任意初始颗粒分布的固定体积分数。

## 文件索引

| 文件 | 内容 |
| --- | --- |
| `equilibrium_target_1000K_xAl014.csv` | 上表的机器可读结果 |
| `gamma_1000K.csv` | 无序 γ 单相，5–30 at.% Al，502 点 |
| `gamma_prime_1000K.csv` | 有序 γ′单相，18–30 at.% Al，242 点 |
| `equilibrium_composition_scan_1000K.csv` | 完整相平衡，5–30 at.% Al，254 点 |
| `thermodynamics_overview.png`、`.pdf` | 自由能、扩散势和目标平衡体积分数图 |
| `metadata.json` | 软件、数据库、方法、单位、输入 SHA256 和适用范围 |
| `validation.json` | 数值检查结果 |

单相扫描最大成分步长为 `Δx_Al=0.0005`，即 **0.05 at.%**。两条网格由各自平衡点向两侧扫描，网格点不完全相同。
全部原始文件在 `../raw/`，计算命令在 `../macros/`，整理脚本为 `../process_exports.py`。

## 两条自由能分支如何取得

无序 γ：仅从 NIDEMO 读取 `DIS_FCC_A1`，计算均匀无序面心立方固溶体。该分支在 γ 平衡成分处与完整平衡中的 `FCC_L12#1 (DISORD)` 一致。

有序 γ′：从目标两相平衡恢复 `FCC_L12#2` 的有序状态，挂起其他相及组元集，关闭全局极小化，沿成分方向连续求单相局部平衡。子晶格占位仍参与极小化，并非固定化学计量 Ni₃Al。全范围中相摩尔分数保持 1，且 Al 两子晶格占位差大于 0.6，确认未跳回无序分支。此占位差只是本数据集中的有序性检查指标，不直接等同于 MOOSE 相场序参量。

单相表包含对其他相而言可能亚稳的分支。完整平衡表是各相混合后的最低 Gibbs 自由能，**不能直接替代 KKS 模型中的单相自由能**。扫描高 Al 端还出现 BCC_B2#2；完整平衡 CSV 中保留了相区名称及其他相总分数。缺席相的成分、占位和体积记为空值，避免误读求解器留存的起始值。

## 单相 CSV 列与单位

| 列 | 含义 |
| --- | --- |
| `x_Al` | Al 原子/摩尔分数，0–1 |
| `Gm_raw_J_mol` | Thermo-Calc 直接输出的摩尔 Gibbs 自由能 |
| `mu_Al_raw_J_mol`、`mu_Ni_raw_J_mol` | 直接输出的化学势，参考态为默认 SER |
| `Vm_m3_mol` | 该成分下的摩尔体积 |
| `dGm_dx_raw_J_mol` | 直接输出的 `μ_Al−μ_Ni` |
| `Gm_reduced_J_mol` | 扣除两相共用线性参考后的自由能 |
| `dGm_reduced_dx_J_mol` | 扣除共用参考斜率后的扩散势 |
| `Gm_J_mol`、`dGm_dx_J_mol` | 由平移量还原的自由能、导数 |
| `mu_Al_J_mol`、`mu_Ni_J_mol` | 由还原后的 G 和导数按 Euler 关系计算的化学势 |
| `phase_region` | 原始扫描相区标签 |

γ′表额外给出相摩尔分数、Al 在两个置换子晶格上的占位，以及 `y₂−y₁`。

设 `c=x_Al`，采用两相共用的参考直线，单位 J/mol：

```text
L(c) = -47425.823 - 135156.45*c
G_reduced(c) = G(c) - L(c)
dG_reduced/dc = μ_Al - μ_Ni + 135156.45
G(c) = G_reduced(c) + L(c)
μ_Ni = G - c*(dG/dc)
μ_Al = G + (1-c)*(dG/dc)
```

这条参考直线近似目标平衡公切线。对两相减去同一线性函数不会改变公切线成分。
Thermo-Calc 的 TABULATE 通常输出六位有效数字；先在软件内部计算平移后的量，再导出，可保留小自由能差。
还原列中更多的小数位代表数值还原，不代表数据库具有相应物理精度。原始列保留以便核查。

## 校验结果

- 两条单相曲线独立求公切线得到 `cγ≈0.12464387`、`cγ′≈0.23716193`；与完整平衡导出的成分最大偏差约 `1.34×10⁻⁷`。
- 目标总体 Al 分数的杠杆定则误差约 `1.19×10⁻⁷`，与原始输出舍入精度一致。
- 目标摩尔体积的相分数加权误差约 `1.51×10⁻¹² m³/mol`。
- 原始单相表的 `G=(1−c)μ_Ni+cμ_Al` 最大残差约 `0.203 J/mol`，与原始六位有效数字舍入一致。
- 单相表无非有限数、无重复成分点；γ′扫描全程保持有序单相。

公切线校验使用以 G 和 μ 差为节点数据的 Cubic Hermite 插值，仅在各自导出范围内求值。MOOSE 演示另采用平衡成分附近的局部二次自由能，见 [模型说明](../../examples/ni_al_coarsening/MODEL.md)。

## 接入 MOOSE 前的约定

应使用单相表建立 KKS 等模型的两条自由能，并检查拟合/插值能否保持导数和公切线。禁止超出已导出的成分区间直接外推，尤其 γ′表不覆盖 18 at.% 以下。

若采用常摩尔体积近似，应选择两相共用、明确记录的 `Vm_ref`，再用 `f=Gm/Vm_ref` 将 J/mol 转为 J/m³。不要将两相各自的不同 `Vm(c)` 分别代入简单 KKS 自由能而仍沿用同一摩尔分数约束；可变摩尔体积需要相容的守恒和化学势定义。CSV 保留实际 `Vm(c)`；当前 MOOSE 演示使用共同参考体积 `7.0×10⁻⁶ m³/mol`。

这些 CSV 的局部拟合已接入 [MOOSE 主输入](../../examples/ni_al_coarsening/coarsening.i)。定量熟化时间仍需迁移率/扩散系数、界面参数以及可能的弹性能；本次未从 MNIDEMO 计算动力学数据。

## 重算与整理

本次通过已授权的 Thermo-Calc Console Mode 执行 `.TCM` 宏；TC-Python 初始化报许可错误，未用该接口生成数据。

宏文件的依赖顺序为：`probe_equilibrium.TCM` → `export_equilibrium_scan.TCM` → `export_gamma.TCM` → `prepare_gamma_prime.TCM` → `export_gamma_prime.TCM` → `export_gamma_shifted_derivative.TCM`。
在 Console Mode 中使用 `macro_file_open D:/MOOSEcalc/thermodynamics/macros/<文件名>.TCM`。第一个宏从 SYS 模块启动；其他宏在前一个返回交互模式后运行。宏运行时会生成供后续宏读取的 POLY3 工作区；仓库仅保存宏和导出的数值文本，不提供这些二进制工作区。
宏内包含本机绝对路径；重算前应复制宏并将输出路径及后续 READ 路径统一指向新的目录，保留本次原始结果。最后一个宏会恢复完整体系的目标平衡及全局极小化设置。

数据整理源码保存在 `thermodynamics/process_exports.py`；读取仓库内已整理的自由能表不需要再次启动 Thermo-Calc。

软件与数据库来源：[Thermo-Calc 免费教育版及随附数据库说明](https://thermocalc.com/academia/free-educational-package/)。
命令含义另参考本机 2026b 官方 HTML5 帮助中的 `CHANGE_STATUS`、`GLOBAL_MINIMIZATION`、`STEP_WITH_OPTIONS`、`ENTER_SYMBOL` 和 `TABULATE`。
