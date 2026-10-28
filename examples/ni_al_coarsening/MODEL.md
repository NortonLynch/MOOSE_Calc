# Ni–Al γ′熟化演示模型

这是采用参考参数建立的二维、等温、扩散耦合相场模型。软件为 MOOSE 2026.10.14；热力学输入来自 Thermo-Calc Educational 2026b / NIDEMO 2.0 的 1000 K 导出结果。

**本模型能够演示奥斯特瓦尔德熟化，但尚不预测经过实验标定的真实时效时间。** 没有引入颗粒运动、吸引力或人为删除颗粒。颗粒通过守恒浓度场交换 Al，界面由自由能驱动演化。

## 热力学近似与方程

对两相减去同一公切线，在各自平衡成分附近使用二次近似：

\[
G_\gamma(c_\gamma)=\tfrac12 A_\gamma(c_\gamma-c_\gamma^0)^2,\qquad
G_{\gamma'}(c_{\gamma'})=\tfrac12 A_{\gamma'}(c_{\gamma'}-c_{\gamma'}^0)^2.
\]

平衡成分 `cγ⁰=0.124644`、`cγ′⁰=0.237162` 来自完整 Thermo-Calc 平衡。曲率由平衡点附近 `±0.00101` 成分窗口中的化学势差斜率拟合：

| 参数 | 数值 |
| --- | ---: |
| Aγ | 347973 J/mol |
| Aγ′ | 424347.4 J/mol |

拟合和原始曲线对比见结果目录中的 `free_energy_approximation.png`。这不是对整个成分区间的 CALPHAD 精确拟合。

使用五次切换函数 `h(η)=η³(6η²−15η+10)`，η=0 为 γ，η=1 为 γ′。对二次相自由能，KKS 的相浓度和等化学势约束可解析消元：

\[
q(h)=\frac{1-h}{A_\gamma}+\frac{h}{A_{\gamma'}},\quad
\mu=\frac{c-(1-h)c_\gamma^0-hc_{\gamma'}^0}{q(h)},\quad
G_\text{chem}=\frac{[c-(1-h)c_\gamma^0-hc_{\gamma'}^0]^2}{2q(h)}.
\]

相浓度恢复为 `cγ=cγ⁰+μ/Aγ`、`cγ′=cγ′⁰+μ/Aγ′`，同时满足质量约束和等化学势。这使求解变量由显式 KKS 的五个减为 `c、η、μ` 三个。

体自由能使用两相共同的常摩尔体积。总泛函为化学能、`Wη²(1−η)²` 和 `κη|∇η|²/2` 的积分，求解：

\[
\partial_t c=\nabla\cdot(M\nabla\mu),\qquad
\partial_t\eta=-L\,\delta F/\delta\eta.
\]

具体使用 MOOSE 的 [SplitCHParsed](https://mooseframework.inl.gov/source/kernels/SplitCHParsed.html)、SplitCHWRes、CoupledTimeDerivative、[AllenCahn](https://mooseframework.inl.gov/source/kernels/AllenCahn.html) 和 ACInterface。自由能导数由 DerivativeParsedMaterial 自动构造。消元形式与二次自由能下的 KKS 局部约束等价，不是把两相平衡包络作为单相势能。

## 参数及来源

| 项目 | 本次设置 | 性质 |
| --- | --- | --- |
| 温度、总体 Al 分数 | 1000 K、0.14 | 用户目标 |
| 两相平衡成分及局部曲率 | 见上表 | Thermo-Calc 导出及局部拟合 |
| 计算域 | 96×96 nm²，二维周期边界 | 演示设定 |
| 初始半径 | 9.2、9.7、10.3、10.8 nm | 约 ±8% 的初始尺寸差异 |
| 初始圆心 / nm | (24,24)、(72,24)、(24,72)、(72,72) | 演示设定 |
| 共同摩尔体积 | 7.0×10⁻⁶ m³/mol | 参考近似，接近导出的两相体积 |
| 界面能 | 0.020 J/m² | 参考假设，未经本合金实验标定 |
| 界面宽度 | η 从 0.1 到 0.9 的距离为 4 nm | 数值正则化 |
| γ 化学扩散系数 | 1×10⁻¹⁸ m²/s | 参考假设，未从 MNIDEMO 提取 |
| 溶质迁移率 | 常数，按 γ 相曲率与上述 D 换算 | 正的守恒扩散系数 |
| 无量纲界面迁移率 L | 10 | 数值参考设定，未作界面动力学标定 |

网格为 96×96 个 Q1 四边形，间距 1 nm；隐式 Euler 时间推进，Newton + MUMPS，4 MPI 进程。初始 η 由四个平滑圆构造，并用统一初始扩散势调整浓度，使离散积分的总体 Al 分数恰为 0.14。

对于 `Wη²(1−η)²+κη|∇η|²/2` 的平直界面，界面能为 `sqrt(2κηW)/6`；0.1–0.9 宽度为 `ln(9)*sqrt(2κη/W)`。据此一致地确定 W 和 κ。参见 [MOOSE KKS 平衡界面说明](https://mooseframework.inl.gov/modules/phase_field/MultiPhase/KKSAnalytical.html)。

长度标度 `l₀=1 nm`，摩尔能标度 `G₀=1000 J/mol`，体能标度 `f₀=G₀/Vm`，时间标度 `t₀=l₀²/Dref=1 s`。Exodus 中坐标数值以 nm 解释，时间值为 `t*`。虽然按这组参考参数有 `t*=t/s`，**不能把它当作已验证的真实时效秒数**；更换扩散参数会改变物理时间映射。

## 本次结果如何判断

- 检查总体 Al 分数是否守恒，而不是只看图像变化。
- 检查总自由能下降。
- 分别跟踪四个初始位置周围 `h(η)` 的面积，以 `R=sqrt(area/π)` 定义等效半径。
- 以 `η≥0.5` 的连通区域计数，识别小颗粒消失。
- 计算各象限 Al 总量变化和扩散势梯度，核实溶质由小颗粒区域转移到大颗粒区域。

有限元节点量输出后按 Q1 面积权重积分。粒子消失时间是首次保存的无 `η≥0.5` 区域时刻，精度受输出间隔限制。四个颗粒保持分离，故象限跟踪有意义；该跟踪法不适用于颗粒跨区移动或并合后的普遍情况。

## 适用边界

这是二维圆柱截面模型，不是三维球形颗粒的定量熟化速率。未考虑相干弹性、晶格失配、晶向、γ′反相畴或各向异性界面能。界面宽度是数值宽度；尚未完成系统的网格、时间步、界面宽度和界面迁移率收敛研究。

获得可信的定量时效时间和粒径分布，还需要适用于本体系的迁移率/扩散数据、界面参数及弹性数据，并开展上述收敛和实验验证。本次不能用任意时间点的图像与实验时效小时数直接对应。

## 运行与复现

主输入文件为同目录的 `coarsening.i`，参数记录为 `results/ni_al_coarsening/model_parameters.json`。

```powershell
# 新运行：指定新文件基名，保留现有结果
wsl.exe -d Ubuntu-24.04 --cd (Get-Location).Path -- bash scripts/run_ni_al_coarsening.sh 4 Outputs/file_base=results/ni_al_coarsening/rerun
```

初始算例通过 `scripts/prepare_ni_al_coarsening.py` 生成。该脚本会重写主输入文件和模型参数；在自定义参数前请保存自己的副本。本次计算曾从检查点恢复以调整求解步长控制；未改变材料参数或手动更改浓度/颗粒字段。相关日志、完整 Exodus 和检查点仅保留在原计算机器，仓库发布关键结果图表及数值校验记录。干净重算不保证接受时间步序列完全相同。
