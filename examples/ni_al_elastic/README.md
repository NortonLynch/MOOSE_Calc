# Ni–Al 查表自由能 + 共格弹性

基于上一轮 `examples/ni_al_lookup/lookup.i`，保留 Ni–14 at.% Al、1000 K、
γ/γ′、96×96 nm 周期域、96×96 Q1 网格、四颗粒初态及界面/迁移率参数。
新增相依立方弹性与 `0.0047 h(eta) I` 本征应变。

| 参数（GPa） | γ | γ′ |
| --- | ---: | ---: |
| C11 | 199 | 207.8 |
| C12 | 143 | 148.8 |
| C44 | 99 | 102.3 |

数值来自用户附件副本 `plasti_source.txt`。附件中的 eV/nm³ 注释小了 1000 倍；
以 GPa 为准，采用 `f0=1000/(7e-6) J/m³` 后，GPa 数值乘 7 得到 MOOSE 输入刚度。
无塑性、本征应变各向同性、各颗粒共用晶向 [100]/[010]/[001]。

机械边界为二维平面应变 `epsilon_zz=0`，位移波动周期，单点消除平移。
`GlobalStrain` 的三个自由度使平均 xx、yy、xy 应力为零。
这不等同于三维自由膨胀或平面应力。

文件顺序：

1. `elastic_initial.i`：固定原 c/eta/mu，求初始力学平衡。
2. `elastic_relaxed.i`：读入平衡位移和全局应变后进行熟化。
3. `zero_misfit.i`：失配置零的退化极限校验。
4. `elastic.i`：生成器基础输入，位移默认零；正式结果使用 relaxed 输入。

`scripts/prepare_elastic_comparison.py` 生成基础输入；
`scripts/initialize_elastic_transient.py` 在稳态场存在后生成 relaxed 输入。
求解器使用 MOOSE 2026.10.14 官方预编译包，包含 phase_field 和 solid_mechanics，
PJFNK + full SMP + MUMPS + 每步自动缩放，隐式 Euler，正式运行 8 MPI 进程。
热力学来自 Thermo-Calc Educational 2026b、NIDEMO 2.0 的原自由能表。

在已恢复原查表结果的 WSL 项目根目录运行：

```bash
bash scripts/run_elastic_comparison.sh
```

首次使用先按仓库根 README 部署 MOOSE，再创建后处理环境并恢复基准场：

```bash
source ~/miniforge3-moose/etc/profile.d/conda.sh
conda create -n moose-post --file environment/postprocess-linux-64.lock.txt  # 已有环境时跳过
conda activate moose-post
python scripts/download_lookup_results.py
bash scripts/run_elastic_comparison.sh
```

若只查看或重新分析已发布的弹性结果，运行 `python scripts/download_elastic_results.py` 即可恢复所需场文件；该下载包含已完成的 `elastic.e`，因此应在另外的项目副本中执行全新求解。

脚本拒绝覆盖已存在的 `results/ni_al_elastic/elastic.e`；需要全新重算时使用独立项目副本，
保留 `results/ni_al_lookup/lookup.e` 作为无弹性基准，并使新副本的弹性结果目录为空。
原始 Thermo-Calc 导出数据可直接重用，无需重新取得软件许可证运行后处理。

物理时间由 `scripts/calibrate_elastic_time.py` 计算：Engström–Ågren fcc 迁移率
与当前 γ 自由能曲率组合，在 γ 基体平衡成分 0.1246438664 处校准。
`Dchem≈6.0247e-18 m²/s`、`t0≈0.16598177 s`、终点约 2.76636 h。
这是文献模型换算，尚未实验标定；没有把 Al 示踪扩散系数直接等同化学互扩散系数。
公式、原始参数、替代评估敏感性及来源保存在 `time_calibration.json` 和报告中。

后处理入口：

```bash
source ~/miniforge3-moose/etc/profile.d/conda.sh
conda activate moose-post
python scripts/calibrate_elastic_time.py
python scripts/check_elastic_controls.py
python scripts/analyze_elastic_comparison.py
python scripts/export_elastic_physical.py
python scripts/write_elastic_report.py
python scripts/analyze_single_particle_growth.py
python scripts/write_single_particle_report.py
```

产物在 `results/ni_al_elastic/`：`REPORT.html`、对比图/GIF、半径历史 CSV、
有效速率表和 `elastic_seconds.e` / `no_elastic_seconds.e`。
物理时间 Exodus 的时间单位为秒，坐标/位移波动为 nm，应力为 MPa；
能量密度为 J/m³，整体二维积分能量为每单位面外厚度的 J/m。
数平均半径立方与各颗粒立方的平均分别导出，主图使用前者。

追加的固定单颗粒统计位于 `results/ni_al_elastic/single_particle/`，报告为其中的
`REPORT.html`。它从原始 Exodus 重新提取终态最大颗粒 P4 的面积，固定跟踪初始
中心 `(72,72) nm` 的同一颗粒，输出 `A^(3/2)` 与 `r_eq^3=(A/pi)^(3/2)` 随秒/小时变化。
独立 CSV、PNG/SVG、时间窗口拟合和轮廓面积复核均在此目录中。只更新此统计时运行
上述最后两条 Python 命令，无需重新求解相场。

四颗粒二维体系中的数平均值会在颗粒消失时跳变，不能仅根据跨事件直线拟合
宣称验证了三维 LSW 定律。报告保留全部预定义窗口及固定颗粒数区间的拟合。

MOOSE 初始输出注意事项：`ComputeEigenstrainBase` 在第零步跳过本征应变。
正式瞬态的节点相场和位移正确读入，但原始第零帧机械辅助输出仍不含该本征应变。
后处理验证节点场与稳态相同，仅在物理时间副本/分析中将第零帧的机械辅助量
替换为对应稳态末帧；原始 `elastic.e` 和 CSV 保留。
