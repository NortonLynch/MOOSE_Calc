"""Build a self-contained Chinese report from validated numerical results."""
from pathlib import Path
import json,base64,html
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/ni_al_elastic'
m=json.loads((OUT/'comparison_metrics.json').read_text())
assert m['status']=='COMPLETE','Do not present an incomplete simulation as final.'
cal=m['time_calibration'];cv=json.loads((OUT/'control_validation.json').read_text())
assert cv['status']=='PASS'
manifest=json.loads((OUT/'physical_exodus_manifest.json').read_text())
names={'no_elastic':'无弹性','elastic':'有弹性'}
def picture(name,caption):
    data=base64.b64encode((OUT/name).read_bytes()).decode()
    return f'<figure><img src="data:image/png;base64,{data}" alt="{html.escape(caption)}"><figcaption>{caption}</figcaption></figure>'
def table(head,rows):
    return '<div class="scroll"><table><thead><tr>'+''.join(f'<th>{v}</th>' for v in head)+'</tr></thead><tbody>'+''.join('<tr>'+''.join(f'<td>{v}</td>' for v in row)+'</tr>' for row in rows)+'</tbody></table></div>'
summary=[]
for k in names:
    z=m[k];lost=z['dissolution_brackets_star']['1']
    bracket='本时间段未消失' if lost is None else '–'.join(f'{v*cal["time_scale_s_per_model_unit"]/3600:.4f}' for v in lost)+' h'
    summary.append([names[k],z['final_count'],f'{z["final_mean_R_nm"]:.5f}',f'{z["final_mean_R_cubed_nm3"]:.3f}',f'{z["radii_final_nm"][3]:.5f}',bracket])
fits=[]
for k in names:
    for f in m['fits'][k]:
        fits.append([names[k],f'{f["lo_star"]*cal["time_scale_s_per_model_unit"]/3600:.3f}–{f["hi_star"]*cal["time_scale_s_per_model_unit"]/3600:.3f}',f'{f["K_nm3_s"]:+.6g}',f'{f["K_m3_s"]:+.4e}',f'{f["R_squared_time_weighted"]:.4f}',','.join(map(str,f['particle_counts'])),'跨越颗粒消失' if f['crosses_particle_loss'] else '颗粒数固定',html.escape(f['window'])])
whole={k:next(f for f in m['fits'][k] if f['window']=='5000-60000') for k in names}
fit_sentence='；'.join(f'{names[k]} K<sub>eff</sub>={whole[k]["K_nm3_s"]:.5g} nm³/s，R²={whole[k]["R_squared_time_weighted"]:.3f}' for k in names)
alt_ratio=cal['alternative_zhang_2010_time_scale_s']/cal['time_scale_s_per_model_unit']
big_change=(m['elastic']['radii_final_nm'][3]/m['no_elastic']['radii_final_nm'][3]-1)*100
loss_text=''
le=m['elastic']['dissolution_brackets_star']['1'];ln=m['no_elastic']['dissolution_brackets_star']['1']
if le and ln:
    change=((sum(le)/2)/(sum(ln)/2)-1)*100
    loss_text=f'P1 的消失时间相对无弹性改变约 {change:+.1f}%（取保存帧区间中点，仅作事件时间估计）。'
second_loss=m['elastic']['dissolution_brackets_star']['2']
if second_loss:
    loss_text+='弹性组 P2 在 '+'–'.join(f'{v*cal["time_scale_s_per_model_unit"]/3600:.4f}' for v in second_loss)+' h 内消失。'
page=f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Ni–Al 弹性耦合熟化与物理时间分析</title><style>
body{{font:16px/1.8 "Microsoft YaHei","Noto Sans CJK SC",sans-serif;color:#263443;background:#edf1f5;margin:0}}main{{max-width:1120px;margin:32px auto;padding:42px 50px;background:white}}h1{{font-size:30px;line-height:1.4}}h2{{margin-top:2em;border-bottom:2px solid #c7dce7;padding-bottom:8px;font-size:23px}}h3{{font-size:18px}}a{{color:#176c9a}}.tag{{color:#477286;font-size:13px;letter-spacing:2px}}.note{{border-left:5px solid #cf5b36;background:#fff4ed;padding:14px 22px}}.facts{{display:flex;gap:18px;flex-wrap:wrap}}.fact{{flex:1;min-width:180px;background:#edf5f9;padding:16px}}.fact b{{font-size:25px;display:block;color:#186488}}img{{width:100%;height:auto}}figure{{margin:24px 0}}figcaption{{color:#596877;font-size:14px}}table{{border-collapse:collapse;width:100%;font-size:14px;margin:16px 0}}td,th{{padding:9px 12px;border-bottom:1px solid #dae2e8;text-align:left}}th{{background:#edf5f9}}.scroll{{overflow:auto}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;background:#f1f4f7;padding:20px;font:14px/1.6 Consolas,monospace}}code{{font-family:Consolas,monospace}}.small{{font-size:14px;color:#596877}}@media(max-width:700px){{main{{margin:0;padding:22px}}h1{{font-size:25px}}}}@media print{{main{{margin:0;padding:0}}figure,table{{break-inside:avoid}}}}
</style><main>
<div class="tag">MOOSE × THERMO-CALC · 实际计算结果</div>
<h1>Ni–14 at.% Al、1000 K：<br>弹性耦合熟化、物理时间与 ⟨r⟩³ 增长检验</h1>
<p>在已完成的查表自由能算例上加入用户提供的 γ/γ′立方弹性与共格失配。两组保持相同成分、温度、初始四颗粒、网格、界面参数和动力学无量纲参数，均计算到 t*=60000。</p>
<div class="facts"><div class="fact"><b>{cal['end_time_hours']:.3f} h</b>本次终点的文献模型时间</div><div class="fact"><b>{cal['gamma_chemical_interdiffusivity_m2_s']:.3e} m²/s</b>γ 基体化学互扩散系数</div><div class="fact"><b>{cal['time_scale_s_per_model_unit']:.6f} s</b>每个模型时间单位</div></div>
<div class="note"><b>计算完成；本算例不能确认普适的 LSW 线性增长。</b> 已输出形貌、应力、⟨r⟩³–物理时间曲线和有效拟合斜率。四颗粒有限体系中，小颗粒消失会改变数平均的分母，造成跳变；不能把跨越此事件的回归斜率当作材料粗化常数。时间轴来自文献迁移率模型与当前热力学表的组合，尚未以本合金的 1000 K 实测时效数据标定。</div>
<h2>1. 同一算例的形貌对照</h2>
{picture('morphology_comparison.png','蓝色为 γ′，黑线为 η=0.5。各帧标注实际保存时间；坐标为 nm。平面应变模型中两个面内平均正应力及剪应力为零，εzz=0。')}
{table(['算例','终态颗粒数','终态 ⟨r⟩ / nm','终态 ⟨r⟩³ / nm³','最大初始颗粒 P4 终态半径 / nm','最小初始颗粒 P1 消失区间'],summary)}
<p>弹性项改变颗粒形状与化学—界面—弹性能之间的平衡。立方弹性允许出现四重对称轮廓；共格失配提高保留析出相的能量代价，因此初期相分数和平均半径也会调整。判断熟化需同时观察每个颗粒：小颗粒缩小、较大颗粒长大，与所有颗粒一致增大或直接接触合并是不同过程。实际时间历史见下图与 CSV。</p>
<p>在相同终点，弹性模型的 P4 半径相对无弹性改变 {big_change:+.2f}%；{loss_text}终态 γ′面积分数分别为 {m['no_elastic']['final_area_fraction']*100:.4f}%（无弹性）与 {m['elastic']['final_area_fraction']*100:.4f}%（有弹性）。这些是本模型的实际对照结果，不直接推广为所有镍基合金的规律。</p>
<p>P4 的 η=0.5 轮廓上，四个对角方向半径的平均值与 x/y 轴方向半径平均值之比分别为 {m['no_elastic']['P4_final_shape']['diagonal_axis_ratio']:.4f}（无弹性）和 {m['elastic']['P4_final_shape']['diagonal_axis_ratio']:.4f}（有弹性）；圆形对应 1。该指标围绕相分数加权形心测量，用于定量补充形貌图。</p>
<p>终态基体核心平均 Al 成分分别为 {m['no_elastic']['c_gamma_final_core_mean']*100:.5f} at.% 与 {m['elastic']['c_gamma_final_core_mean']*100:.5f} at.%；γ′核心分别为 {m['no_elastic']['c_gamma_prime_final_core_mean']*100:.5f} at.% 与 {m['elastic']['c_gamma_prime_final_core_mean']*100:.5f} at.%。这里统计 η&lt;0.01 或 η&gt;0.99 的单元核心，属于有限颗粒、有限时间的成分，不能当作新的平直界面体相平衡成分。</p>
{picture('particle_radii_comparison.png','P1–P4 对应初始名义半径 9.2、9.7、10.3、10.8 nm；扩散界面导致初始积分半径略有差异。曲线以相分数积分得到等面积半径，颗粒溶解后的微小残余相场积分不计入颗粒数量平均。')}
{picture('elastic_stress_final.png','终态面内应力，单位 MPa。零平均应力不表示局部应力为零；平面应变约束仍允许非零 σzz。')}
<h2 id="kinetics">2. ⟨r⟩³–时间与速率提取</h2>
<pre>rᵢ = √[ (∫Ωᵢ h(η) dA) / π ]
h(η) = η³(6η² − 15η + 10)
⟨r⟩ = (1/N) Σᵢ rᵢ          （仅统计仍存在的颗粒）
y = ⟨r⟩³                  （不是 ⟨r³⟩）
y ≈ A + K_eff t            t 以秒计；K_eff 单位 nm³/s
K [m³/s] = K [nm³/s] × 10⁻²⁷</pre>
<p>每个初始颗粒占一个象限，积分使用节点场与梯形权重。以 max(η)≥0.5 判断颗粒存在，并检查连通区域数与四象限分配一致、颗粒未穿越边界。无额外平滑；拟合采用时间加权最小二乘，避免自适应小步在溶解事件附近产生过多输出点而支配回归。截距自由拟合。</p>
{picture('mean_radius_cubed_vs_time.png','上图是数平均半径的立方；下图同步显示颗粒数。虚线仅为 t*=5000–60000 的有效回归，不能跨越颗粒消失事件解释为稳态粗化定律。')}
<p><b>全窗口结果：</b>{fit_sentence}。下表同时保留预先指定窗口和所有具有足够样本的固定颗粒数区间，避免只挑选看似线性的部分。负斜率表示该窗口的数平均半径立方下降；它不是“负的 LSW 粗化常数”。即使局部 R² 很高，也需检查斜率符号、颗粒数和自相似分布。</p>
{table(['算例','物理时间 / h','K_eff / nm³·s⁻¹','K_eff / m³·s⁻¹','时间加权 R²','N','事件','模型时间窗口'],fits)}
<p><b>结论边界：</b>本次完成了线性增长检验与有效速率提取，但没有验证三维统计意义上的 LSW 常数。四颗粒、二维等面积半径、约 13% 相分数、有限盒子与共格弹性均限制这种解释。后续若要检验稳态标度，应在保持本次结果作为对照的基础上增加独立颗粒数与空间尺寸，优先做三维统计体系，并检查多段拟合、粒径分布自相似性和网格/时间步收敛。</p>
<p>四颗粒体系的一个直接解释是面积恒等式：<b>⟨r⟩² = (ΣAᵢ)/(πN) − Var(r)</b>。当 N 固定、总析出面积近似不变时，大颗粒长大、小颗粒缩小使粒径方差增加，数平均半径反而可能减小。N 由 4 降到 3、再由 3 降到 2 时则发生离散跳升。这个有限数量效应说明，观察到熟化并不保证本次 ⟨r⟩³ 曲线已经进入稳定线性区间；更不能改用 ⟨r³⟩ 来替代用户要求的统计量。</p>
<h2 id="calibration">3. 文献动力学与物理时间换算</h2>
<p>主时间标尺采用 Engström–Ågren 的 fcc 迁移率评估（1996），参数依据 <a href="https://doi.org/10.3139/146.110428">Zhang 等（2010）表 3 的既有评估行 [16,22]</a>核对；热力学因子取本项目 NIDEMO 的 γ 单相查表曲率。这是本文计算的有效化学扩散系数，不能把它写成 Janssen 在 1000 K 的直接测量值。</p>
<pre>c = x_Al^γ = {cal['calibration_gamma_Al_atomic_fraction']:.10f},  T = 1000 K
U_Al = c(−142000−72.12T) + (1−c)(−284000−59.83T)
       + c(1−c)(−41300−91.2T)
U_Ni = c(−145900−64.26T) + (1−c)(−287000−69.8T)
       + c(1−c)(−113000+65.5T)              [J/mol atoms]
D_i^tracer = exp(U_i / RT)
Φ = c(1−c) Gγ″ / RT
D̃γ = [(1−c) D_Al^tracer + c D_Ni^tracer] Φ</pre>
{table(['量','计算值'],[['D_Al^tracer',f'{cal["Al_tracer_diffusivity_m2_s"]:.6e} m²/s'],['D_Ni^tracer',f'{cal["Ni_tracer_diffusivity_m2_s"]:.6e} m²/s'],['Gγ″',f'{cal["gamma_Gm_second_derivative_J_mol"]:.6f} J/mol'],['热力学因子 Φ',f'{cal["thermodynamic_factor"]:.6f}'],['D̃γ',f'{cal["gamma_chemical_interdiffusivity_m2_s"]:.6e} m²/s']])}
<pre>l₀ = 1 nm, G₀ = 1000 J/mol, M* = {cal['dimensionless_mobility']:.16g}
D̃γ = M* (Gγ″/G₀) l₀² / t₀
t₀ = {cal['time_scale_s_per_model_unit']:.12f} s
t_physical = t* × t₀
60000 model units = {cal['end_time_s']:.3f} s = {cal['end_time_hours']:.6f} h</pre>
<p>没有把无量纲时间直接当秒。此处校准在 γ 基体平衡成分约 12.464 at.% Al，而非两相混合的总成分 14 at.% Al。保持原 M*、L* 不变，整个时间轴统一缩放；没有同时引入空间成分相关迁移率或应力相关扩散系数。高温评估向 1000 K 的使用及与当前热力学数据库组合都带来模型误差。</p>
<p>原模型的同一个常数迁移率贯穿两相；本次只将时间标尺对应到 γ 基体扩散，没有另外标定有序 γ′的迁移率或界面动力学。界面能 20 mJ/m²与有限界面宽度仍沿用原参考设置。因此带秒/小时的时间轴具有明确的模型依据，但整个熟化过程仍不能视为已完成实验验证的定量预测。</p>
<p>作为参数敏感性对照，采用 Zhang 2010 表 3 本文评估参数可得 D̃γ≈{cal['alternative_zhang_2010_Dchem_m2_s']:.3e} m²/s、t₀≈{cal['alternative_zhang_2010_time_scale_s']:.6f} s，时间轴变长约 {(alt_ratio-1)*100:.1f}%，所有物理速率相应乘 {1/alt_ratio:.3f}。这只是两种评估的差异，不是统计置信区间。</p>
<h2 id="model">4. 弹性数据、单位与求解配置</h2>
{table(['参数','γ','γ′'],[['C11 / GPa','199','207.8'],['C12 / GPa','143','148.8'],['C44 / GPa','99','102.3'],['失配本征应变','0','+0.0047 I']])}
<pre>C(η) = (1−h) Cγ + h Cγ′
ε⁰(η) = 0.0047 h I
f_el = ½(ε−ε⁰) : C(η) : (ε−ε⁰)
f₀ = G₀/Vm = 1.4285714286×10⁸ J/m³,  Vm = 7×10⁻⁶ m³/mol
C* = C/f₀；若 C 输入以 GPa 计，数值乘 7。</pre>
<p>附件的 GPa 数值按原样采用；其 eV/nm³ 注释存在 1000 倍单位差异：199 GPa 对应 1242.0603 eV/nm³，不是 1.242。计算以 GPa 为准。共用晶向 [100]∥x、[010]∥y、[001]∥z。二维小应变、εzz=0，位移波动与相场均周期；单点去除刚体平移，三个全局应变自由度使面内体积平均应力为零。</p>
<p>先固定原始 c/η/μ 求静态力学平衡，再将相同节点场与平衡位移读入瞬态；这一静态计算不消耗扩散时效时间。ElasticEnergyMaterial 提供 ∂f_el/∂η，加入 Allen–Cahn 方程，并同步求应力平衡。KKS 两相化学自由能表保持原样；弹性能没有显式 c 依赖，因此固定 η 时的 KKS 化学分配关系未改动。</p>
<p>MOOSE 官方包 2026.10.14（提交 c140bfd5a6），phase_field + solid_mechanics；Thermo-Calc Educational 2026b，NIDEMO 2.0。96×96 nm，96×96 Q1 单元；隐式 Euler；PJFNK、全耦合 SMP 预条件与 MUMPS，自动方程缩放每步更新，8 个 MPI 进程。时间序列继承原算例；发生不收敛时允许步长细分。原无弹性结果为上一轮查表算例，而非二次近似。</p>
<p>两组在 t*≤30000 的实际接受时间网格是否一致：<b>{m['integration_time_grid_comparison']['identical_through_30000']}</b>；全程是否一致：<b>{m['integration_time_grid_comparison']['identical_accepted_grid']}</b>。若颗粒消失附近发生额外步长细分，事件区间和晚期差异仍含时间离散误差。本轮未进行完整网格与步长加密研究，不能把守恒、能量下降及零失配检查等同于收敛性证明。</p>
<h2 id="validation">5. 校验与输出说明</h2>
{picture('conservation_energy_fraction.png','质量、总自由能及 γ′面积分数历史。二维积分能量换算为每单位面外厚度的 J/m。')}
{table(['校验','结果'],[['零失配对照最大 η 差异',f'{cv["zero_misfit_max_nodal_difference"]["eta"]:.3e}'],['零失配对照最大 c 差异',f'{cv["zero_misfit_max_nodal_difference"]["c"]:.3e}'],['稳态到瞬态 c/η/μ/位移初场差异',f'{max(cv["steady_to_transient_initial_nodal_difference"].values()):.3e}'],['无弹性 / 有弹性最大总 Al 误差',f'{m["no_elastic"]["max_mass_error"]:.3e} / {m["elastic"]["max_mass_error"]:.3e}'],['有弹性最大平均面内应力残差',f'{m["elastic"]["max_average_inplane_stress_MPa"]:.3e} MPa'],['无弹性 / 有弹性最大能量增加',f'{m["no_elastic"]["maximum_energy_increase_scaled"]:.3e} / {m["elastic"]["maximum_energy_increase_scaled"]:.3e}（模型单位）']])}
<p class="small">初始输出处理：MOOSE 的 ComputeEigenstrainBase 在 step=0 跳过本征应变，导致瞬态原始 INITIAL 应力/弹性能辅助输出与已导入的平衡位移不一致。报告及物理时间 Exodus 的第 0 帧仅将这些机械辅助量替换为独立稳态解的末帧数值；c、η、μ、位移节点场完全不改，已经逐点验证一致。原始 elastic.e、CSV 保留以便追溯。之后所有帧均来自实际瞬态求解。</p>
<p>ParaView 打开 <a href="elastic_seconds.e">elastic_seconds.e</a> 与 <a href="no_elastic_seconds.e">no_elastic_seconds.e</a>，勾选变量后 Apply，按 η 着色并播放时间。文件时间单位为秒，坐标及位移波动单位为 nm，stress_*_MPa 为 MPa；energy_density 为 J/m³，μ 仍以 1000 J/mol 为单位。全局应变引起的均匀变形没有加到位移波动变量中，比较形貌时使用参考坐标即可。</p>
<p>下载/再分析：<a href="coarsening_rate_fits.csv">速率拟合表</a> · <a href="elastic_particle_history.csv">弹性颗粒历史</a> · <a href="no_elastic_particle_history.csv">无弹性颗粒历史</a> · <a href="time_calibration.json">时间校准公式与数值</a> · <a href="physical_exodus_manifest.json">场文件校验与单位</a> · <a href="elastic_comparison.gif">演化动画</a>。</p>
<h2>6. 重算与来源</h2>
<p>在 WSL 项目根目录运行 <code>bash scripts/run_elastic_comparison.sh</code>。该脚本要求主结果目录为空，避免覆盖已完成场文件。仅重新作图、换算或更新报告：</p>
<pre>source ~/miniforge3-moose/etc/profile.d/conda.sh
conda activate moose-post
python scripts/calibrate_elastic_time.py
python scripts/check_elastic_controls.py
python scripts/analyze_elastic_comparison.py
python scripts/export_elastic_physical.py
python scripts/write_elastic_report.py</pre>
<ul>
<li>Engström &amp; Ågren (1996), <a href="https://doi.org/10.1515/ijmr-1996-870205">Assessment of diffusional mobilities in face-centered cubic Ni-Cr-Al alloys</a>，主动力学评估。</li>
<li>Zhang, Du, Chen, Steinbach &amp; Huang (2010), <a href="https://doi.org/10.3139/146.110428">Atomic mobilities and diffusivities in the fcc, L1₂ and B2 phases of the Ni–Al system</a>，表 3 用于参数核对与敏感性比较。</li>
<li>Campbell, Boettinger &amp; Kattner (2002), <a href="https://doi.org/10.1016/S1359-6454(01)00383-4">Development of a diffusion mobility database for Ni-base superalloys</a>，用户提供的数据库路线文献；本次没有直接运行商业 MOBNI 数据库。</li>
<li>Zhu et al. (2004), <a href="https://doi.org/10.1016/j.actamat.2004.02.032">Three-dimensional phase-field simulations of coarsening kinetics of γ′ particles in binary Ni–Al alloys</a>，用户提供的三维定量相场参考；本次是原二维四颗粒算例的弹性扩展，不等同于复现该论文。</li>
<li>Janssen (1973), <a href="https://doi.org/10.1007/BF02668017">Diffusion in the nickel-rich part of the Ni–Al system at 1000° to 1300°C</a>，温度为摄氏度范围，未当作 1000 K 实测数值使用。</li>
<li><a href="https://mooseframework.inl.gov/modules/phase_field/Mechanics_Coupling.html">MOOSE Mechanics Coupling</a>；<a href="https://mooseframework.inl.gov/source/materials/ElasticEnergyMaterial.html">ElasticEnergyMaterial</a>；<a href="https://github.com/idaholab/moose/blob/master/modules/solid_mechanics/src/materials/ComputeEigenstrainBase.C">初始本征应变处理源码</a>。</li></ul>
<p class="small">本报告的曲线、拟合与图像全部由本地 MOOSE 场数据计算生成。脚本、原始输入、附件副本与校准参数均保存在项目中。</p></main></html>'''
(OUT/'REPORT.html').write_text(page,encoding='utf-8')
print(OUT/'REPORT.html')
