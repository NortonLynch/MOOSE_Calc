"""Report fixed-particle area and radius-cubed kinetics from computed field data."""
from pathlib import Path
import json, base64, html
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/ni_al_elastic/single_particle'
m = json.loads((OUT/'metrics.json').read_text())
assert m['status'] == 'COMPLETE' and m['validation'].startswith('PASS')
fits = pd.read_csv(OUT/'P4_growth_fits.csv')
names = {'no_elastic': '无弹性', 'elastic': '有弹性'}
factor = m['pi_power_3_2']


def pic(name, caption):
    data = base64.b64encode((OUT/(name+'.png')).read_bytes()).decode()
    return f'<figure><img src="data:image/png;base64,{data}" alt="{html.escape(caption)}"><figcaption>{caption} <a href="{name}.svg">矢量图 SVG</a></figcaption></figure>'


def table(head, rows):
    return '<div class="scroll"><table><thead><tr>'+''.join(f'<th>{x}</th>' for x in head)+'</tr></thead><tbody>'+''.join('<tr>'+''.join(f'<td>{x}</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table></div>'


summary = []
main_fits = []
area_check = []
for key, title in names.items():
    z = m['cases'][key]
    summary.append([title, f'{z["initial_radius_nm"]:.6f}', f'{z["final_radius_nm"]:.6f}',
                    f'{z["final_area_nm2"]:.6f}', f'{z["final_r3_nm3"]:.6f}', f'{z["final_A32_nm3"]:.6f}'])
    for field in ['primary_fit', 'whole_fit']:
        f = z[field]
        main_fits.append([title, f'{f["lo_h"]:.3f}–{f["hi_h"]:.3f}', f'{f["K_nm3_s"]:.7f}',
                          f'{f["K_m3_s"]:.5e}', f'{f["K_nm3_s"]*factor:.7f}', f'{f["R2"]:.6f}',
                          f'{f["RMSE_percent_of_observed_range"]:.2f}%', f't*={f["window"]}'])
    a, b = z['primary_fit'], z['contour_primary_fit']
    area_check.append([title, f'{a["K_nm3_s"]:.7f}', f'{b["K_nm3_s"]:.7f}',
                       f'{100*(b["K_nm3_s"]/a["K_nm3_s"]-1):+.3f}%', f'{b["R2"]:.6f}'])
all_rows = []
for _, f in fits[fits.metric == 'radius_cubed_nm3'].iterrows():
    all_rows.append([names[f['case']], html.escape(f['window']), f'{f.lo_h:.3f}–{f.hi_h:.3f}',
                     f'{f.K_nm3_s:.7f}', f'{f.R2:.6f}', f'{f.RMSE_percent_of_observed_range:.2f}%', f'{f.quadratic_R2:.6f}'])
no, el = m['cases']['no_elastic'], m['cases']['elastic']
ratio = m['primary_K_elastic_to_no_elastic_ratio']
whole_ratio = el['whole_fit']['K_nm3_s']/no['whole_fit']['K_nm3_s']

page = f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Ni–Al：固定大颗粒 P4 的 r³ 与物理时间</title>
<style>body{{font:16px/1.85 "Microsoft YaHei","Noto Sans CJK SC",sans-serif;color:#273748;background:#edf2f6;margin:0}}main{{max-width:1120px;margin:30px auto;padding:40px 48px;background:white}}h1{{font-size:30px;line-height:1.45}}h2{{font-size:23px;margin-top:2em;border-bottom:2px solid #d7e3ec;padding-bottom:8px}}.tag{{font-size:13px;letter-spacing:1px;color:#47758e}}.note{{border-left:5px solid #277da8;background:#edf6fb;padding:15px 22px}}.limit{{border-left:5px solid #cf5b36;background:#fff4ed;padding:15px 22px}}img{{width:100%;height:auto}}figure{{margin:24px 0}}figcaption,.small{{font-size:14px;color:#5d6b78}}a{{color:#176d9c}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;background:#f1f5f8;padding:20px;font:15px/1.7 Consolas,monospace}}table{{border-collapse:collapse;width:100%;font-size:14px;margin:16px 0}}td,th{{padding:9px 12px;border-bottom:1px solid #dae3eb;text-align:left}}th{{background:#edf5fa}}.scroll{{overflow:auto}}@media(max-width:700px){{main{{margin:0;padding:22px}}h1{{font-size:25px}}}}</style>
<main><div class="tag">Ni–14 at.% Al · 1000 K · 固定单颗粒重新统计</div>
<h1>终态最大颗粒 P4：<br>面积的 3/2 次方、半径立方与物理时间</h1>
<div class="note"><b>本次统计的是同一个 P4 颗粒。</b> 从原始 MOOSE Exodus 场重新提取面积，得到 A<sup>3/2</sup>(t) 和 r<sub>P4</sub><sup>3</sup>(t)。两组均固定跟踪初始中心位于 (72,72) nm 的 P4，它也是两组终态最大颗粒；没有按每一帧重新选择颗粒，也没有对全部颗粒取平均。</div>
<p><b>结果：</b>P4 的三次方增长可以在限定时间窗内近似为直线，有弹性组的增长更快。整个计算时段仍存在系统性弯曲、其他小颗粒溶解后的动力学变化，以及弹性组晚期趋于平台，因此不能认定一个恒定斜率覆盖全过程。</p>
<p>本次只改变后处理统计对象，重新使用同一组已完成的 MOOSE 场数据；没有改变成分、温度、自由能、弹性参数或重新求解相场。总时间保持 {m['time_calibration']['end_time_s']:.3f} s（{m['time_calibration']['end_time_hours']:.6f} h）。</p>

<h2>1. 颗粒身份与面积定义</h2>
{pic('P4_particle_identity', '箭头标识被持续跟踪的 P4。左右两组的初始场完全相同，终态形状可不同。')}
<pre>A(t) = ∫Ω_P4 h(η) dA，h(η) = η³(6η²−15η+10)
r_eq(t) = √[A(t)/π]
r_eq³(t) = [A(t)/π]^(3/2)
A^(3/2)(t) = π^(3/2) r_eq³(t) = {factor:.9f} r_eq³(t)
K_A = π^(3/2) K_r；两种纵轴的线性拟合 R² 完全相同。</pre>
<p>面积单位 nm²，A<sup>3/2</sup> 与 r³ 单位均为 nm³，但数值相差 π<sup>3/2</sup>。采用 r<sub>eq</sub> 表示等面积半径，即使弹性颗粒变成圆角方形也可比较。它不是由二维截面直接测得的三维体积等效半径。本次主要面积定义沿用前一轮的相分数积分与节点梯形积分法，只改变统计对象；另用 η=0.5 闭合轮廓的几何面积复核。</p>
{table(['算例','初始 r / nm','终态 r / nm','终态 A / nm²','终态 r³ / nm³','终态 A³ᐟ² / nm³'], summary)}

<h2 id="curves">2. 单颗粒 r³ 和 A³ᐟ² 随物理时间变化</h2>
{pic('P4_radius_cubed_vs_time', '实线为 P4 的面积等效半径立方；绿色虚线拟合 t*=5000–30000，黑色点线拟合 t*=5000–60000。竖线标记其他小颗粒的消失时间区间中点，不代表 P4 身份发生变化。左右纵轴范围不同，应使用斜率数值比较速率。')}
{pic('P4_area_3_2_vs_time', '按用户提出的面积 3/2 次方直接绘图。曲线形状与 r³ 相同，仅纵轴与斜率乘以 π³ᐟ²。')}
<p>最初短时间内的半径下降对应原初态的界面与成分松弛。沿用前次分析预先使用的 t*=5000 截止点，主窗口为 t*=5000–30000（0.231–1.382 h），该段两组接受的时间网格相同，且尚未发生其他颗粒消失；全程对照使用 t*=5000–60000（0.231–2.766 h）。同时保留从 t*=0 开始的拟合及全部固定分段结果，没有按最大 R² 搜索或挑选窗口。</p>

<h2 id="rates">3. 线性拟合、速率与窗口依赖性</h2>
<pre>r_P4³(t) ≈ b + K_r t，t 以秒计，截距 b 自由拟合。
K_r 的单位：nm³/s；乘 10⁻²⁷ 转换为 m³/s。
时间加权最小二乘：端点权重 Δt/2，内部权重 (tᵢ₊₁−tᵢ₋₁)/2。
此加权避免自适应输出点在溶解事件附近聚集而支配拟合。</pre>
{table(['算例','物理时间 / h','K_r / nm³·s⁻¹','K_r / m³·s⁻¹','K_A / nm³·s⁻¹','R²','RMSE/纵轴变化范围','窗口'], main_fits)}
<p>早期共同窗口中，弹性/无弹性单颗粒速率比为 <b>{ratio:.3f}</b>；全程对照窗口中为 <b>{whole_ratio:.3f}</b>。这些是指定窗口下的 P4 有效增长速率，不能据此给材料指定一个与时间窗口无关的“弹性加速倍数”。</p>
<p>早期无弹性拟合为 r³≈{no['primary_fit']['intercept_nm3']:.4f}+{no['primary_fit']['K_nm3_s']:.7f}t；有弹性为 r³≈{el['primary_fit']['intercept_nm3']:.4f}+{el['primary_fit']['K_nm3_s']:.7f}t。r³ 以 nm³、t 以秒计，截距是该窗口回归的外推量，不应当作 t=0 的实际半径立方。</p>
{pic('P4_linearity_diagnostics', '上排给出线性拟合残差，下排给出四个事先指定的不重叠时间段的有效斜率。横线表示拟合时间段，圆点标记其中点，均不是置信区间。残差的系统性弯曲和斜率变化说明全程速率并非常数。')}
<p>主窗口线性拟合的 R² 分别为 {no['primary_fit']['R2']:.4f} 和 {el['primary_fit']['R2']:.4f}；作为曲率诊断，同一窗口增加二次项后 R² 分别为 {no['primary_fit']['quadratic_R2']:.6f} 和 {el['primary_fit']['quadratic_R2']:.6f}。这说明较高线性 R² 仍可能伴随可分辨的弯曲。二次拟合仅用于检查曲率，没有用于更改相场模型或自由能。本数据是确定性、相邻帧相关的模拟轨迹，没有独立实验噪声模型，因此没有把普通回归置信区间或 p 值当作统计证据。</p>
{table(['算例','t* 窗口','实际物理时间 / h','K_r / nm³·s⁻¹','线性 R²','归一化 RMSE','二次 R²'], all_rows)}
<p class="small">population_N 开头的窗口仅按周围体系的颗粒消失事件划分，拟合对象始终为 P4。弹性组仅剩两颗粒后的末段趋于平台，尤其不适合解释成全程恒定速率。</p>

<h2>4. 与奥斯特瓦尔德熟化规律的关系</h2>
<div class="limit"><b>结论：局部近似三次方线性增长，与扩散熟化相容；全程恒定速率尚未成立。</b> 跟踪单颗粒消除了全体数平均分母变化造成的跳变，但不能消除基体过饱和度、邻近颗粒溶解及弹性相互作用随时间改变带来的真实动力学变化。</div>
<p>在理想三维球形、扩散控制的熟化描述中，单颗粒增长式可写为 d(rᵢ³)/dt = B(rᵢ/r_c−1)，B 为正的动力学系数。即使 B 固定，rᵢ/r_c 随时间变化仍会改变单颗粒 rᵢ³ 的斜率；只有该比值在所考察窗口近似恒定时，才自然得到近似直线。这里是依据 <a href="https://doi.org/10.1631/jzus.2005.B0705">Wang 等（2005），式 (6)–(8)</a>作出的解释，不是对本二维弹性计算直接套用的拟合方程。</p>
<p>本项目的二维、四颗粒有限体系以及非球形弹性颗粒限制了三维稳态 LSW 解释。现有数据确实表现出小颗粒溶解、大颗粒通过扩散获得溶质的熟化过程；P4 的 r³ 曲线可用于报告单颗粒窗口速率，但不能单凭其高 R² 证明已经进入普适的稳态粗化区间。</p>

<h2>5. 面积提取复核与时间依据</h2>
{pic('P4_area_definition_check', '相分数积分面积与 η=0.5 闭合轮廓面积得到的 P4 半径立方对照。两种面积定义均直接来自原始场。')}
{table(['算例','积分面积 K_r','轮廓面积 K_r','斜率相对变化','轮廓面积 R²'], area_check)}
<p>上表采用同一主窗口，K_r 单位均为 nm³/s。两种定义的半径全程最大相对差异约 {no['contour_max_radius_relative_difference_percent']:.3f}%；本次主窗口斜率变化小于 0.5%，因此“局部近似线性、全程斜率变化”的判断不依赖于这两种面积定义的选择。重新提取的 P4 半径与之前场数据提取结果最大差异小于 5×10⁻¹⁴ nm。</p>
<p>时间换算保持 t = {m['time_scale_s']:.11f} t* 秒，D̃γ≈{m['time_calibration']['gamma_chemical_interdiffusivity_m2_s']:.6e} m²/s。它来自 Engström–Ågren fcc 迁移率评估与当前 NIDEMO 热力学因子的组合，参数核对来源为 <a href="https://doi.org/10.3139/146.110428">Zhang 等（2010）表 3</a>。这是带有明确物理单位的文献模型时间，尚未经过该合金在 1000 K 的实验时效验证；该不确定性会统一缩放时间与所有速率，不影响 R²。</p>

<h2>6. 文件与复现</h2>
<p><a href="no_elastic_P4_history.csv">无弹性 P4 历史</a> · <a href="elastic_P4_history.csv">有弹性 P4 历史</a> · <a href="P4_growth_fits.csv">全部拟合结果</a> · <a href="metrics.json">定义、参数与校验</a>。CSV 包含模型时间、秒、小时、积分/轮廓面积、A³ᐟ²、半径及 r³。</p>
<p>本次图表均可从原始场重新生成。在 WSL 项目根目录执行：</p>
<pre>source ~/miniforge3-moose/etc/profile.d/conda.sh
conda activate moose-post
python scripts/analyze_single_particle_growth.py
python scripts/write_single_particle_report.py</pre>
<p>场文件仍使用已有 <a href="../elastic_seconds.e">有弹性结果</a> 与 <a href="../no_elastic_seconds.e">无弹性结果</a>；它们的时间已经是秒，本次没有更改这些文件。完整机械模型和前次不同统计量的结果见 <a href="../REPORT.html">弹性模型报告</a>。</p>
<p class="small">校验：两组终态最大颗粒均为 P4；所有保存帧内 P4 存在，轮廓闭合并留在原象限；重新积分与原半径一致；A³ᐟ²=π³ᐟ²r³、斜率换算及 R² 一致性检查通过。未进行新的网格/时间步收敛研究。</p>
</main></html>'''
(OUT/'REPORT.html').write_text(page, encoding='utf-8')
print(OUT/'REPORT.html')
