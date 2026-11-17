"""Archive verified experimental coarsening slopes and their original definitions.

Numerical entries were transcribed from the original author's 1969 Table I
and the 2011 reanalysis Table 6, with visual inspection of both tables.
No thermodynamic/diffusion prediction is labeled as an experimental slope.

The third-party PDFs are optional local references and are not distributed with
the repository. Their pinned hashes identify the copies visually checked in
2026-11-05. Missing PDFs do not prevent regenerating the transcribed tables;
present PDFs must match those historical identities before output is written.
"""
from pathlib import Path
import csv, json, hashlib, math

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/ni_al_elastic/experimental_benchmarks'
OUT.mkdir(parents=True, exist_ok=True)
M_AL, M_NI = 26.9815385, 58.6934


def at_percent(wt):
    return 100*(wt/M_AL)/(wt/M_AL+(100-wt)/M_NI)


sources = {
    'Ardell1969': {
        'title': 'Experimental Confirmation of the Lifshitz-Wagner Theory of Particle Coarsening',
        'author': 'A. J. Ardell', 'year': 1969,
        'publication': 'The Mechanism of Phase Transformations in Crystalline Solids, Institute of Metals, Monograph and Report Series 33, pp.111-116',
        'source_kind': 'Original experimenter reporting a summary of his experimental measurements',
        'url': 'https://www.phase-trans.msm.cam.ac.uk/2008/MPT/Part_3.pdf',
        'local_pdf': 'sources/Ardell_1969_MPT_Part3.pdf',
        'sha256': 'b95e147d8ea2ec24709bf71430b55ecd78d616a14b1c7657764790ac62706871',
        'bytes': 35045906,
        'location': 'Printed p112 Table I and preceding paragraph; PDF page22 (1-based). Composition assignment also checked in printed p113 Fig1/2, PDF page23.',
        'quantity': 'd[(mean cube edge length / 2)^3]/dt',
        'measurement_method': 'Transmission electron microscopy particle sizes',
        'fit_time_window': 'Table I does not give exact start/end times of each regression. Fig2 PSD aging times are not asserted to be the regression window.',
        'visual_verification': True,
        'visual_verification_date': '2026-11-05',
        'notes': '625C and 715C entries are averages for 6.35 and 6.71 wt.% Al alloys (table footnote). 750C and 775C are the 6.71 wt.% alloy. No uncertainty is printed for k in this table.'},
    'Ardell2011': {
        'title': 'A1-L12 interfacial free energies from data on coarsening in five binary Ni alloys, informed by thermodynamic phase diagram assessments',
        'author': 'A. J. Ardell', 'year': 2011,
        'publication': 'Journal of Materials Science 46, 4832-4849',
        'source_kind': 'Published n=3 regression of earlier experimental particle-size data, not an independent new experiment',
        'doi': 'https://doi.org/10.1007/s10853-011-5395-x',
        'url': 'https://d-nb.info/1239497911/34',
        'local_pdf': 'sources/Ardell_2011.pdf',
        'sha256': '9f7754a3856aa56cdf5bf5abf24d205d4a72e9168c9997424d2d604d6113800d',
        'bytes': 484750,
        'location': 'Printed p4844 Table6, first two rows; PDF page13 (1-based)',
        'quantity': 'd[mean radius^3]/dt, Table6 n=3 (LSW analysis)',
        'fit_time_window': 'Not separately specified in Table6',
        'visual_verification': True,
        'visual_verification_date': '2026-11-05',
        'notes': 'Do not substitute Table5 non-integer-exponent coefficients. Initial alloy composition is not tabulated for each row in Table6; keep the source label Ni-Al rather than assigning exact 14 at.%.'},
    'Zhu2004': {
        'title': 'Three-dimensional phase-field simulations of coarsening kinetics of gamma-prime particles in binary Ni-Al alloys',
        'authors': 'J. Z. Zhu, T. Wang, A. J. Ardell, S. H. Zhou, Z. K. Liu, L. Q. Chen',
        'year': 2004, 'publication': 'Acta Materialia 52, 2837-2845',
        'doi': 'https://doi.org/10.1016/j.actamat.2004.02.032',
        'url': 'https://citeseerx.ist.psu.edu/document?doi=497199be7d4e6c92bb78cac9d55c995838742321&repid=rep1&type=pdf',
        'source_kind': 'Phase-field article quoting an experimental result; distinguished from its simulated result',
        'location': 'p2843 section3.3; Eq12 and Fig4 (p2844); experimental reference39',
        'quoted_experimental_source': 'Ardell and Nicholson, J Phys Chem Solids 27 (1966), starting p1793, DOI10.1016/0022-3697(66)90110-7',
        'reported_composition_at_percent_Al': 13.8,
        'reported_temperature_K': 1023,
        'quantity': 'd[mean cube edge length^3]/dt',
        'experimental_value_nm3_h': 4600,
        'simulation_value_nm3_h_excluded': 5570,
        'notes': 'The 13.8 at.% label differs from 6.71 wt.% (~13.53 at.%) in the original author summary. Preserve both labels; do not silently relabel as exact Ni-14 at.% Al.'}
}


def verify_optional_pdf(source):
    """Keep historical source identity distinct from this run's file check."""
    source['file_identity_basis'] = 'Pinned SHA256 and size of the PDF visually inspected on 2026-11-05; not computed from an absent file.'
    source['visual_verification_scope'] = 'Historical visual inspection of the identified PDF and cited table; not repeated by this report generator.'
    p = OUT/source['local_pdf']
    if not p.exists():
        source['local_file_check'] = {
            'status': 'not_performed_local_pdf_absent',
            'local_pdf_present': False,
            'sha256_checked': False,
            'observed_sha256': None,
            'observed_bytes': None,
        }
        return
    if not p.is_file():
        raise ValueError(f'Optional PDF path is not a regular file: {p}')
    size = p.stat().st_size
    digest = hashlib.sha256(p.read_bytes()).hexdigest()
    if size != source['bytes'] or digest != source['sha256']:
        raise ValueError(
            f'Local PDF differs from the historically verified source: {p}. '
            f'Expected {source["bytes"]} bytes, SHA256 {source["sha256"]}; '
            f'observed {size} bytes, SHA256 {digest}. '
            'No report was regenerated; check this local copy against the public source.'
        )
    source['local_file_check'] = {
        'status': 'matches_historically_verified_pdf',
        'local_pdf_present': True,
        'sha256_checked': True,
        'observed_sha256': digest,
        'observed_bytes': size,
    }


for s in sources.values():
    if 'local_pdf' in s:
        verify_optional_pdf(s)

rows = []
for temp_C, raw, comp in [(625, 2.07e-24, '6.35 and 6.71 wt.% Al; averaged'),
                           (715, 7.01e-23, '6.35 and 6.71 wt.% Al; averaged'),
                           (750, 1.60e-22, '6.71 wt.% Al'),
                           (775, 3.61e-22, '6.71 wt.% Al')]:
    rows.append(dict(source_id='Ardell1969', classification='original_author_experimental_summary',
                     temperature_C=temp_C, temperature_K=temp_C+273.15, alloy_as_reported=comp,
                     original_quantity='(mean_cube_edge/2)^3', original_K=raw, original_unit='cm^3/s',
                     K_radius_convention_nm3_s=raw*1e21, K_radius_convention_m3_s=raw*1e-6,
                     uncertainty_nm3_s='', uncertainty_note='not reported in TableI',
                     single_particle=False, exact_Ni14atAl_1000K_match=False,
                     fit_window='not specified in TableI', location='printed p112 TableI; PDF p22'))
for temp_K, raw, err in [(898, 2.1379e-30, .2565e-30), (988, 7.2526e-29, .4011e-29)]:
    rows.append(dict(source_id='Ardell2011', classification='experimental_data_reanalysis_n3',
                     temperature_C=temp_K-273.15, temperature_K=temp_K,
                     alloy_as_reported='Ni-Al; composition not repeated for each Table6 row',
                     original_quantity='mean_radius^3', original_K=raw, original_unit='m^3/s',
                     K_radius_convention_nm3_s=raw*1e27, K_radius_convention_m3_s=raw,
                     uncertainty_nm3_s=err*1e27, uncertainty_note='reported uncertainty; no 95% CI asserted',
                     single_particle=False, exact_Ni14atAl_1000K_match=False,
                     fit_window='not specified in Table6', location='printed p4844 Table6; PDF p13'))
rows.append(dict(source_id='Zhu2004', classification='quoted_experiment_with_explicit_geometry_conversion',
                 temperature_C=1023-273.15, temperature_K=1023,
                 alloy_as_reported='Ni-13.8 at.% Al in Zhu; original summary6.71wt.%Al',
                 original_quantity='mean_cube_edge^3', original_K=4600, original_unit='nm^3/h',
                 K_radius_convention_nm3_s=4600/3600/8,
                 K_radius_convention_m3_s=4600/3600/8*1e-27,
                 uncertainty_nm3_s='', uncertainty_note='no experimental uncertainty quoted for this number',
                 single_particle=False, exact_Ni14atAl_1000K_match=False,
                 fit_window='not independently established here', location='p2843 section3.3; Eq12, Fig4'))
with (OUT/'experimental_rates.csv').open('w', newline='', encoding='utf-8-sig') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)
index = dict(search_date='2026-11-05', target={'Al_at_percent':14, 'temperature_K':1000,
            'phases':'gamma matrix / gamma-prime precipitates'}, exact_match_status='not located in the verified sources; not a claim that none exists',
            source_file_policy='Third-party PDFs remain optional local references and are not distributed in Git. Public links and historical file identities are retained. Missing local files are explicitly not checked; mismatching local files stop generation.',
            sources=sources, wt_to_at_conversion={'Al_molar_mass':M_AL, 'Ni_molar_mass':M_NI,
            '6.35wt_at_percent':at_percent(6.35), '6.71wt_at_percent':at_percent(6.71)},
            conversion_checks={'cm3_to_nm3':1e21, 'cm3_to_m3':1e-6, 'nm3_h_to_nm3_s':1/3600,
            'cube_edge_cubed_to_half_edge_cubed':1/8,
            'Zhu_experiment_half_edge_cubed_nm3_s':4600/3600/8},
            exclusions=['No diffusion-based theoretical K is recorded as an experimental measurement.',
                        'No interpolation or extrapolation to1000K is labeled as a measured point.',
                        'No population-average experimental slope is labeled as a fixed-single-particle slope.',
                        'Non-integer growth exponent coefficients have different dimensions and are not mixed with n=3 slopes.'])
(OUT/'sources.json').write_text(json.dumps(index, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')


def table(head, entries):
    return '<div class="scroll"><table><thead><tr>'+''.join(f'<th>{v}</th>' for v in head)+'</tr></thead><tbody>'+''.join('<tr>'+''.join(f'<td>{v}</td>' for v in row)+'</tr>' for row in entries)+'</tbody></table></div>'


original = [[f'{r["temperature_C"]} °C / {r["temperature_K"]:.2f} K',
             f'{r["original_K"]:.3e}', f'{r["K_radius_convention_nm3_s"]:.6g}',
             f'{r["K_radius_convention_m3_s"]:.3e}',
             '两合金平均值' if 'averaged' in r['alloy_as_reported'] else '6.71 wt.% Al'] for r in rows[:4]]
reanalyzed = [[r['temperature_K'], f'{r["K_radius_convention_nm3_s"]:.7f} ± {r["uncertainty_nm3_s"]:.7f}',
               f'{r["K_radius_convention_m3_s"]:.5e}'] for r in rows[4:6]]
page = f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Ni–Al 实验熟化率：来源与单位核查</title><style>body{{font:16px/1.85 "Microsoft YaHei",sans-serif;color:#263748;background:#eef3f6;margin:0}}main{{max-width:1080px;background:white;margin:30px auto;padding:38px 46px}}h1{{font-size:29px;line-height:1.5}}h2{{font-size:22px;margin-top:1.8em;border-bottom:2px solid #d3e1e9;padding-bottom:8px}}a{{color:#176f99}}.note{{border-left:5px solid #cb642e;background:#fff4eb;padding:15px 20px}}table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{border-bottom:1px solid #dce5eb;padding:9px 12px;text-align:left}}th{{background:#edf5f9}}.scroll{{overflow:auto}}pre{{white-space:pre-wrap;background:#f1f5f8;padding:18px;font:14px/1.7 Consolas,monospace}}.small{{font-size:14px;color:#566c7c}}</style><main>
<h1>Ni–Al γ′实验熟化率：<br>原作者数据、半径定义与邻近温度</h1>
<p>检索目标：Ni–14 at.% Al、1000 K（726.85°C），γ 基体中 γ′析出相熟化。检索与核对日期：2026-11-05。</p>
<div class="note"><b>已找到可核实的邻近条件实验斜率，尚未找到精确匹配 14 at.% Al、1000 K 的直接实测点。</b> 下列数值来自实验粒径数据及其已发表拟合；没有用扩散公式或模拟结果冒充实验值，也没有将温度插值当作实测。</div>
<h2>1. 原作者直接汇总的实验数据</h2>
<p>A. J. Ardell（1969），<i>Experimental Confirmation of the Lifshitz–Wagner Theory of Particle Coarsening</i>，收于 <i>The Mechanism of Phase Transformations in Crystalline Solids</i>，Institute of Metals，Monograph and Report Series 33，pp.111–116。<a href="https://www.phase-trans.msm.cam.ac.uk/2008/MPT/Part_3.pdf#page=22">公开全文及原表（PDF 第22页，印刷 p112，Table I）</a>。</p>
<p>原文使用 TEM 测量颗粒尺寸，表中斜率为 <b>K=d[(ā/2)³]/dt</b>；ā 是平均立方体边长。625°C 与 715°C 的数值是 6.35、6.71 wt.% Al 两合金的平均值；较高两温度对应 6.71 wt.% Al。表中未列 K 的误差或每条回归的精确起止时间。</p>
{table(['温度','原表 K / cm³·s⁻¹','K / nm³·s⁻¹','K / m³·s⁻¹','成分说明'], original)}
<p>按原子量换算，6.35 wt.% Al≈{at_percent(6.35):.2f} at.% Al；6.71 wt.% Al≈{at_percent(6.71):.2f} at.% Al。它们接近本项目成分，但不是精确 14 at.% Al。</p>
<h2>2. 后续对实验数据的 n=3 重拟合</h2>
<p>Ardell（2011），<a href="https://doi.org/10.1007/s10853-011-5395-x">J. Mater. Sci. 46, 4832–4849</a>，<a href="https://d-nb.info/1239497911/34#page=13">公开全文 Table 6（印刷 p4844，PDF 第13页）</a>列出了平均半径立方的实验数据拟合系数。它们是相关历史数据的重新分析，不是两套新增独立实验；Table 6 未逐行重复标称合金成分。</p>
{table(['文献温度 / K','K / nm³·s⁻¹（文献所报不确定度）','K / m³·s⁻¹'], reanalyzed)}
<p>988 K 的 0.072526±0.004011 nm³/s 与原作者汇总的约 988 K、0.0701 nm³/s 接近。这里使用 n=3 的 Table 6，不能混入该文 Table 5 的非整数指数系数。所列误差沿用文献表值，不额外解释为 95% 置信区间。</p>
<h2>3. Zhu 2004 实验对照的核对</h2>
<p><a href="https://doi.org/10.1016/j.actamat.2004.02.032">Zhu 等（2004），Acta Materialia 52, 2837–2845</a>，p2843 §3.3 引用的实验条件写作 Ni–13.8 at.% Al、1023 K，实验斜率为 <b>4600 nm³/h</b>，但统计量是 <b>ā³</b>。其原实验出处为 Ardell 与 Nicholson（1966），<a href="https://doi.org/10.1016/0022-3697(66)90110-7">J. Phys. Chem. Solids 27，起始页1793</a>。原始1966论文全文本次未直接逐页审查；已通过原作者1969原表交叉核实。</p>
<pre>K_a = 4600 / 3600 = 1.277778 nm³/s
K_(a/2) = K_a / 8 = 0.159722 nm³/s ≈ 原作者表中的 0.160 nm³/s</pre>
<p>Zhu 文中的 5570 nm³/h 是模拟斜率，未列入实验数据表。Zhu 的 13.8 at.% 成分标签与原作者 6.71 wt.%（约13.53 at.%）也有差异，保留来源差异，不将它们改写成精确 Ni–14 at.% Al、1000 K。</p>
<h2>4. 与本项目 P4 结果如何对照</h2>
<p>当前 P4 统计量是固定单颗粒的面积等效半径 r_A=√(A/π)，实验表统计的是许多颗粒的平均半边长。二者的几何定义和统计对象都不同；不能仅因斜率接近，就宣称模型已获实验验证，或直接用一个比例重新标定扩散时间。</p>
<p>仅在理想正方形投影 A=a²、并保持相同总体平均方式时，r_A=a/√π，才有 K_rA=8K_(a/2)/π³ᐟ²。以 1023 K 的原表为例，得到约 {8*.160/math.pi**1.5:.5f} nm³/s。这个数是几何换算示例，不是额外的实验测量，更不是单颗粒 P4 在1000 K 的实测斜率。真实颗粒圆角、截面位置和成像投影都会影响这种换算。</p>
<p>用于当前项目的邻近条件参考，优先保留约 988 K 的 0.0701 nm³/s 和约 1023 K 的 0.160 nm³/s，连同原成分与半径定义一起使用。若需要1000 K数值，必须另行作有假设的温度插值或补充实验；本次没有把这种估计记入实测表。</p>
<h2>5. 数据文件与边界</h2>
<p><a href="experimental_rates.csv">可复用数据表</a> · <a href="sources.json">来源、定位、单位与文件校验信息</a>。保存原始单位和原始统计量，同时提供显式换算列。未更改现有 MOOSE 输入、时间标尺或模拟结果。</p>
<p>第三方论文 PDF 作为可选本地参考保留，不随 Git 仓库分发。原表公开链接已在上文给出；<a href="sources.json">来源清单</a>保留 2026-11-05 已视觉核查副本的 SHA256 和文件大小，以标识当时的依据。重新生成报告不需要本地 PDF：缺失时明确记录本次文件校验未执行；存在时必须与历史副本的大小及 SHA256 一致，否则停止生成。脚本不会把空文件、缺失文件或重新计算单位表述成重新审查了论文。复现命令见<a href="README.md">本目录说明</a>。</p>
<p>另需避免将 r<sup>2.4</sup> 对时间的系数混入 r³ 数据。例如 <a href="https://doi.org/10.1007/s11661-021-06440-0">Ardell（2021）</a>的大部分再分析采用 n=2.4，单位为 m<sup>2.4</sup>/s，不能直接当作 m³/s。</p>
<p class="small">2026-11-05 的核验记录：原作者 Table I 的单位、脚注、粒径定义和成分；2011 Table 6 的指数与不确定度；Zhu 实验/模拟系数区分；cm³/s、m³/s、nm³/h 的换算。报告生成脚本重现已核对的数值与单位换算，当前本地文件的校验状态另列于 sources.json。准确同条件实验未检索到，不代表断言不存在。</p>
</main></html>'''
(OUT/'REPORT.html').write_text(page, encoding='utf-8')
assert abs(rows[1]['K_radius_convention_nm3_s']-.0701) < 1e-12
assert abs(rows[2]['K_radius_convention_nm3_s']-.160) < 1e-12
assert abs(rows[5]['K_radius_convention_nm3_s']-.072526) < 1e-12
assert abs(rows[-1]['K_radius_convention_nm3_s']/.160-1) < .002
print(json.dumps({'rows':len(rows), 'unit_checks':'PASS', 'Al_at_percent_from_6.35wt':at_percent(6.35),
                   'Al_at_percent_from_6.71wt':at_percent(6.71),
                   'source_file_checks':{key:s['local_file_check']['status'] for key,s in sources.items() if 'local_pdf' in s},
                   'report':str(OUT/'REPORT.html')}, indent=2))
