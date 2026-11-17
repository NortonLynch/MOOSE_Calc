# Ni–Al 实验熟化率资料

[REPORT.html](REPORT.html) 汇总邻近 Ni–14 at.% Al、1000 K 条件的已核查实验数据。
[experimental_rates.csv](experimental_rates.csv) 保留原始温度、成分、统计量和单位，并列出明确换算后的系数。
[sources.json](sources.json) 记录公开来源、原文定位、历史核查信息及可选本地文件的校验状态。

本目录的数据并不等同于精确 14 at.% Al、1000 K 的单颗粒实测斜率。原作者的主要统计量为
平均立方体半边长的三次方；不能直接当作本项目固定 P4 的面积等效半径三次方。

## 重现报告和单位换算

在仓库根目录使用 Python 3 运行，无需 MOOSE、Thermo-Calc、第三方 Python 包或本地论文 PDF：

```bash
python scripts/write_experimental_benchmark_report.py
```

脚本重现已经人工核对的表格转录及单位换算，输出本目录的 CSV、JSON 和 HTML。
它不会自动下载论文，也不会在生成报告时重新进行文献检索或视觉审查。历史检索日期保持为 2026-11-05。

## 原始论文与来源校验

第三方论文 PDF 仅作为可选本地参考保留，不随 Git 仓库分发。公开来源为：

- [Ardell 1969，Table I，PDF 第 22 页](https://www.phase-trans.msm.cam.ac.uk/2008/MPT/Part_3.pdf#page=22)。
- [Ardell 2011，Table 6，PDF 第 13 页](https://d-nb.info/1239497911/34#page=13)。
- [Zhu 等 2004，DOI](https://doi.org/10.1016/j.actamat.2004.02.032)。

来源清单中的 `sha256` 和 `bytes` 固定标识 2026-11-05 已经视觉核查的两份 PDF，
即使当前计算机没有这些文件，也会保留该历史依据。`visual_verification` 描述历史核查，
`local_file_check` 单独描述本次脚本执行时的文件状态：

- 本地 PDF 缺失：报告仍可生成，标记 `not_performed_local_pdf_absent`，`sha256_checked=false`，观测值为 `null`。
- 本地 PDF 存在且匹配历史大小与 SHA256：标记 `matches_historically_verified_pdf`。
- 本地 PDF 是空文件、内容不同或不是普通文件：脚本报错，并在写入数据与报告之前停止，不会更新历史校验值。

可选本地文件位置为 `sources/Ardell_1969_MPT_Part3.pdf` 和 `sources/Ardell_2011.pdf`。
公开来源日后可能更换文件版本；新的字节内容需要另行核查，不能仅因为 URL 相同就覆盖历史来源身份。

988 K 的 2011 系数是历史实验数据的重新拟合；Zhu 2004 引用的 4600 nm³/h 是平均边长立方系数，
换为平均半边长立方须除以 8。来源表及报告保持这些区分。
