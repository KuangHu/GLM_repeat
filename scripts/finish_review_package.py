#!/usr/bin/env python3
"""Complete the offline review package: entry page, report, manifests, checksums.

build_review_package.py writes the candidate pages, figures, sequences and tables.
This stage adds everything a reviewer needs to open the package cold on a laptop
with no network and no cluster access:

  START_HERE.html   Chinese entry page; relative links only
  README_zh.md      the same orientation in plain Markdown
  report/           this round's report plus the frozen upstream reports
  manifests/        the frozen data-selection and threshold manifests
  MANIFEST.sha256   checksum of every file in the package

Nothing here re-derives or re-tunes anything: it copies frozen artefacts and
reports numbers already computed by the pipeline.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import csv
import shutil
from pathlib import Path

SECTION8 = '## 8. 对浏览册 69 个候选的审阅（助手判断，非人工确认）\n\n判断写在 `review/assistant_review.tsv` 的 `AUTO_REVIEW_*` 列。\n**`HUMAN_*` 列仍然全部留空** —— 机器判断不应占据留给人工阅读的字段。\n每条判断都由明确规则产生，并附带触发它的证据（`AUTO_REVIEW_reasoning`），\n因此任何一条都可以被逐项质疑。`cannot_judge` 是一等结论，不会为了好看而回避。\n\n{rnd}\n{enr}\n**两层不可合并。** 只有随机层描述"审阅者会遇到什么"；富集层是按规则挑选的导览。\n\n### 8.1 本轮审阅得到的主要结论\n\n**(1) 宿主路线上约一半的候选是 tRNA 基因。**\n\n| 运行 | 注释 tRNA 数 | 候选中与 tRNA 相交 |\n|---|---:|---:|\n{at}\n\n两轮宿主面板分别为 52.4% 与 52.6%，两轮噬菌体为 17.0% 与 11.9%，各自高度一致。\n这意味着在当前阈值下，**宿主路线的 detector 在很大程度上是一个 tRNA / 结构性 ncRNA 检测器**。\ntRNA 的 masked predictability 高，是因为它们受到强烈的结构与序列约束，\n**这与"存在序列重复"是两回事**。\n\n**(2) `unresolved` 主要不是失败，而是模型类别不匹配。**\n宿主 round-2 的 89 个 `no_model` 候选中有 60 个是 tRNA。\nungapped repeat explainer 无法表示结构性 ncRNA 基因，这是类别错配，不是漏解释。\n报告第 2 节把 `unresolved` 当作"未解释"统计，这个口径需要在下一轮修正。\n\n**(3) CRISPR 阵列被完整召回，并额外提示 4 个未被专用工具调用的位点。**\npilercr 在宿主 round-2 注释出 9 个 CRISPR 阵列，**9 个全部被候选覆盖**。\n另有 4 个 *Aquifex aeolicus* 候选携带该基因组的 CRISPR repeat consensus\n（`GTTCCTAATGTACCGTGTGGAGTTGAAACC`）却落在所有 pilercr 阵列之外，\n其中一个（`NC_000918.1:263146-263290`）出现在随机层。\n这类位点可能是退化或残缺的阵列，是本轮最值得单独核查的结果。\n**这不是"发现了新的生物学系统"**：它们是专用工具未调用、而本流程标记出来的已知类型位点。\n\n**(4) 12-nt 种子的偶然性必须逐例校准。**\n对每个候选区域做保二核苷酸组成的洗牌（300 次）后比较重复 12-mer 数量：\n`candidate_associated` 的模型几乎全部显著（p<0.05），\n而 `no_model` 候选中只有 1/30 显著 —— 它们确实没有超出偶然的短词重复。\n但 `flank_only` 的 `repeat_pair` 只有一半显著，且其 copies 全部落在候选之外，\n**这类模型不构成关于该候选本身的证据**，一律记为 `cannot_judge`。\n\n**(5) 低复杂度必须同时看 consensus，不能只看包络。**\n包络常有数千 nt，会把低复杂度的 consensus 稀释到看不见。\n例如 `NC_048045.1:37076-40240` 的 consensus 是 `GGCGGCGGCGGC`（熵 0.92），\n而整个包络的熵是 1.86。本轮已把 consensus 熵加入低复杂度判定。\n\n### 8.2 这轮审阅没有做的事\n\n- 没有重新调参：detector、阈值、explainer 全部保持冻结。\n- 没有把 tRNA / CRISPR 注释反馈进候选判定，它们只用于解释与展示。\n- 没有宣称任何候选是经过确认的生物学发现。\n'

ROUTE_LABEL = {'r1_host': '宿主 round-1', 'r1_phage': '噬菌体 round-1',
               'r2_host': '宿主 round-2', 'r2_phage': '噬菌体 round-2'}


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with Path(path).open('rb') as fh:
        for b in iter(lambda: fh.read(chunk), b''):
            h.update(b)
    return h.hexdigest()


def load(run_dir):
    R = Path(run_dir)
    j = lambda p: json.loads((R / p).read_text()) if (R / p).exists() else {}
    return {'calls': j('calls/summary.json'), 'explain': j('explain_full/summary.json'),
            'policy': j('policy_attribution/summary.json'),
            'probe': j('no_model_probe/summary.json')}


def md_table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |',
           '|' + '|'.join(['---'] * len(header)) + '|']
    out += ['| ' + ' | '.join(str(x) for x in r) + ' |' for r in rows]
    return '\n'.join(out)


def build_review_section(rev, ann):
    if not rev:
        return ''
    def blk(d, title):
        o = [f'**{title}（n={d["n"]}）**', '']
        for field, label in (('repeat_evidence', '重复证据'),
                             ('explanation', '解释完整性'),
                             ('followup', '后续价值')):
            items = ', '.join(f'{k} {v}' for k, v in d[field].items())
            o.append(f'- {label}：{items}')
        o.append(f"- 与注释 tRNA 相交：{d['tRNA_overlapping']}；"
                 f"位于 CRISPR 阵列内：{d['CRISPR_overlapping']}；"
                 f"携带 CRISPR repeat 但不在阵列内：{d['CRISPR_repeat_off_array']}")
        return chr(10).join(o) + chr(10)
    at = chr(10).join(
        f"| {k} | {v['n_trna']} | {v['cand_trna_pct']} |" for k, v in ann.items())
    return (chr(10) + chr(10) + SECTION8).format(
        rnd=blk(rev['random_stratum'], '随机层'),
        enr=blk(rev['enriched_stratum'], '富集层'), at=at)


def build_report(stats, sel_prov, rev=None, ann=None):
    g = lambda lab, sec, k, d='—': stats[lab][sec].get(k, d)
    labs = [l for l in ('r1_host', 'r1_phage', 'r2_host', 'r2_phage') if l in stats]
    fmt = lambda v, n=1: f'{v:,.{n}f}' if isinstance(v, float) else f'{v:,}' if isinstance(v, int) else v

    scan = md_table(
        ['运行', '候选数', '唯一扫描 nt', '候选/Mb', '目标区 候选/Mb', 'CDS 策略丢失包络'],
        [[ROUTE_LABEL[l], fmt(g(l, 'calls', 'n_candidates')),
          fmt(g(l, 'calls', 'unique_scanned_nt')),
          fmt(g(l, 'calls', 'all_candidates_per_scanned_Mb')),
          fmt(g(l, 'calls', 'target_candidates_per_target_Mb')),
          fmt(g(l, 'calls', 'audit_full_without_policy_overlap'))] for l in labs])

    expl = md_table(
        ['运行', 'candidate_associated', 'flank_only', 'unresolved', '主模型', '假设总数'],
        [[ROUTE_LABEL[l],
          g(l, 'explain', 'status_counts', {}).get('candidate_associated', '—'),
          g(l, 'explain', 'status_counts', {}).get('flank_only', '—'),
          g(l, 'explain', 'status_counts', {}).get('unresolved', '—'),
          g(l, 'explain', 'n_primary_models'), g(l, 'explain', 'n_hypotheses')]
         for l in labs])

    pol = md_table(
        ['运行', 'full copies', 'policy copies', '精确保留', '丢失', '新增', '受影响候选'],
        [[ROUTE_LABEL[l], fmt(g(l, 'policy', 'full_copies')), fmt(g(l, 'policy', 'policy_copies')),
          fmt(g(l, 'policy', 'retained_exact')), fmt(g(l, 'policy', 'lost_exact')),
          fmt(g(l, 'policy', 'gained_exact')), g(l, 'policy', 'affected_candidates')]
         for l in labs])

    prb = md_table(
        ['运行', '无模型候选', 'no_repeated_seed', 'seeds_without_extension',
         'ungapped_similarity', 'gapped（explainer 无法表示）', '阳性对照灵敏度'],
        [[ROUTE_LABEL[l], g(l, 'probe', 'n_no_model_candidates'),
          g(l, 'probe', 'verdicts', {}).get('no_repeated_seed', 0),
          g(l, 'probe', 'verdicts', {}).get('seeds_without_extension', 0),
          g(l, 'probe', 'verdicts', {}).get('ungapped_similarity', 0),
          g(l, 'probe', 'verdicts', {}).get('gapped_similarity_explainer_cannot_represent', 0),
          f"{g(l, 'probe', 'positive_control_detected')}/{g(l, 'probe', 'positive_control_n')}"]
         for l in labs])

    total_nt = sum(g(l, 'calls', 'unique_scanned_nt', 0) for l in labs)
    total_c = sum(g(l, 'calls', 'n_candidates', 0) for l in labs)
    r2_nt = sum(g(l, 'calls', 'unique_scanned_nt', 0) for l in labs if l.startswith('r2'))
    r2_c = sum(g(l, 'calls', 'n_candidates', 0) for l in labs if l.startswith('r2'))

    review_section = build_review_section(rev, ann) if rev else ''
    return f"""# GLM_repeat round-1 + round-2 扫描报告（供人工审阅）

**本轮没有对 detector 做任何重新调参。** 阈值、打分器、explainer 全部沿用冻结版本；
round-2 只是在新数据上运行同一套冻结流程。

## 0. 这份报告报告什么、不报告什么

- 报告的是**候选产出率（candidate yield）**，即每 Mb 唯一扫描序列产生多少候选。
- **不是 FDR，不是精确率，不是生物学概率。** 没有任何候选经过人工确认。
- explainer 的 `candidate_associated` 只表示"存在一个能对齐上的重复模型"，
  不表示该重复是真实的生物学结构，也不表示它有功能。
- probe 的阴性结果只表示"在该 probe 的参数下未检出合格的局部相似性"，
  **不表示该区域没有重复**。

## 1. 扫描规模与产出

round-2 新增 {r2_nt/1e6:.2f} Mb（{r2_c} 个候选），四次运行合计
{total_nt/1e6:.2f} Mb、{total_c} 个候选。

{scan}

两次独立的宿主面板给出 {g('r1_host','calls','all_candidates_per_scanned_Mb'):.1f} 与
{g('r2_host','calls','all_candidates_per_scanned_Mb'):.1f} 候选/Mb，两次噬菌体面板给出
{g('r1_phage','calls','all_candidates_per_scanned_Mb'):.1f} 与
{g('r2_phage','calls','all_candidates_per_scanned_Mb'):.1f} 候选/Mb。
噬菌体路线的产出率约为宿主路线的 2.7 倍，这在两轮中都成立。

**这不是一个生物学结论。** 噬菌体基因组更短、编码密度更高、窗口覆盖方式不同，
而阈值是在宿主背景上标定的；产出率差异至少有一部分来自这些非生物学因素。

## 2. explainer（v0.3，冻结）

{expl}

`unresolved` 在两条路线上都占三到六成。这些候选的 detector 信号存在，但
ungapped explainer 没能给出模型。第 4 节的 probe 用于区分其中哪些仍带有
indel 型相似性。

## 3. CDS 策略的影响

{pol}

**CDS 策略的影响是基因组依赖的，不是系统性的。** 宿主 round-2 丢失 0 个包络、
0 条 copy 记录；噬菌体 round-2 丢失 30 个包络、758 条 copy 记录。
归因表（`tables/` 与各运行的 `policy_attribution/attribution.tsv`）把差异拆成
两个独立的轴：`upstream_trigger`（哪一级流水线发生了变化）与
`representation_outcome`（这条 copy 是否仍以某种形式出现）。
只有 `no_corresponding_copy` 表示证据在 policy 输出中确实缺失，
**即便如此，它也是计算证据的差异，不是经过验证的生物学漏检。**

## 4. 无模型候选的 indel 感知 probe（仅诊断）

{prb}

probe 不参与候选判定，也不回写任何冻结输出。阳性对照是强制的：
它在 explainer **已经**建模的候选上运行，用来证明这个 probe 本身有灵敏度。

分层验证（`report/PROBE_FROZEN.md`）已经说明，该 probe 在 8–12 nt 单元上的
灵敏度只有 0.20–0.22，并且在典型区域长度上有 13–26% 的假阳性率。因此：

- **probe 阳性 ≠ 可信重复**；典型长度下相当一部分是偶然相似。
- **probe 阴性 ≠ 假阳性候选**；短单元重复大多会被漏掉。

## 5. 浏览册的抽样方式

浏览册共 {sel_prov['selected_total']} 个候选，分两层，**必须分开解读**：

- **随机层（{sel_prov['selected_random']} 个）**：从 round-2 两条路线各随机抽
  {sel_prov['n_random_per_route']} 个，种子 {sel_prov['seed']}，完整抽样名单保存在
  `manifests/selection_provenance.json`。**只有这一层支持无偏的产出描述。**
  被随机抽中的候选即使同时被富集规则选中，仍然留在随机层的分母里。
- **富集层（{sel_prov['selected_enriched_only']} 个）**：按显式命名的规则挑选
  （清晰结构、多拷贝、不同基因邻域、高分无模型、不同谱系）。
  这一层**刻意**使用了 detector 分数、重复结构和邻近基因信息，
  它是一次导览，不是样本；**不能从中计算任何产出率、精确率或 FDR**。

## 6. 已知局限（必须随结论一起保留）

1. **VIPR 漏检**：已知的 VIPR 面板漏检在本轮仍未解决，round-2 没有针对它做任何修改。
2. **短重复 probe 局限**：见第 4 节，8–12 nt 单元灵敏度 0.20–0.22。
3. **环状基因组原点跨越未检测**：跨原点的重复不会被发现。
4. **CDS 注释无功能信息**：基因是 pyrodigal 从头预测，所有 `product` 都是
   "hypothetical protein"。浏览册里的"基因邻域"分类是**结构性的**
   （是否与 CDS 重叠、两侧基因朝向），不含任何功能判断。
5. **传递性合并**：合并后的包络可能包含多个生物学位点。
6. **阈值标定的迁移性**：阈值在旧的全基因组背景上标定，不保证在 CDS 选择后
   仍给出同样的每 Mb 产出。
7. **未分类候选不是新生物学系统**：没有归入已知结构的候选只是"未解释"，
   不得描述为新发现的生物学系统。

## 7. 本轮发现的一个流水线缺陷（已记录，未改冻结代码）

`calls/audit_full.tsv` 与 `calls/audit_policy.tsv` 带有 `candidate_id` 列，
但该列在所有行中都是空的（例如宿主 round-2 为 156/156 行为空）。
explainer 从坐标重新推导 ID，所以一直没有暴露；邻域导出脚本会校验该列，因此报错。
本轮的处理是改用 `calls/candidates.tsv`（ID 已填充，坐标集合与 `audit_full.tsv`
逐行一致，四次运行均已核对）。**冻结代码没有修改**，该缺陷留待下一轮修复。
{review_section}"""


def start_here(stats, sel_prov, n_book, n_all):
    labs = [l for l in ('r1_host', 'r1_phage', 'r2_host', 'r2_phage') if l in stats]
    rows = ''.join(
        f"<tr><td>{ROUTE_LABEL[l]}</td>"
        f"<td>{stats[l]['calls'].get('n_candidates', '—'):,}</td>"
        f"<td>{stats[l]['calls'].get('unique_scanned_nt', 0)/1e6:.2f}</td>"
        f"<td>{stats[l]['calls'].get('all_candidates_per_scanned_Mb', 0):.1f}</td></tr>"
        for l in labs)
    total_nt = sum(stats[l]['calls'].get('unique_scanned_nt', 0) for l in labs)
    return f"""<p class="mut">GLM_repeat 离线审阅包 · 全部链接为相对路径，无需联网</p>
<h1>从这里开始</h1>

<div class="warn"><b>这个包里全部是计算候选，没有任何一个经过人工确认。</b>
所有机器生成的标签都以 <code>AUTO_</code> 开头；所有需要你判断的字段都以
<code>HUMAN_</code> 开头，并且全部留空。</div>

<h2>这是什么</h2>
<p>用 gLM2_650M 基因组语言模型在 {total_nt/1e6:.2f} Mb 细菌/古菌与噬菌体序列上
扫描出的重复候选，以及对每个候选的结构解释尝试。共 {n_all:,} 个候选，
其中 {n_book} 个做成了可逐页翻阅的浏览册。</p>

<h2>这不是什么</h2>
<ul>
<li>报告的是<b>候选产出率</b>，<b>不是 FDR、不是精确率、不是生物学概率</b>。</li>
<li>explainer 给出模型 ≠ 该重复真实存在；probe 阴性 ≠ 该区域没有重复。</li>
<li>没有归类的候选只是"未解释"，<b>不是新发现的生物学系统</b>。</li>
</ul>

<h2>扫描规模</h2>
<table><thead><tr><th>运行</th><th>候选数</th><th>唯一扫描 Mb</th><th>候选/Mb</th></tr></thead>
<tbody>{rows}</tbody></table>

<h2>建议的审阅顺序</h2>
<ol>
<li>先读 <a href="report/ROUND2_REPORT_zh.md">report/ROUND2_REPORT_zh.md</a>
    —— 本轮做了什么、产出多少、有哪些局限。</li>
<li>再翻 <a href="index.html">浏览册索引</a>（{n_book} 个候选，可搜索、可排序）。
    每一页都有 detector 信号曲线、基因邻域图、全部重复假设与可下载序列。</li>
<li>把判断写进 <a href="review/browse_book_review.tsv">review/browse_book_review.tsv</a>
    的 <code>HUMAN_</code> 列。</li>
<li>需要全量数据时看 <a href="tables/all_candidates.tsv">tables/all_candidates.tsv</a>
    （全部 {n_all:,} 个候选）。</li>
</ol>

<div class="warn"><b>浏览册的两层必须分开解读。</b>
随机层 {sel_prov['selected_random']} 个（round-2 每条路线各
{sel_prov['n_random_per_route']} 个，种子 {sel_prov['seed']}）是<b>唯一</b>支持
无偏产出描述的部分；富集层 {sel_prov['selected_enriched_only']} 个是按规则挑选的导览，
<b>不能从中计算任何比例</b>。完整抽样名单见
<a href="manifests/selection_provenance.json">manifests/selection_provenance.json</a>。</div>

<h2>目录</h2>
<table><tbody>
<tr><th><a href="index.html">index.html</a></th><td>浏览册索引</td></tr>
<tr><th>candidates/</th><td>每个候选一页 HTML</td></tr>
<tr><th>figures/</th><td>每个候选的 detector 信号 SVG 与基因邻域 SVG</td></tr>
<tr><th>sequences/</th><td>每个候选的 FASTA</td></tr>
<tr><th>tables/</th><td>全量候选表</td></tr>
<tr><th>review/</th><td>待填写的人工审查表</td></tr>
<tr><th>report/</th><td>本轮报告与上游冻结报告</td></tr>
<tr><th>manifests/</th><td>数据选择、阈值、抽样的冻结清单</td></tr>
<tr><th>provenance/</th><td>代码校验和、环境、模型固定版本</td></tr>
<tr><th>MANIFEST.sha256</th><td>全部文件校验和</td></tr>
</tbody></table>

<h2>核验完整性</h2>
<p>解压后在包的根目录执行：</p>
<div class="seq">sha256sum -c MANIFEST.sha256</div>

<p class="mut">阈值、打分器与 explainer 均为冻结版本，本轮未做任何重新调参。</p>
"""


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package', type=Path, required=True)
    p.add_argument('--run', action='append', required=True, metavar='LABEL=DIR')
    p.add_argument('--selection-dir', type=Path, required=True)
    p.add_argument('--repo', type=Path, required=True)
    p.add_argument('--db', type=Path, required=True)
    a = p.parse_args()
    P = a.package
    stats = {}
    for spec in a.run:
        lab, _, d = spec.partition('=')
        stats[lab] = load(d)
    sel_prov = json.loads((a.selection_dir / 'selection_provenance.json').read_text())

    (P / 'report').mkdir(parents=True, exist_ok=True)
    (P / 'manifests').mkdir(parents=True, exist_ok=True)

    # ---- frozen upstream reports and manifests
    copied = []
    for src, dst in [
        (a.repo / 'results/pilot/calibration/thresholds/FROZEN.md', 'report/THRESHOLDS_FROZEN.md'),
        (a.repo / 'results/pilot/MERGING_AUDIT.md', 'report/MERGING_AUDIT.md'),
        (a.repo / 'results/explainer_v3/EXPLAINER_V3.md', 'report/EXPLAINER_V3.md'),
        (a.repo / 'results/probe_validation/PROBE_FROZEN.md', 'report/PROBE_FROZEN.md'),
        (a.repo / 'results/region_round1/FINDINGS.md', 'report/ROUND1_FINDINGS.md'),
        (a.repo / 'results/scorer/PEAK_CLUSTER_DESIGN.md', 'report/PEAK_CLUSTER_DESIGN.md'),
        (a.repo / 'results/pilot/calibration/thresholds/thresholds.json', 'manifests/thresholds.json'),
        (a.db / 'reports/ROUND1_REPORT.md', 'report/DB_ROUND1_REPORT.md'),
        (a.db / 'manifests/selection_config.json', 'manifests/db_selection_config.json'),
        (a.db / 'manifests/manifest.tsv', 'manifests/db_manifest_round1.tsv'),
        (a.db / 'manifests/candidates_round2.tsv', 'manifests/db_candidates_round2.tsv'),
        (a.db / 'manifests/download_ledger_round2.tsv', 'manifests/db_download_ledger_round2.tsv'),
        (a.db / 'manifests/exclusions.tsv', 'manifests/db_exclusions.tsv'),
        (a.selection_dir / 'selection_provenance.json', 'manifests/selection_provenance.json'),
        (a.selection_dir / 'browse_selection.tsv', 'manifests/browse_selection.tsv'),
    ]:
        if Path(src).exists():
            shutil.copy2(src, P / dst)
            copied.append(dst)
        else:
            print(f'  MISSING (not copied): {src}')

    # ---- per-run frozen provenance and summaries
    (P / 'provenance/runs').mkdir(parents=True, exist_ok=True)
    for spec in a.run:
        lab, _, d = spec.partition('=')
        for rel in ('calls/provenance.json', 'calls/summary.json',
                    'explain_full/summary.json', 'policy_attribution/summary.json',
                    'no_model_probe/summary.json', 'prepared/provenance.json',
                    'neighborhoods/provenance.json'):
            src = Path(d) / rel
            if src.exists():
                dst = P / 'provenance/runs' / lab / rel.replace('/', '__')
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)

    # ---- report and entry page
    n_all = sum(1 for _ in (P / 'tables/all_candidates.tsv').open()) - 1
    n_book = sum(1 for _ in (P / 'review/browse_book_review.tsv').open()) - 1
    rev = ann = None
    rp = P / 'review/assistant_review.summary.json'
    if rp.exists():
        rev = json.loads(rp.read_text())
        ann = {}
        for spec in a.run:
            lab, _, dd = spec.partition('=')
            t = list(csv.DictReader((Path(dd) / 'calls/trna_overlap.tsv').open(),
                                    delimiter=chr(9)))
            n_ov = sum(1 for r in t if int(r['n_tRNA_overlapping'] or 0))
            n_trna = len({x for r in t for x in (r['tRNA_coords'] or '').split(',') if x})
            ann[ROUTE_LABEL[lab]] = {'n_trna': n_trna,
                                     'cand_trna_pct': f'{100*n_ov/len(t):.1f}%'}
    report = build_report(stats, sel_prov, rev, ann)
    (P / 'report/ROUND2_REPORT_zh.md').write_text(report)

    css = (P / 'index.html').read_text().split('<style>')[1].split('</style>')[0]
    body = start_here(stats, sel_prov, n_book, n_all)
    (P / 'START_HERE.html').write_text(
        '<!doctype html><html lang="zh"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f'<title>GLM_repeat 审阅包</title><style>{css}</style></head>'
        f'<body>{body}</body></html>')
    (P / 'README_zh.md').write_text(
        '# GLM_repeat round-1 + round-2 离线审阅包\n\n'
        '浏览器打开 `START_HERE.html` 开始。全部链接为相对路径，无需联网。\n\n'
        '**包内全部是计算候选，没有任何一个经过人工确认。**\n'
        '机器标签以 `AUTO_` 开头；需要你判断的字段以 `HUMAN_` 开头，全部留空。\n\n'
        '核验完整性：在包的根目录执行 `sha256sum -c MANIFEST.sha256`。\n\n'
        '---\n\n' + report)

    # ---- checksums (MANIFEST.sha256 itself excluded)
    files = sorted(f for f in P.rglob('*') if f.is_file() and f.name != 'MANIFEST.sha256')
    with (P / 'MANIFEST.sha256').open('w') as fh:
        for f in files:
            fh.write(f'{sha256(f)}  {f.relative_to(P)}\n')
    print(f'package completed: {len(files)} files, {n_book} browse-book pages, '
          f'{n_all} candidates in the full table')
    print(f'copied {len(copied)} frozen artefacts')


if __name__ == '__main__':
    main()
