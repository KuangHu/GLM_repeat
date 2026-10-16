# 第一轮区域扫描：CDS 筛选、审计对照与基因邻域

版本：region-scan v0.1，2026-10-16。此版本增加可选扫描范围选择，不修改冻结的 detector 或 v0.3 repeat explainer。CPU 代码及集成测试已在本地运行；GPU 推理、Slurm 调度和大规模吞吐需要在 cluster 验证。

## 1. 设计与适用范围

流程：原始 FASTA → CDS 预测 → 可审查的扫描窗口选择 → 原来的 gLM2 scorer → 严格检查得分完整性 → 候选及原始成员 → v0.3 repeat hypotheses → 基因邻域、GenBank 与图。

输入面向细菌、古菌、噬菌体 contig。第一轮建议每批约 5–10 Mb；各 CPU 阶段会把这一批序列和表载入内存。一个 FASTA 内的 contig ID 必须全局唯一。对不同样品都叫 `contig_1` 的情况，应先统一重命名 FASTA、GFF 和标签。无需蛋白 anchor 或已知 repeat motif。基因功能只参与后续解释。

**这里的“非编码”是未被所选 CDS 策略排除的区域，不是实验确认的非编码区。** Prodigal 的高 confidence 也可能误报 ORF。不存在一个经过本项目验证的 confidence 值，可以保证不切断长 ncRNA/repeat。默认参数是明确记录的启发式起点。

### 默认参数

| 参数 | 默认值 | 含义 |
|---|---|---|
| caller / mode | pyrodigal / meta | 可显式切换 pyrodigal-gv；single 只用于同一个基因组的 contig |
| policy | conservative | 下述完整、长且高分的 CDS 才可以缩小目标范围 |
| min-cds | 300 nt | CDS 各部分总长度门槛，不是普适生物学界限 |
| cds-confidence | 99 | 需要 GFF 的 `conf>=99`、`cscore>0`、`partial=00` 同时成立 |
| min-noncoding | 40 nt | 对可排除 CDS 的补集应用，先判断长度，后加上下文 |
| boundary-pad | 75 nt | 为保留的目标区域增加两侧边界缓冲 |
| window / stride | 3000 / 1500 nt | 使用与原 pipeline 相同的原始网格窗口 |
| min-segment | 40 nt | 非 ACGT 分割后的更短片段不评分，会计入报告 |
| audit-fraction | 0.02 | 按 contig 数抽完整 contig 对照，向上取整，固定种子 |
| seed | 20261016 | 仅用于对照抽样；GPU scorer 仍使用冻结的 20261014 |

`conservative` 保留短 CDS、partial CDS、缺少分数的 CDS 以及低分 CDS；它们不会切碎长的保留区域。若不确定 CDS 与另一个可排除 CDS 重叠，不确定 CDS 的区域仍会被恢复为扫描目标。joined CDS 只排除实际 CDS parts，不排除 parts 之间的间隙。

另外有两个明确的对照策略：`--policy none` 保留所有可评分序列，忽略非编码长度阈值；`--policy all` 将所有预测 CDS 用于排除，适合研究过度过滤的影响，默认不使用。外部 GFF 若无 `conf/cscore/partial`，conservative 会保留这些 CDS；不要为增加过滤量自行补造置信度。

所有原始 DNA 坐标和碱基保留。基因预测把非 ACGT 字符等长转换成 N，并禁止跨这些位置预测 CDS；评分窗口在非 ACGT 处断开。不同片段不拼接。CDS 不替换成 N，也不从模型上下文中删除。

目标区域和窗口是两层：任何与目标区域相交的原始窗口都可评分，模型因此仍可能看到 CDS；窗口内产生的 CDS/context 候选会输出并标记。3 kb 上下文会削弱节省算量的比例，这是设计上的取舍。

当前按线性 contig 处理。环状序列在原点两侧的基因、repeat 和邻域不会被跨原点重建，邻域边缘有标记。不能用本版本声称完整覆盖环状原点的 repeat。只有单个大染色体时，默认至少抽一个完整 contig 可能等于全扫描；请先查看 summary，而不是假定“2% 对照就是 2% 碱基”。

## 2. 文件与接口

| 文件 | 职责 |
|---|---|
| `scripts/predict_cds.py` | 调用 Pyrodigal/Pyrodigal-gv，输出 genes.gff3、proteins.faa、contigs.tsv 和版本/输入哈希 |
| `scripts/prepare_region_scan.py` | 生成目标、CDS 决策、所有窗口决策、所选窗口及实际算量 |
| `scripts/benchmark_region_selection.py` | 固定预测结果的 0/40/60/100 nt 长度敏感性实验；标签只用于事后覆盖统计 |
| `scripts/score_region_chunk.py` | 安全分片后调用原来的 timing_score.py；空分片不加载模型 |
| `slurm/region_score_lrc.sbatch` | 可配置路径的 Lawrencium GPU array 示例 |
| `scripts/call_region_candidates.py` | 严格得分检查、原 detector 调用、候选成员及对照组反事实比较 |
| `scripts/repeat_neighborhoods.py` | 按候选提取邻域、关联 CDS/repeat 证据、导出 GenBank 与 pyGenomeViz 图 |
| `scripts/region_io.py` | 统一坐标、GFF3、区间并集、输入检查和 provenance |
| `tests/test_region_scan.py` | 误报 ORF、区间、审计、缺失得分、整条处理链等行为测试 |

所有 TSV 坐标为 **0-based half-open** `[start0,end0)`；GFF3 自动从 1-based inclusive 转换。`candidate_id=contig:start0-end0` 与 v0.3 explainer 一致。图中横轴是相对于提取区域起点的距离，图注给出原始 contig 区间；GenBank 使用局部坐标，`neighborhoods.tsv` 给出局部到原始坐标的映射。

## 3. Cluster 安装与路径

推荐将压缩包解压到一个新目录，或将包内文件同步到现有仓库。旧的模型目录和结果保留。源码包不带模型或真实序列。以下命令都从代码目录执行；用同一个 shell 保留变量。

```bash
cd /path/to/GLM_repeat_region_scan_v0.1
export PROJECT="$PWD"

# CPU 预处理/邻域独立环境，不升级已经验证的 GPU 环境。
python3 -m venv "$PROJECT/.venv-region"
export REGION_PYTHON="$PROJECT/.venv-region/bin/python"
"$REGION_PYTHON" -m pip install -r requirements-region-scan.txt

# 使用现有、已验证的 GPU 环境和已审查、固定 revision 的模型。
export GLMREP_PYTHON="$HOME/.conda/envs/glm_repeat/bin/python"
export GLMREP_SNAPSHOT=/global/home/users/kh36969/tools/GLM_repeat/models/glm2_650M

export GENOME_FASTA=/absolute/path/to/input.fa
export RUN=/absolute/path/to/new_region_round1
mkdir -p "$RUN" logs

"$REGION_PYTHON" -m pytest -q tests/test_region_scan.py tests/test_explainer_v3.py \
  tests/test_windows.py tests/test_peak_cluster.py
```

FASTA 用未压缩 `.fa`，以兼容后面的 explainer；基因预测和窗口准备本身也支持 `.gz`。上述路径需要替换。Slurm 文件保留原项目的账户/分区/A40 设置，若集群配置不同需调整。CPU 大任务也应按站点规范提交 CPU 作业。

## 4. 扫描前：预测与选择

混合 contig / 初次探索：

```bash
"$REGION_PYTHON" scripts/predict_cds.py \
  --fasta "$GENOME_FASTA" --out "$RUN/cds" --caller pyrodigal --mode meta

"$REGION_PYTHON" scripts/prepare_region_scan.py \
  --fasta "$GENOME_FASTA" --gff "$RUN/cds/genes.gff3" \
  --out "$RUN/prepared" --policy conservative \
  --min-noncoding 40 --boundary-pad 75 --audit-fraction 0.02

cat "$RUN/prepared/summary.json"
```

单个细菌/古菌基因组的充足序列可使用 `--caller pyrodigal --mode single --translation-table 11`，必须确认同一个文件的所有 contig 来自同一个基因组。single 训练不足会报错，不悄悄改成另一种模式。

噬菌体可在新的输出目录重复 `--caller pyrodigal-gv --mode meta`，比较 CDS 与窗口差异。meta 自动选择模型及遗传密码，结果保存在 `cds/contigs.tsv`；此时不能用 `--translation-table` 假装强制另一种密码。

主要输出：

- `cds_decisions.tsv`：每个 CDS 是否满足当前排除策略；实际窗口仍可能覆盖它。
- `regions.tsv`：保留补集、太短补集、不确定 CDS 恢复区域。
- `targets.tsv`：包含边界缓冲的最终目标区间并集。
- `windows_all.tsv`：所有原始窗口，包括跳过原因和对照标记。
- `windows.tsv` / `windows.fa` / `records.txt`：真正送给 scorer 的输入。
- `summary.json`：完整网格、策略网格、含对照网格三套窗口数和总 nt；减幅是输入长度代理，不是实测 GPU 时长。

`--audit-contigs /path/to/ids.txt` 可额外指定完整 contig 对照，每行一个 ID。指定对照是加法；若只要显式指定的对照，可加 `--audit-fraction 0`。不依赖标签挑选对照，避免只审计熟悉的正例。

若有开发集位点/单元标签，可先运行 CPU 敏感性比较：

```bash
"$REGION_PYTHON" scripts/benchmark_region_selection.py \
  --fasta "$GENOME_FASTA" --gff "$RUN/cds/genes.gff3" \
  --loci /path/to/loci.tsv --units /path/to/repeat_units.tsv \
  --out "$RUN/length_comparison" --lengths 0 40 60 100
```

标签至少含 `contig start0 end0` 三个制表符分隔列；`--units` 可省略。此脚本自动关闭随机对照来测量纯策略的覆盖率；不会根据结果挑选阈值。参数比较使用开发集，最终泛化评估仍需独立基因组。`min-noncoding=0` 仍使用 CDS 策略，不等于 `policy=none`。

## 5. GPU 扫描与候选调用

先核对 `summary.json` 的 `n_selected_windows` 和实际输入规模。若为 0，不提交 GPU 作业；这是过滤结果，不代表没有 repeat。

```bash
# 根据窗口数和站点配额设置 array 大小；空分片会安全退出。
N=16
sbatch --array=1-"$N" slurm/region_score_lrc.sbatch \
  "$RUN/prepared" "$RUN/scores" "$N" random_partition
```

该 wrapper 使用原 scorer 的双方向、20 groups × 2 repeats、seed 20261014、attention-mask-policy=always，匹配已有 pilot 标定。保留原始 scorer 实现；不在这轮改变模型上下文。单个 shard 失败时，保留其失败记录并使用新的输出位置重跑，汇总时只能纳入每个窗口/方向的一份完成结果。不要把缺少的 shard 当作阴性。

所有 array 任务完成后：

```bash
"$REGION_PYTHON" scripts/call_region_candidates.py \
  --prepared "$RUN/prepared" --scores "$RUN/scores" \
  --threshold-file results/pilot/calibration/thresholds/thresholds.json \
  --threshold-key bg100_per_Mb --out "$RUN/calls"
```

调用前必须同时满足：每个窗口有正反两个得分；长度和序列哈希一致；全部值有限；坐标约定正确；模型 revision、mode、groups/repeats/seed/attention policy 与阈值文件一致。缺失、重复、混合模型都会报错。

当前随包阈值来自旧的全基因组 pilot。**`bg100_per_Mb` 是历史标定名称，不保证筛选后每 Mb 仍为 100 个候选，也不是 FDR。** 程序同时报告每 scanned Mb 的全部候选密度和每 target Mb 的关联候选密度。若改变模型或评分模式，请用同一选择策略的独立背景重新标定；`--threshold 0.8` 只适合显式标记为探索性的运行，不带密度保证。

产物：

- `candidates.tsv`：所有候选，分类为 target_associated、context_only、audit_only。
- `members.tsv`：原窗口、方向、原始 cluster 区间和分数。
- `audit_full.tsv` / `audit_policy.tsv`：对完整对照 contig，用同一份得分分别运行全窗口和仅策略窗口两种调用，统计哪些全扫描候选在策略结果里没有任何重叠。
- `summary.json`：完整性、算量、密度和审计计数。

原 detector 保留 ≤500 nt peak gap 规则和跨窗口重叠合并。这里的合并仅表示候选包络，不代表确定的单个 repeat locus。链式重叠标记为 `transitive_chain`，原始成员保留供检查；后续不因 ±10 kb 邻域重叠继续合并。对照中的 overlap recovery 也不等于逐个 repeat copy 的召回率。

## 6. Repeat 解释与 gene architecture

```bash
"$REGION_PYTHON" scripts/explain_repeats_v3.py \
  --candidates "$RUN/calls/candidates.tsv" --genomes "$GENOME_FASTA" \
  --out "$RUN/explainer"

"$REGION_PYTHON" scripts/repeat_neighborhoods.py \
  --fasta "$GENOME_FASTA" --gff "$RUN/cds/genes.gff3" \
  --proteins "$RUN/cds/proteins.faa" \
  --candidates "$RUN/calls/candidates.tsv" --members "$RUN/calls/members.tsv" \
  --copies "$RUN/explainer/copies.tsv" --models "$RUN/explainer/models.tsv" \
  --flank 10000 --nearest-genes 5 --plot --out "$RUN/neighborhoods"
```

提取范围同时覆盖候选 ±10 kb 与左右最近各 5 个基因；稀疏区域可能大于 20 kb。包含候选重叠 CDS。其他基因若仅部分落在提取窗口内，会标记 clipped，并省略其完整 protein translation，避免把完整蛋白附到截短 CDS。序列边界会截断邻域；没有跨 contig 或环状原点拼接。

`--models models.tsv` 只限制图和导出的 repeat 图层使用展示代表模型，不删除 explainer 的替代假设；省略它可以显示全部模型。相同坐标只画一次，原始 copy/model 表继续保留。这里展示的是 core/periodic copies；v0.3 的可选延伸仍在 extensions.tsv 中，未混成额外独立 copies。

查看 `neighborhoods/index.html`；`plots/` 下是 SVG/PNG；`genes.tsv` 含每个邻域基因的方向、左/右/重叠、距离、product 与 GFF attributes；`genbank/` 可作为后续 clinker 输入。

CDS caller 只能预测 ORF，默认 product 是 hypothetical protein。若需要功能架构，先用 Bakta（细菌）或 Pharokka（噬菌体）取得全 contig 的功能 GFF3/蛋白，再把对应文件传给 **邻域步骤**。不必因此改变前面的扫描选择。GFF 中的 CDS ID 与 protein FASTA ID 必须匹配；程序报告每个邻域实际附带 translation 的数量。缺少注释或 translation 不应解读为没有功能或同源性。

可用 clinker 比较挑选后的同类候选邻域（单独安装，不是此包的必要依赖）：

```bash
# 将有比较意义的少量 GenBank 文件放入 selected_neighborhoods；避免一次比较所有候选。
clinker selected_neighborhoods/*.gbk -p neighborhood_comparison.html
```

pyGenomeViz 已集成用于自己的候选/重复/CDS 图层；clinker 用于蛋白同源及基因排列比较。没有在此版本自动推断 operon、功能模块或 repeat 的转录方向。邻近 Cas 等蛋白可以成为后验线索，但不作为命中的前置条件。

## 7. 已运行的验证及限制

本地安装并调用了 Pyrodigal 3.7.1、Pyrodigal-gv 0.3.2、pyGenomeViz 1.6.1。旧测试 70 项与新测试 25 项共 **95 项通过**（CPU 集成 70 项；需要 PyTorch 的原测试另在已有 PyTorch 环境运行 25 项）。包括实际基因预测 smoke test，以及模拟得分→候选→v0.3 explainer→GenBank 的整条连接测试。Slurm 文件通过 Bash 语法检查。

使用已有 8 条公开控制序列（5,172,072 nt）的 CPU 预实验：

| 配置 | 完整网格窗口 | 策略窗口 | 输入窗口总 nt 减幅 | 完整覆盖检查区域 | 完整覆盖标注单元 |
|---|---:|---:|---:|---:|---:|
| min-noncoding 0 | 3441 | 3279 | 4.71% | 15/15 | 230/230 |
| min-noncoding 40 | 3441 | 2884 | 16.20% | 15/15 | 230/230 |
| min-noncoding 60 | 3441 | 2636 | 23.40% | 15/15 | 230/230 |
| min-noncoding 100 | 3441 | 2179 | 36.68% | 15/15 | 230/230 |

这些行固定 conservative、boundary-pad 75、关闭审计对照。默认 40 nt 加抽中的完整 contig 对照后为 2903 个窗口、减幅 **15.65%**。Pyrodigal 与 Pyrodigal-gv 分别产生 5540 / 5527 个 CDS，但在这批数据的默认配置下最终所选窗口数相同；不据此宣布哪一个 caller 更好。

这 15 个区域包含历史 CRISPR/ART 控制、VIPR 开发位点及一个 off-label repeat pair。现有单元标签仅覆盖其中一部分，不是独立 VIPR 验证。100 nt 参数仍覆盖 60 nt VIPR 位点，是因为所在保留区域/窗口足够大；不能据此认为所有短 repeat 都能通过 100 nt 过滤。

**上述统计是几何覆盖和计算量代理，不是新 GPU 检测结果，不是运行时间加速比，不是全体 repeat 的召回率。** CPU 测试也不保证新输入分布的 CDS 假阳性率。第一轮应保持对照并检查 audit_full_without_policy_overlap，必要时进一步用人工和 repeat copy 证据解释遗漏。

已有 v0.3 限制仍存在：motif 分支依赖 exact 12-mer seed；只处理 ungapped substitutions；period 搜索范围为 2–300；得分和 support 不是经过校准的生物学概率。不得把当前 pipeline 描述成完整的 general repeat detector。

## 8. 依据与复现

- Pyrodigal API / 坐标 / confidence：[GeneFinder](https://pyrodigal.readthedocs.io/en/stable/api/gene_finder.html)、[Genes](https://pyrodigal.readthedocs.io/en/stable/api/genes.html)。
- 病毒扩展模型：[Pyrodigal-gv](https://github.com/althonos/pyrodigal-gv)。
- 邻域图接口：[pyGenomeViz](https://moshi4.github.io/pyGenomeViz/getting_started/)。
- 蛋白同源基因簇比较：[clinker](https://github.com/gamcil/clinker)。
- 后验功能注释：[Bakta](https://github.com/oschwengers/bakta)、[Pharokka](https://github.com/gbouras13/pharokka)。

各阶段写入 provenance.json，保存参数、输入文件哈希、Python 源码哈希。基因预测另保存 caller 版本。新输出目录要求为空，以免混合旧结果。源码包附 MANIFEST.sha256，可在解压目录运行 `sha256sum -c MANIFEST.sha256` 验证。
