# tables/ — 原始明细

浏览册与报告用到的全部底层表。较大的原始明细表以 `.tsv.gz` 存放
（`gunzip` 或 `zcat` 即可读取，python `pandas.read_csv` 直接支持 `.gz`）。

| 文件 | 内容 |
|---|---|
| `yield_by_run.tsv` | 每轮产出率、基因组级 bootstrap 区间、单基因组分布 |
| `yield_by_genome.tsv` | **每个基因组**的候选数、唯一扫描 nt、/Mb、tRNA 占比 |
| `yield_method.json` | 统计方法与报告规则（为什么必须按基因组重抽） |
| `all_candidates.tsv` | 全部 4712 个候选及其浏览册归属 |
| `<run>_candidates.tsv` | 该轮候选表（坐标、分数、peak 数、来源窗口） |
| `<run>_peaks.tsv.gz` | 冻结 detector 实际调用的 peak（坐标、summit、均值、方向、来源窗口） |
| `<run>_models.tsv.gz` / `<run>_copies.tsv.gz` | explainer v0.3 的模型与逐条 copy 证据 |
| `<run>_neighbor_genes.tsv.gz` | 候选邻域内的预测 CDS（位置、朝向、距离、重叠） |
| `<run>_trna_overlap.tsv` | aragorn tRNA 注释与候选的相交情况 |
| `<run>_crispr_overlap.tsv` | pilercr CRISPR 阵列与候选的相交情况 |
| `review_digest.tsv` | 浏览册 115 个候选的证据摘要 |
| `seed_chance.tsv` | 12-mer 种子的保二核苷酸洗牌偶然性校准 |

`<run>` ∈ {region_round1, phage_round1, round2_host, round2_phage,
round3_host, round3_phage}。坐标一律为 0-based half-open。
