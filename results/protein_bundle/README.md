# GLM_repeat protein bundle v1

七个 panel 的 pyrodigal 预测蛋白，**ID 已全局唯一化**，附坐标与候选邻域索引。
供 domain/function 注释与蛋白聚类使用。

## 为什么必须用这份，而不是直接 cat 原始文件

pyrodigal 按 run 编号，**每个 panel 的 CDS ID 都从 `cds00000001` 重新开始**。
七个 `proteins.faa` 直接拼接会让 **24,119 个原始 ID 同时出现在多个 panel**，
把互不相关的蛋白合并到同一个标识符下 —— 大多数聚类工具会静默保留第一条或去重，
**不会报错**。本 bundle 的 ID 为 `<panel>|<original_cds_id>`，122,704 条全部唯一。

## 内容

| 路径 | 说明 |
|---|---|
| `proteins/all_panels.faa` | 全部 122,704 条蛋白，ID 已唯一化 |
| `proteins/<panel>.faa` | 按 panel 分开的同一批序列 |
| `protein_index.tsv` | 每条蛋白 → panel、contig、坐标、链、长度、partial、pyrodigal 置信度，以及候选邻域关联 |
| `provenance/` | 七个 panel 的 cds / prepared / calls / neighborhoods 运行来源记录 |
| `bundle_manifest.json` | 计数、唯一性断言、口径说明 |

FASTA header 形如：
```
>r3_host|cds00012345 panel=r3_host contig=NZ_AP024685.1 start0=12345 end0=13456 strand=+ len_aa=370
```

## protein_index.tsv 的候选关联列

- `in_candidate_neighborhood` / `n_candidate_neighborhoods` / `candidate_ids` /
  `relations` / `max_overlap_nt_with_candidate`
  来自冻结的邻域导出（候选 ±flank）。
- 122,704 条中有 **88,225 条**落在至少一个候选邻域内。
- **"靠近候选"只是位置事实，不是关联证据。** 邻域窗口较宽，
  因此多数蛋白都会落入某个邻域；若要问"某个蛋白家族是否在候选附近富集"，
  必须与一个合适的零模型比较，不能直接用这个比例。

## 两条必须随结论保留的限制

1. **没有任何功能注释。** 这些是 pyrodigal 从头预测，源 GFF 中每个
   `product` 都是 "hypothetical protein"。所有 domain / function 判断
   只能来自你即将运行的分析，不能引用本流程的上游。
2. **坐标为 0-based half-open**，contig 名见 header 与索引表。

## 原始文件位置（未改动）

```
/global/scratch/users/kh36969/GLM_repeat/<dir>/cds/proteins.faa
```
`<dir>` ∈ region_round1, phage_round1, round2_host, round2_phage,
round3_host, round3_phage, db10_phage
