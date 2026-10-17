# GLM_repeat round-1 + round-2 离线审阅包

浏览器打开 `START_HERE.html` 开始。全部链接为相对路径，无需联网。

**包内全部是计算候选，没有任何一个经过人工确认。**
机器标签以 `AUTO_` 开头；需要你判断的字段以 `HUMAN_` 开头，全部留空。

核验完整性：在包的根目录执行 `sha256sum -c MANIFEST.sha256`。

---

# GLM_repeat round-1 + round-2 扫描报告（供人工审阅）

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

round-2 新增 9.92 Mb（628 个候选），四次运行合计
19.09 Mb、1197 个候选。

| 运行 | 候选数 | 唯一扫描 nt | 候选/Mb | 目标区 候选/Mb | CDS 策略丢失包络 |
|---|---|---|---|---|---|
| 宿主 round-1 | 170 | 4,616,499 | 36.8 | 117.4 | 6 |
| 噬菌体 round-1 | 399 | 4,551,396 | 87.7 | 236.0 | 17 |
| 宿主 round-2 | 156 | 4,819,813 | 32.4 | 138.5 | 0 |
| 噬菌体 round-2 | 472 | 5,100,950 | 92.5 | 245.0 | 30 |

两次独立的宿主面板给出 36.8 与
32.4 候选/Mb，两次噬菌体面板给出
87.7 与
92.5 候选/Mb。
噬菌体路线的产出率约为宿主路线的 2.7 倍，这在两轮中都成立。

**这不是一个生物学结论。** 噬菌体基因组更短、编码密度更高、窗口覆盖方式不同，
而阈值是在宿主背景上标定的；产出率差异至少有一部分来自这些非生物学因素。

## 2. explainer（v0.3，冻结）

| 运行 | candidate_associated | flank_only | unresolved | 主模型 | 假设总数 |
|---|---|---|---|---|---|
| 宿主 round-1 | 74 | 18 | 78 | 92 | 459 |
| 噬菌体 round-1 | 215 | 36 | 148 | 251 | 1385 |
| 宿主 round-2 | 58 | 9 | 89 | 67 | 348 |
| 噬菌体 round-2 | 275 | 36 | 161 | 311 | 1507 |

`unresolved` 在两条路线上都占三到六成。这些候选的 detector 信号存在，但
ungapped explainer 没能给出模型。第 4 节的 probe 用于区分其中哪些仍带有
indel 型相似性。

## 3. CDS 策略的影响

| 运行 | full copies | policy copies | 精确保留 | 丢失 | 新增 | 受影响候选 |
|---|---|---|---|---|---|---|
| 宿主 round-1 | 2,079 | 1,807 | 1,792 | 287 | 15 | 8 |
| 噬菌体 round-1 | 5,837 | 5,477 | 5,445 | 392 | 32 | 18 |
| 宿主 round-2 | 2,641 | 2,641 | 2,641 | 0 | 0 | 0 |
| 噬菌体 round-2 | 6,201 | 5,482 | 5,443 | 758 | 39 | 34 |

**CDS 策略的影响是基因组依赖的，不是系统性的。** 宿主 round-2 丢失 0 个包络、
0 条 copy 记录；噬菌体 round-2 丢失 30 个包络、758 条 copy 记录。
归因表（`tables/` 与各运行的 `policy_attribution/attribution.tsv`）把差异拆成
两个独立的轴：`upstream_trigger`（哪一级流水线发生了变化）与
`representation_outcome`（这条 copy 是否仍以某种形式出现）。
只有 `no_corresponding_copy` 表示证据在 policy 输出中确实缺失，
**即便如此，它也是计算证据的差异，不是经过验证的生物学漏检。**

## 4. 无模型候选的 indel 感知 probe（仅诊断）

| 运行 | 无模型候选 | no_repeated_seed | seeds_without_extension | ungapped_similarity | gapped（explainer 无法表示） | 阳性对照灵敏度 |
|---|---|---|---|---|---|---|
| 宿主 round-1 | 78 | 6 | 50 | 21 | 1 | 18/20 |
| 噬菌体 round-1 | 148 | 26 | 89 | 31 | 2 | 17/20 |
| 宿主 round-2 | 89 | 26 | 49 | 14 | 0 | 20/20 |
| 噬菌体 round-2 | 161 | 28 | 91 | 41 | 1 | 18/20 |

probe 不参与候选判定，也不回写任何冻结输出。阳性对照是强制的：
它在 explainer **已经**建模的候选上运行，用来证明这个 probe 本身有灵敏度。

分层验证（`report/PROBE_FROZEN.md`）已经说明，该 probe 在 8–12 nt 单元上的
灵敏度只有 0.20–0.22，并且在典型区域长度上有 13–26% 的假阳性率。因此：

- **probe 阳性 ≠ 可信重复**；典型长度下相当一部分是偶然相似。
- **probe 阴性 ≠ 假阳性候选**；短单元重复大多会被漏掉。

## 5. 浏览册的抽样方式

浏览册共 69 个候选，分两层，**必须分开解读**：

- **随机层（40 个）**：从 round-2 两条路线各随机抽
  20 个，种子 20261016，完整抽样名单保存在
  `manifests/selection_provenance.json`。**只有这一层支持无偏的产出描述。**
  被随机抽中的候选即使同时被富集规则选中，仍然留在随机层的分母里。
- **富集层（29 个）**：按显式命名的规则挑选
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


## 8. 对浏览册 69 个候选的审阅（助手判断，非人工确认）

判断写在 `review/assistant_review.tsv` 的 `AUTO_REVIEW_*` 列。
**`HUMAN_*` 列仍然全部留空** —— 机器判断不应占据留给人工阅读的字段。
每条判断都由明确规则产生，并附带触发它的证据（`AUTO_REVIEW_reasoning`），
因此任何一条都可以被逐项质疑。`cannot_judge` 是一等结论，不会为了好看而回避。

**随机层（n=40）**

- 重复证据：no_repeat_found 22, cannot_judge 7, suspected_repeat 6, clear_repeat 5
- 解释完整性：explainer_wrong_model_class 14, nothing_to_explain 12, boundary_or_phase_problem 5, partial 5, model_adequate 4
- 后续价值：known_common_structure 18, not_worth_followup 10, function_unknown 8, worth_followup 4
- 与注释 tRNA 相交：15；位于 CRISPR 阵列内：0；携带 CRISPR repeat 但不在阵列内：1

**富集层（n=29）**

- 重复证据：clear_repeat 12, suspected_repeat 8, no_repeat_found 7, cannot_judge 2
- 解释完整性：model_adequate 7, partial 6, structure_unsupported 5, nothing_to_explain 4, explainer_wrong_model_class 4, model_is_fragment_of_known_array 2, no_model_but_signal_present 1
- 后续价值：known_common_structure 12, worth_followup 9, not_worth_followup 4, function_unknown 4
- 与注释 tRNA 相交：6；位于 CRISPR 阵列内：2；携带 CRISPR repeat 但不在阵列内：0

**两层不可合并。** 只有随机层描述"审阅者会遇到什么"；富集层是按规则挑选的导览。

### 8.1 本轮审阅得到的主要结论

**(1) 宿主路线上约一半的候选是 tRNA 基因。**

| 运行 | 注释 tRNA 数 | 候选中与 tRNA 相交 |
|---|---:|---:|
| 宿主 round-1 | 124 | 52.4% |
| 噬菌体 round-1 | 241 | 17.0% |
| 宿主 round-2 | 102 | 52.6% |
| 噬菌体 round-2 | 267 | 11.9% |

两轮宿主面板分别为 52.4% 与 52.6%，两轮噬菌体为 17.0% 与 11.9%，各自高度一致。
这意味着在当前阈值下，**宿主路线的 detector 在很大程度上是一个 tRNA / 结构性 ncRNA 检测器**。
tRNA 的 masked predictability 高，是因为它们受到强烈的结构与序列约束，
**这与"存在序列重复"是两回事**。

**(2) `unresolved` 主要不是失败，而是模型类别不匹配。**
宿主 round-2 的 89 个 `no_model` 候选中有 60 个是 tRNA。
ungapped repeat explainer 无法表示结构性 ncRNA 基因，这是类别错配，不是漏解释。
报告第 2 节把 `unresolved` 当作"未解释"统计，这个口径需要在下一轮修正。

**(3) CRISPR 阵列被完整召回，并额外提示 4 个未被专用工具调用的位点。**
pilercr 在宿主 round-2 注释出 9 个 CRISPR 阵列，**9 个全部被候选覆盖**。
另有 4 个 *Aquifex aeolicus* 候选携带该基因组的 CRISPR repeat consensus
（`GTTCCTAATGTACCGTGTGGAGTTGAAACC`）却落在所有 pilercr 阵列之外，
其中一个（`NC_000918.1:263146-263290`）出现在随机层。
这类位点可能是退化或残缺的阵列，是本轮最值得单独核查的结果。
**这不是"发现了新的生物学系统"**：它们是专用工具未调用、而本流程标记出来的已知类型位点。

**(4) 12-nt 种子的偶然性必须逐例校准。**
对每个候选区域做保二核苷酸组成的洗牌（300 次）后比较重复 12-mer 数量：
`candidate_associated` 的模型几乎全部显著（p<0.05），
而 `no_model` 候选中只有 1/30 显著 —— 它们确实没有超出偶然的短词重复。
但 `flank_only` 的 `repeat_pair` 只有一半显著，且其 copies 全部落在候选之外，
**这类模型不构成关于该候选本身的证据**，一律记为 `cannot_judge`。

**(5) 低复杂度必须同时看 consensus，不能只看包络。**
包络常有数千 nt，会把低复杂度的 consensus 稀释到看不见。
例如 `NC_048045.1:37076-40240` 的 consensus 是 `GGCGGCGGCGGC`（熵 0.92），
而整个包络的熵是 1.86。本轮已把 consensus 熵加入低复杂度判定。

### 8.2 这轮审阅没有做的事

- 没有重新调参：detector、阈值、explainer 全部保持冻结。
- 没有把 tRNA / CRISPR 注释反馈进候选判定，它们只用于解释与展示。
- 没有宣称任何候选是经过确认的生物学发现。
