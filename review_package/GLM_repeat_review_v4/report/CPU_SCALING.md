# CPU 分配实测（决定本轮与后续的计费成本）

前三轮的扫描作业都申请 `-c 4`，但 `sacct` 显示 CPU 效率恒为 **24.9%**
（round3_host 与 round3_phage 完全一致），即只有 1 个核在工作。
由于 Slurm 在本集群按 **CPU 数**计费（`es1` 未设 `TRESBillingWeights`，
`billing=4` 就等于 `cpu=4`，A40 本身不计入 billing），多申请的核是纯浪费。

## 实测

同一个 prepared 窗口集、同一个 chunk（29 个窗口 / 58 条 track），只改 `-c`：

| `-c` | 状态 | 耗时 | 每窗口 | 相对 |
|---|---|---:|---:|---:|
| 1 | **无法运行** | — | — | QOS `es_normal` 规定 `MinTRES cpu=2`，提交后永远 PENDING（`QOSMinCpuNotSatisfied`） |
| 2 | COMPLETED | 623 s | 21.5 s | 基准 |
| 4 | COMPLETED | 622 s | 21.4 s | −0.2% |

**4 个核与 2 个核的运行时间差 0.2%，在噪声范围内。** 工作负载受 GPU 限制，
CPU 侧基本单线程，因此 `-c 2` 是这个 QOS 下能用的最省配置。

## 影响

| | `-c 4`（旧） | `-c 2`（新） |
|---|---:|---:|
| db10（36.73 Mb, 24,305 窗口, 145 job-h） | 580 billing-h | **290 billing-h** |
| 整库剩余（366 Mb） | ~5,860 billing-h | **~2,930 billing-h** |

前三轮已按 `-c 4` 运行，合计 796.3 billing-h；其中约一半是空转的核。
本轮起改用 `-c 2`。**这只改变资源申请，不改变任何科学参数**：
模型、阈值、打分器、explainer、窗口化全部保持冻结。

## 注意

集群门户报告的 SU 与 Slurm 的 billing-hour **不是同一个单位**
（2026-10-01 起全项目 57,398 billing-h 对应门户 25,330.28 SU），
换算在 Slurm 之外、按分区计费，无法从集群内部推导。
上表一律是 billing-hour，不是 SU。
