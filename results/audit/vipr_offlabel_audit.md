# VIPR off-label audit (SUSP1, window w0003) — NOT benchmark ground truth

Status: sequence-supported repeat candidate; vrRNA identity UNCONFIRMED.
This file is kept apart from data/bench labels. Benchmark metrics continue to use
only the original label SUSP1_vrRNA_35TR_2 [29135,29229) minus strand.

Candidate: NC_028808.2:[27682,27906) (0-based half-open). It begins at the
window's left edge, so it may extend further upstream; that flank was not scored.

Evidence (sequence, independent of gLM2):
- 39-nt exact direct repeat AGTCTCATTTTAATGTAGCTGAGACTAACTACCCATACC at 27714 and
  27857 (2 copies, both + strand; no other copy in the genome). It contains an
  internal inverted repeat (GTCTCA ... TGAGAC).
- 16 YNNGG pentamer matches on the minus strand within 240 nt; genome-wide
  expectation ~3.7 per 240 nt. The labelled 94-nt vrRNA has 12.

Evidence (model, exact track from timing6_exact, orientation-averaged):
- mean 0.922 nats above uniform over [0,224) vs 0.437 for the labelled vrRNA and
  0.279 for the whole window; the six highest 64-nt blocks of w0003 all lie here.

Not done / needed for confirmation: comparison with the SUSP1 VIPR system
annotation or author-deposited RNA set (other vrRNAs), RNA-seq evidence,
scoring of the upstream flank.
