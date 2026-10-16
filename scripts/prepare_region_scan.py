#!/usr/bin/env python3
"""Select original scoring windows using CDS evidence; keep audit contigs whole.

No DNA masking, stitching, motif gate or label-dependent selection is performed.
The conservative policy is an explicit heuristic, NOT verified coding truth.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

from make_windows import valid_segments, window_starts
from region_io import (SpanIndex, complement, fasta, interval, new_dir, provenance,
                       read_gff, table, union, write_table)


def excludes_gene(g, policy, min_cds, confidence):
    if policy == 'none':
        return False
    if policy == 'all':
        return True
    attrs = g['attrs']
    # Missing confidence/partial metadata is uncertainty, never a zero-gene call.
    try:
        conf, coding = float(attrs['conf']), float(attrs['cscore'])
    except (KeyError, ValueError):
        return False
    return (attrs.get('partial') == '00' and g['length_nt'] >= min_cds
            and math.isfinite(conf) and conf >= confidence and math.isfinite(coding)
            and coding > 0)


def regions_for(seq, genes, policy, min_noncoding, boundary, min_cds, confidence):
    excluded, uncertain = [], []
    decisions = []
    for g in genes:
        skip = excludes_gene(g, policy, min_cds, confidence)
        decisions.append({**g, 'exclude_by_policy': skip,
                          'reason': 'CDS_policy' if skip else 'uncertain_or_policy_none'})
        (excluded if skip else uncertain).extend(g['parts'])
    excl_index, rescue_index = SpanIndex(excluded), SpanIndex(uncertain)
    rows, targets = [], []
    for ss, ee in valid_segments(seq):
        # Length applies to complement of excludable CDSs, before padding.
        # Weak/partial CDSs thus cannot split the long retained complement.
        for s, e in complement(excl_index.spans, ss, ee):
            keep = policy == 'none' or e-s >= min_noncoding
            rows.append(dict(start0=s, end0=e, reason='eligible_complement' if keep
                             else 'short_complement', eligible=keep))
            if keep:
                targets.append((max(ss, s-boundary), min(ee, e+boundary)))
        if policy == 'conservative':
            # Uncertain CDSs remain searchable even if overlapped by another call.
            for s, e in rescue_index.spans:
                s, e = max(ss, s), min(ee, e)
                if s < e:
                    targets.append((max(ss, s-boundary), min(ee, e+boundary)))
                    rows.append(dict(start0=s, end0=e, reason='uncertain_CDS_rescue', eligible=True))
    return union(targets), rows, decisions


def choose_audit(contigs, fraction, seed, explicit):
    unknown = set(explicit)-set(contigs)
    if unknown:
        raise ValueError(f'Unknown audit contigs: {sorted(unknown)}')
    ordered = sorted(contigs, key=lambda c: hashlib.sha256(f'{seed}:{c}'.encode()).hexdigest())
    return set(ordered[:math.ceil(fraction*len(ordered))]) | set(explicit)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta', type=Path, required=True)
    p.add_argument('--gff', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--policy', choices=['none','conservative','all'], default='conservative')
    p.add_argument('--min-noncoding', type=int, default=40)
    p.add_argument('--boundary-pad', type=int, default=75)
    p.add_argument('--min-cds', type=int, default=300)
    p.add_argument('--cds-confidence', type=float, default=99.0)
    p.add_argument('--window', type=int, default=3000)
    p.add_argument('--stride', type=int, default=1500)
    p.add_argument('--min-segment', type=int, default=40)
    p.add_argument('--audit-fraction', type=float, default=0.02,
                   help='Fraction of whole contigs; ceil means at least one if >0')
    p.add_argument('--audit-contigs', type=Path, help='Additional contig IDs, one per line')
    p.add_argument('--seed', type=int, default=20261016)
    p.add_argument('--loci', type=Path, help='Optional evaluation ONLY: contig,start0,end0 TSV')
    a = p.parse_args()
    if (min(a.min_noncoding, a.boundary_pad) < 0 or min(a.min_cds,a.min_segment,a.stride) < 1
        or not a.stride <= a.window <= 3000 or not 0 <= a.audit_fraction <= 1
        or not 0 <= a.cds_confidence <= 100):
        p.error('Invalid lengths, confidence, fraction, or window (maximum 3000 nt)')
    seqs = fasta(a.fasta)
    genes = read_gff(a.gff, seqs)
    by = defaultdict(list)
    for g in genes:
        by[g['contig']].append(g)
    explicit = a.audit_contigs.read_text().split() if a.audit_contigs else []
    audit = choose_audit(seqs, a.audit_fraction, a.seed, explicit)
    all_windows, selected, region_rows, gene_rows, target_rows, stats = [], [], [], [], [], []
    indexes, coverage, policy_coverage = {}, {}, {}
    for c, seq in seqs.items():
        targets, rows, decisions = regions_for(seq, by[c], a.policy, a.min_noncoding,
                                             a.boundary_pad, a.min_cds, a.cds_confidence)
        ti = indexes[c] = SpanIndex(targets)
        region_rows.extend(dict(contig=c, **r) for r in rows)
        gene_rows.extend({**g, 'parts': json.dumps(g['parts'])} for g in decisions)
        target_rows.extend(dict(contig=c, start0=s, end0=e) for s,e in targets)
        wins, policy_wins, picked = [], [], []
        for ss, ee in valid_segments(seq):
            if ee-ss < a.min_segment:
                continue
            for s,e in window_starts(ee-ss, a.window, a.stride):
                s,e = s+ss,e+ss
                overlap_nt = ti.overlap(s,e)
                keep = overlap_nt > 0 or c in audit
                row = dict(window_id=f'{c}:{s}-{e}', contig=c, start0=s, end0=e,
                    length_nt=e-s, segment_start0=ss, segment_end0=ee,
                    window_sha256=hashlib.sha256(seq[s:e].lower().encode()).hexdigest(),
                    target_overlap_nt=overlap_nt, policy_selected=overlap_nt>0,
                    audit_contig=c in audit, selected=keep,
                    reason='target' if overlap_nt else 'audit_only' if keep else 'filtered')
                all_windows.append(row)
                wins.append((s,e))
                if overlap_nt:
                    policy_wins.append((s,e))
                if keep:
                    selected.append(row)
                    picked.append((s,e))
        coverage[c], policy_coverage[c] = SpanIndex(picked), SpanIndex(policy_wins)
        stats.append(dict(contig=c, total_nt=len(seq), n_CDS=len(by[c]), audit_contig=c in audit,
            n_full_windows=len(wins), n_policy_windows=len(policy_wins), n_selected_windows=len(picked),
            full_window_nt=sum(e-s for s,e in wins), policy_window_nt=sum(e-s for s,e in policy_wins),
            selected_window_nt=sum(e-s for s,e in picked),
            target_nt=ti.prefix[-1], unique_full_scannable_nt=SpanIndex(wins).prefix[-1],
            unique_selected_nt=coverage[c].prefix[-1],
            short_unscorable_nt=sum(e-s for s,e in valid_segments(seq) if e-s<a.min_segment)))
    out = new_dir(a.out)
    fields = ['window_id','contig','start0','end0','length_nt','segment_start0','segment_end0',
              'window_sha256','target_overlap_nt','policy_selected','audit_contig','selected','reason']
    write_table(out/'windows.tsv', selected, fields)
    write_table(out/'windows_all.tsv', all_windows, fields)
    write_table(out/'regions.tsv', region_rows, ['contig','start0','end0','reason','eligible'])
    write_table(out/'targets.tsv', target_rows, ['contig','start0','end0'])
    write_table(out/'cds_decisions.tsv', gene_rows,
                ['contig','gene_id','start0','end0','strand','parts','exclude_by_policy','reason'])
    write_table(out/'contigs.tsv', stats, list(stats[0]))
    with (out/'windows.fa').open('w', encoding='utf-8') as fh:
        for r in selected:
            fh.write(f'>{r["window_id"]}\n{seqs[r["contig"]][r["start0"]:r["end0"]].lower()}\n')
    (out/'records.txt').write_text(''.join(r['window_id']+'\n' for r in selected), encoding='utf-8')
    sums = {k: sum(r[k] for r in stats) for k in stats[0] if k not in ('contig','audit_contig')}
    sums['window_nt_reduction_fraction'] = (1-sums['selected_window_nt']/sums['full_window_nt']
                                           if sums['full_window_nt'] else None)
    sums['policy_window_nt_reduction_fraction'] = (1-sums['policy_window_nt']/sums['full_window_nt']
                                                  if sums['full_window_nt'] else None)
    sums.update(audit_contigs=sorted(audit), n_selected_windows=len(selected),
                compute_note='Summed input nt is a workload proxy, not measured GPU time; two orientations.',
                no_candidates_note='No selected windows is a selection result, not evidence of absence.')
    (out/'summary.json').write_text(json.dumps(sums, indent=2)+'\n', encoding='utf-8')
    if a.loci:
        checks = []
        for r in table(a.loci):
            c,s,e = interval(r,seqs)
            checks.append({**r, 'policy_target_overlap_nt': indexes[c].overlap(s,e),
                'policy_window_covered_nt': policy_coverage[c].overlap(s,e),
                'selected_window_covered_nt': coverage[c].overlap(s,e),
                'policy_fully_covered': policy_coverage[c].overlap(s,e)==e-s,
                'selected_fully_covered': coverage[c].overlap(s,e)==e-s})
        fields_check = list(checks[0]) if checks else ['contig','start0','end0','policy_fully_covered']
        write_table(out/'locus_retention.tsv', checks, fields_check)
    provenance(out, a, [a.fasta,a.gff,a.audit_contigs,a.loci],
               limitations=['CDS confidence thresholds are uncalibrated heuristics',
                            'locus retention measures coverage, not detector recall',
                            'linear contigs only; no wrap-around detection',
                            'labels are evaluated after selection and never change it'])
    print(json.dumps(sums, indent=2))


if __name__ == '__main__':
    main()
