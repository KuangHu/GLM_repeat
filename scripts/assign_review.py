#!/usr/bin/env python3
"""Apply explicit review rules to the browse book and record the reasoning.

These are ASSISTANT judgements, written to AUTO_REVIEW_* columns. The HUMAN_*
columns stay blank: a machine verdict must not occupy the field reserved for the
reviewer's own reading.

Every verdict is produced by a stated rule and carries the evidence that fired
it, so a reviewer can disagree with any single call by checking the same numbers.
Three independent questions are answered separately, because a candidate can have
solid repeat evidence and a useless explanation, or no repeat at all and still be
a biologically interesting locus.

  repeat_evidence   does the sequence itself show a repeat?
  explanation       does the explainer's model account for THIS candidate?
  followup          is it worth looking at further?

`cannot_judge` is a first-class outcome and is never avoided for tidiness.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


def tsv(path):
    with Path(path).open() as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def num(v, d=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return d


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--digest', type=Path, required=True)
    ap.add_argument('--chance', type=Path, required=True)
    ap.add_argument('--trna', action='append', required=True, metavar='LABEL=FILE')
    ap.add_argument('--crispr', action='append', default=[], metavar='LABEL=FILE')
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()

    chance = {(r['run'], r['candidate_id']): r for r in tsv(a.chance)}
    trna = {}
    for spec in a.trna:
        lab, _, f = spec.partition('=')
        for r in tsv(f):
            trna[(lab, r['candidate_id'])] = r

    crispr = {}
    for spec in a.crispr:
        lab, _, f = spec.partition('=')
        for r in tsv(f):
            crispr[(lab, r['candidate_id'])] = r

    rows = []
    for d in tsv(a.digest):
        k = (d['run'], d['candidate_id'])
        ch, tr, cr = chance[k], trna.get(k, {}), crispr.get(k, {})
        n_cr = int(cr.get('n_CRISPR_arrays_overlapping', 0) or 0)
        cr_frac = num(cr.get('CRISPR_coverage_fraction'), 0.0)
        cr_off = cr.get('CRISPR_repeat_without_called_array', '') == 'True'
        surprising = ch['surprising_at_0.05'] == 'True'
        p = num(ch['p_empirical'], 1.0)
        n_trna = int(tr.get('n_tRNA_overlapping', 0) or 0)
        trna_frac = num(tr.get('tRNA_coverage_fraction'), 0.0)
        status = d['status']
        arch = d['architecture']
        ncop = int(d['primary_n_copies'] or 0)
        nin = int(d['primary_copies_inside'] or 0)
        mid = num(d['median_loo_identity'], 0.0)
        probe = d['probe_verdict']
        # The envelope can be thousands of nt, which dilutes a low-complexity
        # consensus to invisibility, so the consensus is judged on its own too.
        cons_H = num(d['consensus_entropy'], 2.0)
        unit = int(num(d['unit_len'], 0))
        lowc = (num(d['cand_top3_3mer_frac']) > 0.30
                or int(d['cand_max_homopolymer'] or 0) >= 8
                or num(d['cand_entropy_bits']) < 1.70
                or (d['consensus'] and cons_H < 1.20))
        # a model whose copies barely agree is not a supported model
        weak_model = bool(d['consensus']) and (
            (unit <= 4 and mid < 0.70) or (unit >= 6 and mid < 0.60))
        ev, why = [], []

        # ---------------- repeat evidence -------------------------------
        if status == 'candidate_associated' and ncop >= 3 and nin >= 2 and mid >= 0.83 and surprising:
            rep = 'clear_repeat'
            why.append(f'{ncop} copies ({nin} inside), LOO identity {mid}, p={p}')
        elif status == 'candidate_associated' and nin >= 2 and surprising:
            rep = 'suspected_repeat'
            why.append(f'only {ncop} copies ({nin} inside) but excess 12-mers p={p}')
        elif status == 'candidate_associated' and surprising:
            rep = 'suspected_repeat'
            why.append(f'model copies mostly outside the envelope ({nin}/{ncop} inside), p={p}')
        elif status == 'flank_only':
            rep = 'cannot_judge'
            why.append('model sits entirely in the flanks; it is not evidence about '
                       'this candidate either way')
        elif status == 'no_model' and not surprising and probe in (
                'no_repeated_seed', 'seeds_without_extension'):
            rep = 'no_repeat_found'
            why.append(f'no model, 12-mer content at chance (p={p}), probe={probe}')
        elif status == 'no_model' and not surprising:
            rep = 'no_repeat_found'
            why.append(f'no model and 12-mer content at chance (p={p})')
        else:
            rep = 'cannot_judge'
            why.append(f'status={status}, surprising={surprising}, p={p}')
        if lowc and rep in ('clear_repeat', 'suspected_repeat'):
            rep = 'suspected_repeat'
            why.append('LOW COMPLEXITY: entropy '
                       f"{d['cand_entropy_bits']}, top-3 3-mer {d['cand_top3_3mer_frac']}, "
                       f"homopolymer {d['cand_max_homopolymer']} -- a composition artefact "
                       'cannot be excluded')

        # ---------------- explanation completeness ----------------------
        if weak_model:
            exp = 'structure_unsupported'
            why.append(f'unit {unit} nt with median LOO identity {mid}: the copies do '
                       'not agree well enough to support this model')
        elif n_cr and cr_frac >= 0.3:
            exp = 'model_is_fragment_of_known_array'
            why.append(f"candidate lies in a pilercr CRISPR array "
                       f"({cr.get('CRISPR_array_copies')} copies, repeat "
                       f"{cr.get('CRISPR_repeat_len')} nt / spacer "
                       f"{cr.get('CRISPR_spacer_len')} nt); the 12-nt model consensus is a "
                       'fragment of the full CRISPR repeat, so it describes a real array '
                       'but does not represent its true unit')
        elif n_trna and trna_frac >= 0.3:
            exp = 'explainer_wrong_model_class'
            why.append(f'overlaps {n_trna} annotated tRNA ({trna_frac:.0%} of the '
                       'envelope); the model describes conservation shared between tRNA '
                       'genes, not a repeat element')
        elif n_trna and trna_frac >= 0.3 and status in ('no_model', 'flank_only'):
            exp = 'explainer_wrong_model_class'
            why.append(f'overlaps {n_trna} annotated tRNA ({trna_frac:.0%} of the '
                       'envelope); a repeat explainer cannot represent a structured '
                       'ncRNA gene, so this is a category mismatch, not a failure')
        elif status == 'no_model' and rep == 'no_repeat_found':
            exp = 'nothing_to_explain'
            why.append('no repeat evidence, so there is no model to be missing')
        elif status == 'no_model':
            exp = 'no_model_but_signal_present'
        elif status == 'flank_only':
            exp = 'boundary_or_phase_problem'
        elif nin >= max(2, ncop - 1) and rep == 'clear_repeat':
            exp = 'model_adequate'
        elif nin >= 1:
            exp = 'partial'
            why.append(f'{nin}/{ncop} model copies lie inside the envelope')
        else:
            exp = 'structure_unsupported'
        if status == 'candidate_associated' and not surprising:
            exp = 'structure_unsupported'
            why.append(f'model proposed but 12-mer excess is at chance (p={p})')

        # ---------------- worth following up ----------------------------
        if cr_off:
            fup = 'worth_followup'
            why.append('carries this genome\'s CRISPR repeat consensus '
                       f"({cr.get('model_consensus_in_CRISPR_repeat')}) but falls OUTSIDE "
                       'every pilercr-called array: a possible degenerate or fragmentary '
                       'array that the dedicated tool did not call')
        elif n_cr and cr_frac >= 0.3:
            fup = 'known_common_structure'
            why.append('inside an annotated CRISPR array')
        elif n_trna and trna_frac >= 0.3:
            fup = 'known_common_structure'
            why.append('tRNA locus: known, abundant, well characterised')
        elif lowc:
            fup = 'known_common_structure'
            why.append('low-complexity tract')
        elif rep == 'clear_repeat' and not n_trna:
            fup = 'worth_followup'
        elif rep == 'no_repeat_found':
            fup = 'not_worth_followup'
        else:
            fup = 'function_unknown'

        rows.append({
            'run': d['run'], 'candidate_id': d['candidate_id'],
            'in_random_sample': d['in_random_sample'], 'strata': d['strata'],
            'score': d['score'], 'span_nt': d['span_nt'],
            'status': status, 'architecture': arch,
            'primary_n_copies': ncop, 'primary_copies_inside': nin,
            'median_loo_identity': d['median_loo_identity'],
            'seed_chance_p': ch['p_empirical'],
            'seed_excess_surprising': surprising,
            'n_CRISPR_arrays_overlapping': n_cr,
            'CRISPR_coverage_fraction': cr.get('CRISPR_coverage_fraction', ''),
            'CRISPR_repeat_len': cr.get('CRISPR_repeat_len', ''),
            'CRISPR_spacer_len': cr.get('CRISPR_spacer_len', ''),
            'CRISPR_repeat_without_called_array': cr_off,
            'weak_model_flag': weak_model,
            'n_tRNA_overlapping': n_trna,
            'tRNA_coverage_fraction': tr.get('tRNA_coverage_fraction', ''),
            'tRNA_names': tr.get('tRNA_names', ''),
            'low_complexity_flag': lowc,
            'probe_verdict': probe,
            'AUTO_REVIEW_repeat_evidence': rep,
            'AUTO_REVIEW_explanation': exp,
            'AUTO_REVIEW_followup': fup,
            'AUTO_REVIEW_reasoning': '; '.join(why),
            'HUMAN_sequence_evidence': '', 'HUMAN_explanation_coverage': '',
            'HUMAN_biological_context': '', 'HUMAN_notes': '',
        })

    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    def tally(sub, field):
        return dict(Counter(r[field] for r in sub).most_common())
    rnd = [r for r in rows if r['in_random_sample'].lower() == 'true']
    enr = [r for r in rows if r['in_random_sample'].lower() != 'true']
    summary = {
        'note': 'ASSISTANT judgements in AUTO_REVIEW_* columns; HUMAN_* left blank. '
                'Random and enriched strata are reported separately and must never '
                'be pooled: only the random stratum describes what a reviewer would '
                'encounter.',
        'random_stratum': {'n': len(rnd),
                           'repeat_evidence': tally(rnd, 'AUTO_REVIEW_repeat_evidence'),
                           'explanation': tally(rnd, 'AUTO_REVIEW_explanation'),
                           'followup': tally(rnd, 'AUTO_REVIEW_followup'),
                           'tRNA_overlapping': sum(1 for r in rnd if r['n_tRNA_overlapping']),
                           'CRISPR_overlapping': sum(1 for r in rnd if r['n_CRISPR_arrays_overlapping']),
                           'CRISPR_repeat_off_array': sum(1 for r in rnd if r['CRISPR_repeat_without_called_array'])},
        'enriched_stratum': {'n': len(enr),
                             'repeat_evidence': tally(enr, 'AUTO_REVIEW_repeat_evidence'),
                             'explanation': tally(enr, 'AUTO_REVIEW_explanation'),
                             'followup': tally(enr, 'AUTO_REVIEW_followup'),
                             'tRNA_overlapping': sum(1 for r in enr if r['n_tRNA_overlapping']),
                             'CRISPR_overlapping': sum(1 for r in enr if r['n_CRISPR_arrays_overlapping']),
                             'CRISPR_repeat_off_array': sum(1 for r in enr if r['CRISPR_repeat_without_called_array'])},
    }
    a.out.with_suffix('.summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
