#!/usr/bin/env python3
"""Select the ~60-candidate browse book for human review.

Two strata, kept strictly separate because they answer different questions:

  random    20 candidates drawn uniformly at random from EACH round-2 route with a
            fixed seed. This is the only stratum that supports an unbiased
            statement about what a reviewer would see. A candidate drawn here
            stays in the random denominator even if an enriched rule also picks
            it; `in_random_sample` is recorded independently of `rules`.

  enriched  candidates chosen by explicit, named rules so the reviewer sees the
            range of things the pipeline produces: clear structures, high copy
            number, different gene neighbourhoods, high-scoring candidates the
            explainer could not model, and different source lineages. Enriched
            selection deliberately uses detector score, repeat structure and
            neighbouring gene products -- it is a guided tour, NOT a sample, and
            nothing about yield may be estimated from it.

Every selected row records every rule that fired, so the two strata can always be
re-separated downstream.
"""
from __future__ import annotations

import argparse
import csv
import json
import random
from collections import defaultdict
from pathlib import Path

DEFAULT_RANDOM_ROUTES = ('r2_host', 'r2_phage')


def table(path):
    p = Path(path)
    if not p.exists():
        return []
    with p.open() as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def load_run(name, root):
    root = Path(root)
    cand = table(root / 'calls/candidates.tsv')
    models = {}
    for m in table(root / 'explain_full/models.tsv'):
        if m.get('primary', '').lower() in ('true', '1', 'yes'):
            models[m['candidate_id']] = m
    probe = {r['candidate_id']: r for r in table(root / 'no_model_probe/no_model_probe.tsv')}
    genes = defaultdict(list)
    for g in table(root / 'neighborhoods/genes.tsv'):
        genes[g['candidate_id']].append(g)
    for c in cand:
        c['run'] = name
        c['_model'] = models.get(c['candidate_id'])
        c['_probe'] = probe.get(c['candidate_id'])
        c['_genes'] = genes.get(c['candidate_id'], [])
    return cand


def neighborhood_class(genes):
    """Structural neighbourhood class.

    The CDS set is pyrodigal ab initio prediction: every `product` is
    "hypothetical protein", so no functional signature exists in this data. The
    only honest neighbourhood descriptor available is the ARCHITECTURE around the
    candidate -- whether it sits inside a coding region, and how the flanking
    genes are oriented. This is a display/selection label only; it asserts nothing
    about function.
    """
    if not genes:
        return 'no_called_CDS'
    if any(g.get('relation') == 'overlap' for g in genes):
        return 'cds_overlapping'
    dist = lambda g: int(float(g.get('distance_nt') or 0))
    left = sorted((g for g in genes if g.get('relation') == 'left'), key=dist)
    right = sorted((g for g in genes if g.get('relation') == 'right'), key=dist)
    if not left or not right:
        return 'intergenic_one_sided'
    ls, rs = left[0].get('strand', '?'), right[0].get('strand', '?')
    if ls == '-' and rs == '+':
        return 'intergenic_divergent'
    if ls == '+' and rs == '-':
        return 'intergenic_convergent'
    if ls == rs == '+':
        return 'intergenic_tandem_forward'
    if ls == rs == '-':
        return 'intergenic_tandem_reverse'
    return 'intergenic_unclassified'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', action='append', required=True,
                   metavar='NAME=DIR', help='repeatable, e.g. round2_host=/path/to/run')
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--seed', type=int, default=20261016)
    p.add_argument('--n-random', type=int, default=20)
    p.add_argument('--n-lineage', type=int, default=5,
                   help='bounded number of additional source genomes to represent')
    p.add_argument('--random-route', action='append', dest='random_routes',
                   help='run label to draw the random stratum from (repeatable)')
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    random_routes = tuple(a.random_routes or DEFAULT_RANDOM_ROUTES)

    runs = {}
    for spec in a.run:
        name, _, d = spec.partition('=')
        runs[name] = load_run(name, d)
    allc = [c for cs in runs.values() for c in cs]
    key = lambda c: (c['run'], c['candidate_id'])
    index = {key(c): c for c in allc}
    rules = defaultdict(list)
    random_keys = set()

    # ---------- random stratum (fixed seed, deterministic ordering) ----------
    rng = random.Random(a.seed)
    random_lists = {}
    for route in random_routes:
        pool = sorted(runs.get(route, []), key=lambda c: c['candidate_id'])
        drawn = rng.sample(pool, min(a.n_random, len(pool)))
        drawn.sort(key=lambda c: c['candidate_id'])
        random_lists[route] = [c['candidate_id'] for c in drawn]
        for c in drawn:
            random_keys.add(key(c))
            rules[key(c)].append(f'random:{route}')

    # ---------- enriched rules ----------
    num = lambda v, d=0.0: float(v) if v not in (None, '', '.') else d
    integer = lambda v, d=0: int(float(v)) if v not in (None, '', '.') else d
    modelled = [c for c in allc if c['_model']]

    # 1. clear structure: best-supported primary model per architecture
    by_arch = defaultdict(list)
    for c in modelled:
        by_arch[c['_model'].get('architecture', 'unspecified')].append(c)
    for arch, cs in sorted(by_arch.items()):
        cs.sort(key=lambda c: (-num(c['_model'].get('median_loo_consensus_identity')),
                               -integer(c['_model'].get('n_copies')), c['candidate_id']))
        for c in cs[:1]:
            rules[key(c)].append(f'clear_structure:{arch}')

    # 2. multi-copy: highest copy number per run
    for name, cs in sorted(runs.items()):
        m = [c for c in cs if c['_model']]
        m.sort(key=lambda c: (-integer(c['_model'].get('n_copies')), c['candidate_id']))
        for c in m[:2]:
            rules[key(c)].append(f'multi_copy:{name}')

    # 3. different gene neighbourhoods: one representative per structural class
    by_class = defaultdict(list)
    for c in allc:
        by_class[neighborhood_class(c['_genes'])].append(c)
    for cls, cs in sorted(by_class.items()):
        cs.sort(key=lambda c: (-num(c['score']), c['run'], c['candidate_id']))
        rules[key(cs[0])].append(f'gene_neighborhood:{cls}')

    # 4. high-scoring candidates the explainer could not model
    for name, cs in sorted(runs.items()):
        nm = [c for c in cs if not c['_model']]
        nm.sort(key=lambda c: (-num(c['score']), c['candidate_id']))
        for c in nm[:2]:
            rules[key(c)].append(f'high_score_no_model:{name}')

    # 5. lineage spread: a bounded set of source genomes not yet represented.
    # The phage routes contain ~65 genomes each, so "one per contig" would swamp
    # the book. Uncovered contigs are sorted by accession and sampled at evenly
    # spaced indices, which spans the list deterministically without using any
    # detector score, repeat structure or gene content to decide WHICH lineages
    # appear; the score is used only to pick a representative within a lineage.
    covered = {index[k]['contig'] for k in rules}
    by_contig = defaultdict(list)
    for c in allc:
        if c['contig'] not in covered:
            by_contig[c['contig']].append(c)
    uncovered = sorted(by_contig)
    if uncovered:
        n = min(a.n_lineage, len(uncovered))
        picks = [uncovered[round(i * (len(uncovered) - 1) / max(1, n - 1))]
                 for i in range(n)] if n > 1 else uncovered[:1]
        for contig in sorted(set(picks)):
            cs = sorted(by_contig[contig],
                        key=lambda c: (-num(c['score']), c['candidate_id']))
            rules[key(cs[0])].append('lineage_spread')

    # ---------- emit ----------
    rows = []
    for k, rs in rules.items():
        c = index[k]
        m = c['_model'] or {}
        pr = c['_probe'] or {}
        rows.append({
            'run': c['run'], 'candidate_id': c['candidate_id'],
            'contig': c['contig'], 'start0': c['start0'], 'end0': c['end0'],
            'span_nt': c['span_nt'], 'score': c['score'], 'n_peaks': c.get('n_peaks', ''),
            'in_random_sample': k in random_keys,
            'n_rules': len(rs), 'rules': ','.join(sorted(rs)),
            'explainer_status': m.get('status', 'no_model'),
            'architecture': m.get('architecture', ''),
            'unit_len': m.get('unit_len', ''), 'n_copies': m.get('n_copies', ''),
            'median_loo_consensus_identity': m.get('median_loo_consensus_identity', ''),
            'probe_verdict': pr.get('verdict', ''),
            'n_neighbor_CDS': len(c['_genes']),
            'neighborhood_class': neighborhood_class(c['_genes']),
            'HUMAN_sequence_evidence': '', 'HUMAN_explanation_coverage': '',
            'HUMAN_biological_context': '', 'HUMAN_notes': '',
        })
    rows.sort(key=lambda r: (not r['in_random_sample'], r['run'], r['candidate_id']))
    with (a.out / 'browse_selection.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    counts = defaultdict(int)
    for r in rows:
        for s in r['rules'].split(','):
            counts[s.split(':')[0]] += 1
    prov = {
        'seed': a.seed,
        'n_random_per_route': a.n_random,
        'random_routes': list(random_routes),
        'random_sampling_lists': random_lists,
        'random_draw_procedure':
            'per route: candidates sorted by candidate_id, then random.Random(seed).sample '
            'drawn once per route in the order listed in random_routes; one shared RNG, so '
            'the draw is reproducible only for this exact route order and pool',
        'runs': {k: len(v) for k, v in runs.items()},
        'selected_total': len(rows),
        'selected_random': sum(r['in_random_sample'] for r in rows),
        'selected_enriched_only': sum(1 for r in rows if not r['in_random_sample']),
        'rule_counts': dict(counts),
        'n_lineage_cap': a.n_lineage,
        'annotation_limitation':
            'CDS are pyrodigal ab initio predictions; every product is '
            '"hypothetical protein". No functional annotation exists in this data, '
            'so neighbourhood classes are structural (overlap / flank orientation) '
            'and carry no functional claim.',
        'interpretation_rule':
            'Only the random stratum supports an unbiased statement about what a '
            'reviewer would encounter. The enriched stratum is a guided tour selected '
            'using detector score, repeat structure and neighbouring gene products; no '
            'yield, precision or FDR may be computed from it. Random IDs remain in the '
            'random denominator even when an enriched rule also selected them.',
        'human_fields_blank': ['HUMAN_sequence_evidence', 'HUMAN_explanation_coverage',
                               'HUMAN_biological_context', 'HUMAN_notes'],
    }
    (a.out / 'selection_provenance.json').write_text(json.dumps(prov, indent=2) + '\n')
    print(f"{len(rows)} candidates selected "
          f"({prov['selected_random']} random, {prov['selected_enriched_only']} enriched-only)")
    for k, v in sorted(counts.items()):
        print(f'  {k:24s} {v}')


if __name__ == '__main__':
    main()
