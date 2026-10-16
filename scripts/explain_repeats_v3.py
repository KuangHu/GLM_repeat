#!/usr/bin/env python3
"""Explainer v0.3: coordinate-linked hypotheses, retained cores and optional extensions.

No detector inference or annotation is used for discovery. All distinct hypotheses
are retained. A primary row is a DISPLAY representative, not a validated winner.
Alignment remains ungapped; exact-seed and period-range coverage limits remain.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

BASES = 'ACGT'
VERSION = '0.3'
CONFIG = dict(seed_k=12, seed_models=6, core_mismatch_fraction=0.2,
              extension_max_nt=60, extension_min_conservation=0.75,
              min_period=2, max_period=300, periodic_min_identity=0.55,
              periodic_min_units=3, pseudocount_per_base=0.25)


def read_fasta(path):
    result, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            fields = line[1:].split()
            if not fields or fields[0] in result:
                raise ValueError('Empty or duplicate FASTA identifier')
            name = fields[0]
            result[name] = []
        elif line.strip():
            if name is None:
                raise ValueError('Sequence before FASTA header')
            result[name].append(line.strip().upper())
    return {k: ''.join(v) for k, v in result.items()}


def read_table(path):
    with Path(path).open(newline='', encoding='utf-8') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_table(path, rows, fields):
    with Path(path).open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fields, delimiter='\t', extrasaction='raise')
        writer.writeheader()
        writer.writerows(rows)


def mismatches(a, b):
    if len(a) != len(b):
        raise ValueError('Ungapped distance requires equal lengths')
    return sum(x != y for x, y in zip(a, b))


def pwm_from(strings):
    if not strings or len({len(s) for s in strings}) != 1:
        raise ValueError('PWM requires nonempty equal-length copies')
    if any(set(s) - set(BASES) for s in strings):
        raise ValueError('Ambiguous bases are not imputed')
    encoded = np.array([[BASES.index(b) for b in s] for s in strings])
    counts = np.stack([(encoded == b).sum(axis=0) for b in range(4)], axis=1)
    return (counts + 0.25) / (len(strings) + 1)


def consensus(pwm):
    return ''.join(BASES[i] for i in pwm.argmax(axis=1))


def background(seq):
    counts = np.array([seq.count(b) for b in BASES], dtype=float) + 0.25
    return counts / counts.sum()


def profile(strings, bg):
    """Descriptive final distances and leave-one-copy-out evidence, not p-values.

    Copy recruitment was performed on this same region; leave-one-out does not
    remove the bias from searching for the copies in the first place.
    """
    pwm = pwm_from(strings)
    cons = consensus(pwm)
    mm = [mismatches(s, cons) for s in strings]
    loo_identity, loo_enrichment = [], []
    if len(strings) > 1:
        encoded = np.array([[BASES.index(b) for b in s] for s in strings])
        counts = np.stack([(encoded == b).sum(axis=0) for b in range(4)], axis=1)
        for s, indices in zip(strings, encoded):
            other = counts.astype(float).copy()
            other[np.arange(len(s)), indices] -= 1
            other = (other + 0.25) / len(strings)
            loo_identity.append(1 - mismatches(s, consensus(other)) / len(s))
            loo_enrichment.append(float(np.log2(
                other[np.arange(len(s)), indices] / bg[indices]).mean()))
    else:
        loo_identity, loo_enrichment = [None], [None]
    return dict(consensus=cons, pwm=pwm, mismatches=mm,
                consensus_identity=[1 - n / len(cons) for n in mm],
                loo_identity=loo_identity, loo_enrichment=loo_enrichment)


def find_copies(seq, motif, max_mm):
    k, hits, i = len(motif), [], 0
    while i <= len(seq) - k:
        piece = seq[i:i + k]
        if not (set(piece) - set(BASES)) and mismatches(piece, motif) <= max_mm:
            hits.append((i, i + k))
            i += k
        else:
            i += 1
    return hits


def extension_offsets(seq, starts, k):
    offsets = []
    for direction in (-1, 1):
        size = 0
        while size < CONFIG['extension_max_nt']:
            coords = [s-size-1 if direction == -1 else s+k+size for s in starts]
            if any(i < 0 or i >= len(seq) or seq[i] not in BASES for i in coords):
                break
            col = [seq[i] for i in coords]
            if Counter(col).most_common(1)[0][1] / len(col) < CONFIG['extension_min_conservation']:
                break
            size += 1
        offsets.append(size)
    return tuple(offsets)


def anchored_extensions(seq, cores, k):
    """Extend each existing core at its own coordinates; never recruit/delete cores."""
    left, right = extension_offsets(seq, [a for a, _ in cores], k)
    if not (left or right):
        return [], '', None
    intervals = [(a-left, b+right) for a, b in cores]
    strings = [seq[a:b] for a, b in intervals]
    template = consensus(pwm_from(strings))
    budget = max(1, round(CONFIG['core_mismatch_fraction'] * len(template)))
    attempts = []
    for (a, b), s in zip(intervals, strings):
        mm = mismatches(s, template)
        attempts.append(dict(start=a, end=b, sequence=s, search_mismatches=mm,
                             accepted=mm <= budget,
                             reason='accepted' if mm <= budget else 'extension_template_mismatch'))
    return attempts, template, budget


def motif_models(seq, min_copies=2):
    k = CONFIG['seed_k']
    counts = Counter(seq[i:i+k] for i in range(len(seq)-k+1))
    seeds = [(s, n) for s, n in counts.items()
             if n >= min_copies and len(set(s)) > 2 and not (set(s)-set(BASES))]
    # Same seed budget as v0.2; tied seeds retain genomic order.
    seeds.sort(key=lambda item: -item[1])
    result, seen = [], set()
    budget = max(1, round(CONFIG['core_mismatch_fraction'] * k))
    for seed, _ in seeds[:CONFIG['seed_models']]:
        if len(find_copies(seq, seed, 0)) < min_copies:
            continue
        initial = find_copies(seq, seed, budget)
        query = consensus(pwm_from([seq[a:b] for a, b in initial]))
        cores = find_copies(seq, query, budget)
        key = tuple(cores)
        if len(cores) < min_copies or key in seen:
            continue
        seen.add(key)
        attempts, ext_template, ext_budget = anchored_extensions(seq, cores, k)
        result.append(dict(source='motif_core', copies=cores, unit_len=k,
                           search_motif=query, search_budget=budget, seed_motif=seed,
                           extensions=attempts, extension_template=ext_template,
                           extension_budget=ext_budget, interval=None, period=None))
    return result


def period_diagnostics(seq):
    """Keep all offset diagnostics; do not promote the noisy raw-identity argmax.

    The old 0.55 identity floor remains. Excess over the local composition is
    recorded, not interpreted as significance. Harmonic alternatives are retained.
    """
    if set(seq)-set(BASES):
        raise ValueError('Period diagnostics require one unambiguous segment')
    rows = []
    chance = float((background(seq)**2).sum())
    for d in range(CONFIG['min_period'], min(CONFIG['max_period'], len(seq)//3)+1):
        n = len(seq)-d
        identity = sum(a == b for a, b in zip(seq[:-d], seq[d:])) / n
        rows.append(dict(period=d, comparisons=n, identity=identity,
                         composition_identity=chance, identity_excess=identity-chance,
                         passes_identity_floor=identity >= CONFIG['periodic_min_identity']))
    supported = [r['period'] for r in rows if r['passes_identity_floor']]
    for row in rows:
        row['supported_divisors'] = ','.join(str(d) for d in supported
                                           if d < row['period'] and row['period'] % d == 0)
    return rows


def periodic_intervals(seq, d):
    """Return every supported run, not just the longest interval in the flanks."""
    eq = np.array([a == b for a, b in zip(seq[:-d], seq[d:])], dtype=int)
    if not len(eq):
        return []
    k = min(max(2*d, 20), len(eq))
    sums = np.cumsum(np.r_[0, eq])
    good = (sums[k:]-sums[:-k])/k >= CONFIG['periodic_min_identity']
    padded = np.r_[False, good, False].astype(int)
    starts, ends = np.flatnonzero(np.diff(padded) == 1), np.flatnonzero(np.diff(padded) == -1)
    # Good windows start in [s,e); last supported compared base is e+k-2.
    intervals = [(int(s), int(e+k-1+d)) for s, e in zip(starts, ends)]
    # Sliding-window runs can describe overlapping portions of the same interval.
    merged = []
    for s, e in intervals:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(e, merged[-1][1]))
        else:
            merged.append((s, e))
    return merged


def periodic_model(seq, lo, hi, d):
    n = (hi-lo)//d
    if n < CONFIG['periodic_min_units']:
        return None
    # Phase zero is a coordinate convention; all cyclic origins remain unanchored.
    return dict(source='periodic', copies=[(lo+i*d, lo+(i+1)*d) for i in range(n)],
                unit_len=d, search_motif='', search_budget=None, seed_motif='',
                extensions=[], extension_template='', extension_budget=None,
                interval=(lo, hi), period=d)


def overlap(a, b, c, d):
    return max(0, min(b, d)-max(a, c))


def union_length(intervals):
    total, right = 0, None
    for start, end in sorted(intervals):
        total += max(0, end-max(start, right if right is not None else start))
        right = max(end, right if right is not None else end)
    return total


def classify(copies, unit_len):
    gaps = np.diff([a for a, _ in copies])
    if len(copies) < 3:
        return 'repeat_pair', float(gaps[0]) if len(gaps) else None, None
    cv = float(gaps.std()/gaps.mean())
    edge = np.median([copies[i+1][0]-copies[i][1] for i in range(len(copies)-1)])
    kind = ('short_period_tandem' if unit_len <= 10 else 'tandem_repeat') if edge <= 2 else (
        'regular_spaced_array' if cv <= 0.15 else 'variable_spacer_array' if cv <= 0.5 else 'dispersed_repeat')
    return kind, float(np.median(gaps)), cv


def analyze(seq, candidate_start, candidate_end, min_copies=2):
    """Coordinates are relative to the provided region. Ns split, never become A."""
    models, diagnostics = [], []
    for segment in re.finditer('[ACGT]+', seq):
        offset, part = segment.start(), segment.group()
        local = motif_models(part, min_copies)
        cs, ce = max(0, candidate_start-offset), min(len(part), candidate_end-offset)
        if ce > cs:
            ds = period_diagnostics(part[cs:ce])
            for r in ds:
                diagnostics.append(dict(segment_start=offset+cs, segment_end=offset+ce, **r))
                if r['passes_identity_floor']:
                    for lo, hi in periodic_intervals(part, r['period']):
                        m = periodic_model(part, lo, hi, r['period'])
                        if m:
                            local.append(m)
        for m in local:
            m['copies'] = [(a+offset, b+offset) for a, b in m['copies']]
            if m['interval']:
                m['interval'] = tuple(i+offset for i in m['interval'])
            for attempt in m['extensions']:
                attempt['start'] += offset
                attempt['end'] += offset
            models.append(m)
    bg = background(seq)
    for m in models:
        strings = [seq[a:b] for a, b in m['copies']]
        m['profile'] = profile(strings, bg)
        m['n_overlapping'] = sum(overlap(a,b,candidate_start,candidate_end)>0 for a,b in m['copies'])
        m['n_contained'] = sum(candidate_start<=a and b<=candidate_end for a,b in m['copies'])
        m['candidate_coverage_nt'] = union_length([(max(a,candidate_start),min(b,candidate_end))
            for a,b in m['copies'] if overlap(a,b,candidate_start,candidate_end)>0])
        m['association'] = 'candidate_associated' if m['n_overlapping'] else 'flank_only'
        m['supported_divisor'] = None
    for m in models:
        if m['period']:
            divisors = [n['period'] for n in models if n['period'] and n['period']<m['period']
                        and m['period']%n['period']==0 and n['association']==m['association']
                        and overlap(*m['interval'],*n['interval'])>0]
            m['supported_divisor'] = min(divisors) if divisors else None
    # This is display ordering, not an acceptance/confidence score. Every model is emitted.
    models.sort(key=lambda m: (-bool(m['n_overlapping']), bool(m['supported_divisor']), -m['candidate_coverage_nt'],
                              -float(np.mean(m['profile']['loo_enrichment'])),
                              m['unit_len'], m['source'], m['copies'][0][0]))
    status = 'candidate_associated' if any(m['n_overlapping'] for m in models) else (
        'flank_only' if models else 'unresolved')
    return models, diagnostics, status


MODEL_FIELDS = ['candidate_id','model_id','primary','status','association','contig','cand_start0','cand_end0',
    'detector_score','annotation_category','model_source','architecture','consensus','unit_len','n_copies',
    'n_copies_overlapping_candidate','n_copies_contained_in_candidate','candidate_coverage_nt',
    'extended_copies','extension_consensus','extension_unit_len','median_consensus_identity',
    'median_loo_consensus_identity','mean_loo_log2_enrichment','median_start_to_start','spacing_cv',
    'period','supported_divisor_period','periodic_interval_start0','periodic_interval_end0','phase_status','display_phase',
    'n_alternative_models','n_CDS_overlapping','CDS_overlap_fraction','annotation_available']
COPY_FIELDS = ['candidate_id','model_id','copy_id','interval_id','contig','start0','end0','sequence',
    'mismatches_to_consensus','consensus_identity','loo_consensus_identity','loo_log2_enrichment',
    'overlaps_candidate','contained_in_candidate','search_motif','search_mismatches','search_budget',
    'overlapping_peaks']
EXT_FIELDS = ['candidate_id','model_id','core_copy_id','start0','end0','sequence','accepted','reason',
    'search_motif','search_mismatches','search_budget','final_consensus','mismatches_to_consensus',
    'consensus_identity','loo_consensus_identity','loo_log2_enrichment']
CAND_FIELDS = ['candidate_id','contig','start0','end0','detector_score','annotation_category','status','reason',
    'primary_model_id','n_models','n_candidate_associated_models','n_flank_only_models','ambiguous_bases',
    'n_CDS_overlapping','CDS_overlap_fraction','annotation_available']
DIAG_FIELDS = ['candidate_id','segment_start0','segment_end0','period','comparisons','identity',
               'composition_identity','identity_excess','passes_identity_floor','supported_divisors']
PWM_FIELDS = ['model_id','layer','column0','A','C','G','T']
ALN_FIELDS = ['candidate_id','model_id','layer','copy_id','alignment']
REL_FIELDS = ['interval_id_a','interval_id_b','relation','overlap_nt']


def run(candidates, seqs, out, flank=500, min_copies=2, peaks=None):
    if flank < 0 or min_copies < 2:
        raise ValueError('flank must be nonnegative; min_copies must be at least 2')
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rows, primary, copies, extensions, cand_rows, diags, pwms, alns = [], [], [], [], [], [], [], []
    seen = set()
    for c in candidates:
        contig, cs, ce = c['contig'], int(c['start0']), int(c['end0'])
        cid = f'{contig}:{cs}-{ce}'
        if cid in seen or not 0 <= cs < ce:
            raise ValueError(f'Duplicate or invalid candidate: {cid}')
        seen.add(cid)
        g = seqs.get(contig)
        if g is not None and ce > len(g):
            raise ValueError(f'Candidate outside contig: {cid}')
        ann = {k:c.get(k,'') for k in ['n_CDS_overlapping','CDS_overlap_fraction','annotation_available']}
        s, e = (max(0,cs-flank), min(len(g),ce+flank)) if g is not None else (0,0)
        sub = g[s:e] if g is not None else ''
        models, diagnostics, status = analyze(sub,cs-s,ce-s,min_copies) if g is not None else ([],[],'unresolved')
        reason = ('sequence_unavailable' if g is None else 'no_model' if not models else
                  'only_flank_models' if status=='flank_only' else 'model_overlaps_candidate_not_biological_validation')
        n_associated = sum(bool(m['n_overlapping']) for m in models)
        cand_rows.append(dict(candidate_id=cid,contig=contig,start0=cs,end0=ce,
            detector_score=c.get('score',''),annotation_category=c.get('annotation_category',''),status=status,
            reason=reason,primary_model_id=f'{cid}/m001' if models else '',n_models=len(models),
            n_candidate_associated_models=n_associated,n_flank_only_models=len(models)-n_associated,
            ambiguous_bases=sum(b not in BASES for b in sub),**ann))
        for d in diagnostics:
            diags.append(dict(candidate_id=cid,segment_start0=s+d['segment_start'],segment_end0=s+d['segment_end'],
                              **{k:v for k,v in d.items() if k not in ('segment_start','segment_end')}))
        for idx,m in enumerate(models,1):
            mid = f'{cid}/m{idx:03d}'
            p = m['profile']
            kind,pitch,cv = classify(m['copies'],m['unit_len'])
            accepted = [x for x in m['extensions'] if x['accepted']]
            ep = profile([x['sequence'] for x in accepted],background(sub)) if accepted else None
            row = dict(candidate_id=cid,model_id=mid,primary=idx==1,status=status,association=m['association'],
                contig=contig,cand_start0=cs,cand_end0=ce,detector_score=c.get('score',''),
                annotation_category=c.get('annotation_category',''),model_source=m['source'],architecture=kind,
                consensus=p['consensus'],unit_len=m['unit_len'],n_copies=len(m['copies']),
                n_copies_overlapping_candidate=m['n_overlapping'],n_copies_contained_in_candidate=m['n_contained'],
                candidate_coverage_nt=m['candidate_coverage_nt'],extended_copies=len(accepted),
                extension_consensus=ep['consensus'] if ep else '',extension_unit_len=len(ep['consensus']) if ep else '',
                median_consensus_identity=float(np.median(p['consensus_identity'])),
                median_loo_consensus_identity=float(np.median(p['loo_identity'])),
                mean_loo_log2_enrichment=float(np.mean(p['loo_enrichment'])),
                median_start_to_start=pitch,spacing_cv=cv,period=m['period'],supported_divisor_period=m['supported_divisor'],
                periodic_interval_start0=s+m['interval'][0] if m['interval'] else '',
                periodic_interval_end0=s+m['interval'][1] if m['interval'] else '',
                phase_status='unanchored_cyclic' if m['period'] else 'not_applicable',
                display_phase=0 if m['period'] else '',n_alternative_models=len(models)-1,**ann)
            rows.append(row)
            if idx==1:
                primary.append(row)
            for layer,pp in [('core' if m['source']=='motif_core' else 'periodic',p),('extension',ep)]:
                if pp:
                    pwms.extend(dict(model_id=mid,layer=layer,column0=j,**dict(zip(BASES,map(float,col))))
                                for j,col in enumerate(pp['pwm']))
            extidx = 0
            for j,(a,b) in enumerate(m['copies']):
                copyid = f'{mid}/c{j+1:04d}'
                start,end = s+a,s+b
                seq = sub[a:b]
                copies.append(dict(candidate_id=cid,model_id=mid,copy_id=copyid,
                    interval_id=f'{contig}:{start}-{end}',contig=contig,start0=start,end0=end,sequence=seq,
                    mismatches_to_consensus=p['mismatches'][j],consensus_identity=p['consensus_identity'][j],
                    loo_consensus_identity=p['loo_identity'][j],loo_log2_enrichment=p['loo_enrichment'][j],
                    overlaps_candidate=overlap(start,end,cs,ce)>0,contained_in_candidate=cs<=start and end<=ce,
                    search_motif=m['search_motif'],search_mismatches=mismatches(seq,m['search_motif']) if m['search_motif'] else '',
                    search_budget=m['search_budget'],overlapping_peaks=(sum(overlap(start,end,x,y)>0 for x,y in peaks.get(contig,[]))
                    if peaks is not None else '')))
                alns.append(dict(candidate_id=cid,model_id=mid,layer='core' if m['source']=='motif_core' else 'periodic',
                                 copy_id=copyid,alignment=seq))
                if m['extensions']:
                    x = m['extensions'][j]
                    extensions.append(dict(candidate_id=cid,model_id=mid,core_copy_id=copyid,start0=s+x['start'],
                        end0=s+x['end'],sequence=x['sequence'],accepted=x['accepted'],reason=x['reason'],
                        search_motif=m['extension_template'],search_mismatches=x['search_mismatches'],
                        search_budget=m['extension_budget'],final_consensus=ep['consensus'] if ep else '',
                        mismatches_to_consensus=mismatches(x['sequence'],ep['consensus']) if ep else '',
                        consensus_identity=1-mismatches(x['sequence'],ep['consensus'])/len(x['sequence']) if ep else '',
                        loo_consensus_identity=ep['loo_identity'][extidx] if x['accepted'] else '',
                        loo_log2_enrichment=ep['loo_enrichment'][extidx] if x['accepted'] else ''))
                    if x['accepted']:
                        alns.append(dict(candidate_id=cid,model_id=mid,layer='extension',copy_id=copyid,alignment=x['sequence']))
                        extidx += 1
    # Shared coordinates have stable IDs. Shifted/overlapping intervals are relations,
    # not silently collapsed into a presumed single biological copy.
    unique = {r['interval_id']:(r['contig'],r['start0'],r['end0']) for r in copies}
    relations, by_contig = [], defaultdict(list)
    for iid,(contig,a,b) in unique.items():
        by_contig[contig].append((a,b,iid))
    for intervals in by_contig.values():
        intervals.sort()
        for i,(a,b,iid) in enumerate(intervals):
            for c,d,jid in intervals[i+1:]:
                if c>=b:
                    break
                relations.append(dict(interval_id_a=iid,interval_id_b=jid,relation='overlap_identity_unresolved',
                                      overlap_nt=overlap(a,b,c,d)))
    for name,data,fields in [('models.tsv',primary,MODEL_FIELDS),('hypotheses.tsv',rows,MODEL_FIELDS),
        ('copies.tsv',copies,COPY_FIELDS),('extensions.tsv',extensions,EXT_FIELDS),
        ('candidates.tsv',cand_rows,CAND_FIELDS),('period_diagnostics.tsv',diags,DIAG_FIELDS),
        ('pwms.tsv',pwms,PWM_FIELDS),('alignments.tsv',alns,ALN_FIELDS),('copy_relations.tsv',relations,REL_FIELDS)]:
        write_table(out/name,data,fields)
    states = Counter(r['status'] for r in cand_rows)
    summary = dict(version=VERSION,n_candidates=len(candidates),status_counts=dict(states),
        n_primary_models=len(primary),n_hypotheses=len(rows),n_copy_evidence_rows=len(copies),
        n_unique_coordinate_intervals=len(unique),n_extension_attempts=len(extensions),
        n_extensions_accepted=sum(r['accepted'] for r in extensions),
        n_multiple_associated_hypotheses=sum(r['n_candidate_associated_models']>1 for r in cand_rows),
        parameters={**CONFIG,'flank':flank,'min_copies':min_copies},peaks_available=peaks is not None,
        display_selection='association, non-harmonic before supported-divisor alternatives, candidate covered nt, descriptive LOO enrichment, shorter unit; all hypotheses retained',
        support_limit='LOO excludes a copy from its profile, not from discovery/selection; not calibrated probability',
        phase_policy='phase 0 is interval-relative display convention; cyclic biological origin unresolved',
        limitations=['ungapped_substitutions_only','exact_12mer_seed_for_motif_branch','period_range_2_300',
                     'period_identity_floor_not_a_significance_test','candidate_association_not_validation'],
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidates',type=Path,required=True)
    p.add_argument('--genomes',type=Path,nargs='+',required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--flank',type=int,default=500)
    p.add_argument('--min-copies',type=int,default=2)
    p.add_argument('--peaks',type=Path)
    a = p.parse_args()
    seqs, sources = {}, [a.candidates]
    for path in a.genomes:
        for f in ([path] if path.is_file() else sorted(path.glob('*.fna'))):
            new = read_fasta(f)
            if set(new)&set(seqs):
                raise ValueError('Duplicate contig across genome inputs')
            seqs.update(new)
            sources.append(f)
    peaks = None
    if a.peaks:
        peaks = defaultdict(list)
        for r in read_table(a.peaks):
            peaks[r['contig']].append((int(r['start0']),int(r['end0'])))
        sources.append(a.peaks)
    summary = run(read_table(a.candidates),seqs,a.out,a.flank,a.min_copies,peaks)
    inputs = {str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources}
    (a.out/'inputs.sha256.json').write_text(json.dumps(inputs,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:summary[k] for k in ['version','n_candidates','status_counts','n_hypotheses']},indent=2))


if __name__ == '__main__':
    main()
