#!/usr/bin/env python3
"""Build the offline review package for round-1 + round-2.

Produces a self-contained directory that opens from START_HERE.html with relative
links only: no cluster absolute paths, no CDN, no network service. Figures are
written as inline SVG (and PNG where a raster is useful). Every automatic label is
prefixed AUTO_; every human judgement column is left blank.

The package never asserts that a candidate is a confirmed repeat: all candidates
are computational and unconfirmed until a human fills the HUMAN_ columns.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import random
import re
from collections import defaultdict
from pathlib import Path

ESC = html.escape


def read_fasta(path):
    out, name = {}, None
    for line in Path(path).open():
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split()[0]
            out[name] = []
        elif line:
            out[name].append(line)
    return {k: ''.join(v).upper() for k, v in out.items()}


def tsv(path):
    p = Path(path)
    return list(csv.DictReader(p.open(), delimiter='\t')) if p.exists() else []


# ----------------------------------------------------------------- figures
def svg_track(values, width=900, height=120, cand=None, peaks=(), copies=(),
              x0=0, title=''):
    """Detector score track. x axis is genomic, y is nats above uniform."""
    if not values:
        return '<p>no track</p>'
    n = len(values)
    lo, hi = min(values), max(values)
    hi = max(hi, 1.3863)
    rng = (hi - lo) or 1.0
    step = max(1, n // width)
    pts = []
    for i in range(0, n, step):
        x = i / n * width
        y = height - (values[i] - lo) / rng * (height - 24) - 12
        pts.append(f'{x:.1f},{y:.1f}')
    ln4 = height - (1.3863 - lo) / rng * (height - 24) - 12
    zero = height - (0 - lo) / rng * (height - 24) - 12
    bands = ''
    for s, e in copies:
        bands += (f'<rect x="{(s - x0) / n * width:.1f}" y="0" '
                  f'width="{max(1.2, (e - s) / n * width):.1f}" height="{height}" '
                  f'fill="#2563eb" opacity="0.16"/>')
    for s, e in peaks:
        bands += (f'<rect x="{(s - x0) / n * width:.1f}" y="{height - 8}" '
                  f'width="{max(1.2, (e - s) / n * width):.1f}" height="6" fill="#dc2626"/>')
    cb = ''
    if cand:
        cs, ce = cand
        cb = (f'<rect x="{(cs - x0) / n * width:.1f}" y="0" '
              f'width="{max(2, (ce - cs) / n * width):.1f}" height="{height}" '
              f'fill="none" stroke="#111827" stroke-dasharray="4 3" stroke-width="1.4"/>')
    return f'''<svg viewBox="0 0 {width} {height + 34}" width="100%" role="img"
 aria-label="{ESC(title)}" xmlns="http://www.w3.org/2000/svg">
<rect width="{width}" height="{height}" fill="#fbfbfd"/>{bands}{cb}
<line x1="0" y1="{zero:.1f}" x2="{width}" y2="{zero:.1f}" stroke="#9ca3af" stroke-dasharray="2 3"/>
<line x1="0" y1="{ln4:.1f}" x2="{width}" y2="{ln4:.1f}" stroke="#f59e0b" stroke-dasharray="5 3"/>
<polyline fill="none" stroke="#111827" stroke-width="1.1" points="{' '.join(pts)}"/>
<text x="4" y="{height + 13}" font-size="11" fill="#374151">{x0:,} nt</text>
<text x="{width - 4}" y="{height + 13}" font-size="11" text-anchor="end" fill="#374151">{x0 + n:,} nt</text>
<text x="4" y="{height + 29}" font-size="11" fill="#b45309">橙线 = ln(4)=1.386 上限</text>
<text x="{width / 2}" y="{height + 29}" font-size="11" text-anchor="middle" fill="#6b7280">红条 = peak，蓝带 = repeat copy，虚线框 = 候选范围</text>
</svg>'''


def svg_genes(genes, lo, hi, width=900, height=92, cand=None):
    """Gene neighbourhood with explicit coordinates and strand arrows."""
    if hi <= lo:
        return ''
    span = hi - lo
    rows = []
    if cand:
        cs, ce = cand
        rows.append(f'<rect x="{(cs - lo) / span * width:.1f}" y="8" '
                    f'width="{max(2, (ce - cs) / span * width):.1f}" height="{height - 34}" '
                    f'fill="#2563eb" opacity="0.10" stroke="#2563eb" stroke-dasharray="3 2"/>')
    for g in genes:
        s, e, strand, name = g
        x = (max(s, lo) - lo) / span * width
        w = max(3.0, (min(e, hi) - max(s, lo)) / span * width)
        y = 30 if strand == '+' else 52
        tip = 6 if w > 10 else 0
        if strand == '+':
            pts = f'{x},{y} {x + w - tip},{y} {x + w},{y + 8} {x + w - tip},{y + 16} {x},{y + 16}'
        else:
            pts = f'{x + w},{y} {x + tip},{y} {x},{y + 8} {x + tip},{y + 16} {x + w},{y + 16}'
        rows.append(f'<polygon points="{pts}" fill="#64748b" opacity="0.85"/>')
        if w > 54:
            rows.append(f'<text x="{x + w / 2:.1f}" y="{y + 12}" font-size="9" '
                        f'text-anchor="middle" fill="#ffffff">{ESC(name[:14])}</text>')
    return f'''<svg viewBox="0 0 {width} {height}" width="100%" role="img"
 aria-label="gene neighbourhood" xmlns="http://www.w3.org/2000/svg">
<rect width="{width}" height="{height}" fill="#ffffff"/>
<line x1="0" y1="46" x2="{width}" y2="46" stroke="#cbd5e1"/>{''.join(rows)}
<text x="4" y="14" font-size="11" fill="#374151">{lo:,}</text>
<text x="{width - 4}" y="14" font-size="11" text-anchor="end" fill="#374151">{hi:,}</text>
<text x="4" y="{height - 4}" font-size="10" fill="#6b7280">上行 = + 链，下行 = - 链；灰框为预测 CDS</text>
</svg>'''


HEAD = '''<!doctype html><html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{t}</title>
<style>
:root{{--fg:#111827;--mut:#6b7280;--line:#e5e7eb;--bg:#ffffff;--card:#f9fafb;--acc:#2563eb}}
*{{box-sizing:border-box}}
body{{margin:0;padding:24px;font:15px/1.65 -apple-system,"Segoe UI","Noto Sans CJK SC","Microsoft YaHei",sans-serif;color:var(--fg);background:var(--bg);max-width:1180px;margin-inline:auto}}
h1{{font-size:24px;margin:0 0 6px}} h2{{font-size:18px;margin:28px 0 10px;padding-bottom:6px;border-bottom:1px solid var(--line)}}
h3{{font-size:15px;margin:18px 0 8px}}
a{{color:var(--acc)}} code{{background:var(--card);padding:1px 5px;border-radius:4px;font-size:13px;word-break:break-all}}
table{{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0}}
th,td{{border:1px solid var(--line);padding:6px 9px;text-align:left;vertical-align:top}}
th{{background:var(--card);position:sticky;top:0;cursor:pointer}}
.mut{{color:var(--mut)}} .card{{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:14px;margin:12px 0}}
.warn{{background:#fff7ed;border:1px solid #fed7aa;border-radius:8px;padding:12px;margin:12px 0}}
.seq{{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:12px;white-space:pre-wrap;word-break:break-all;background:var(--card);padding:10px;border-radius:6px;line-height:1.5}}
.tag{{display:inline-block;background:#eef2ff;color:#3730a3;border-radius:999px;padding:1px 9px;font-size:12px;margin-right:4px}}
input,select{{padding:6px 8px;border:1px solid var(--line);border-radius:6px;font-size:13px;margin:0 6px 6px 0}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:10px}}
.kv{{font-size:13px}} .kv b{{display:block;color:var(--mut);font-weight:500;font-size:11px}}
@media(max-width:700px){{body{{padding:14px}}}}
</style></head><body>'''


def page(title, body):
    return HEAD.format(t=ESC(title)) + body + '</body></html>'


def internal_sampling(a, allc):
    rng = random.Random(a.seed)
    chosen = defaultdict(list)
    for lab in a.random_rounds:
        pool = [c for c in allc if c['run'] == lab]
        take = pool if len(pool) <= a.n_random else rng.sample(pool, a.n_random)
        for c in take:
            chosen[c['uid']].append('random')
        print(f'random stratum {lab}: {len(take)} of {len(pool)} (denominator {len(pool)})')
    enriched = []
    rest = [c for c in allc if 'random' not in chosen.get(c['uid'], [])]
    buckets = {
        'clear_structure': [c for c in rest if c['AUTO_architecture'] in
                            ('regular_spaced_array', 'variable_spacer_array', 'tandem_repeat')],
        'many_copies': [c for c in rest if (c['AUTO_n_copies'] or '0').isdigit()
                        and int(c['AUTO_n_copies'] or 0) >= 6],
        'high_no_model': [c for c in rest if c['AUTO_explainer_status'] == 'no_model'],
        'short_period': [c for c in rest if c['AUTO_architecture'] == 'short_period_tandem'],
        'repeat_pair': [c for c in rest if c['AUTO_architecture'] == 'repeat_pair'],
    }
    per = max(1, a.n_enriched // max(1, len(buckets)))
    seen_lin = defaultdict(int)
    for name, pool in buckets.items():
        pool = sorted(pool, key=lambda c: -c['AUTO_score'])
        picked = 0
        for c in pool:
            if picked >= per:
                break
            if seen_lin[c['contig']] >= 2:      # avoid one genome dominating
                continue
            chosen[c['uid']].append(f'enriched:{name}')
            seen_lin[c['contig']] += 1
            enriched.append(c)
            picked += 1
    return chosen


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs', nargs='+', required=True,
                   help='label=path, e.g. r1_host=/path/region_round1')
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--n-random', type=int, default=20)
    p.add_argument('--n-enriched', type=int, default=20)
    p.add_argument('--seed', type=int, default=20261016)
    p.add_argument('--random-rounds', nargs='*', default=['r2_host', 'r2_phage'])
    p.add_argument('--existing-review', type=Path)
    p.add_argument('--selection', type=Path,
                   help='browse_selection.tsv from select_browse_book.py. When given, '
                        'it is the authoritative browse-book membership and the '
                        'internal sampling below is not used, so the sampling list is '
                        'frozen in exactly one place.')
    a = p.parse_args()
    for d in ('report', 'manifests', 'tables', 'review', 'candidates', 'figures',
              'sequences', 'provenance'):
        (a.out / d).mkdir(parents=True, exist_ok=True)

    runs = {}
    for spec in a.runs:
        lab, path = spec.split('=', 1)
        runs[lab] = Path(path)

    # ---------------------------------------------------------- gather
    allc, tracks_index = [], {}
    for lab, R in runs.items():
        route = 'phage' if 'phage' in lab else 'host'
        # label is rN_route; the round comes from that prefix, never from a
        # two-way guess, which silently mislabelled round-3 as round-2.
        if not (len(lab) > 1 and lab[0] == 'r' and lab[1].isdigit()):
            raise ValueError(f'run label must look like r<N>_<route>, got {lab!r}')
        rnd = f'round-{lab[1]}'
        cands = tsv(R / 'calls/audit_full.tsv')
        models = {m['candidate_id']: m for m in tsv(R / 'explain_full/models.tsv')}
        probe = {m['candidate_id']: m for m in tsv(R / 'no_model_probe/no_model_probe.tsv')}
        copies = defaultdict(list)
        for c in tsv(R / 'explain_full/copies.tsv'):
            copies[c['candidate_id']].append(c)
        hyp = defaultdict(list)
        for h in tsv(R / 'explain_full/hypotheses.tsv'):
            hyp[h['candidate_id']].append(h)
        members = defaultdict(list)
        for m in tsv(R / 'calls/members.tsv'):
            members[m.get('candidate_id', '')].append(m)
        peaks = defaultdict(list)
        for pk in tsv(R / 'calls/peaks.tsv'):
            peaks[pk['candidate_id']].append(pk)
        ngenes = defaultdict(list)
        for gn in tsv(R / 'neighborhoods/genes.tsv'):
            ngenes[gn['candidate_id']].append(gn)
        trna = {x['candidate_id']: x for x in tsv(R / 'calls/trna_overlap.tsv')}
        crispr = {x['candidate_id']: x for x in tsv(R / 'calls/crispr_overlap.tsv')}
        for c in cands:
            cid = f"{c['contig']}:{c['start0']}-{c['end0']}"
            m = models.get(cid, {})
            pr = probe.get(cid, {})
            allc.append({
                'uid': f'{lab}__{cid}'.replace(':', '_'), 'run': lab, 'round': rnd,
                'route': route, 'candidate_id': cid, 'contig': c['contig'],
                'start0': int(c['start0']), 'end0': int(c['end0']),
                'span_nt': int(c['end0']) - int(c['start0']),
                'AUTO_score': round(float(c['score']), 4),
                'AUTO_n_peaks': c.get('n_peaks', ''),
                'AUTO_explainer_status': m.get('status', 'no_model'),
                'AUTO_architecture': m.get('architecture', ''),
                'AUTO_consensus': m.get('consensus', ''),
                'AUTO_unit_len': m.get('unit_len', ''),
                'AUTO_n_copies': m.get('n_copies', '0'),
                'AUTO_n_hypotheses': len(hyp.get(cid, [])),
                'AUTO_probe_verdict': pr.get('verdict', ''),
                'AUTO_source_windows': c.get('source_windows', ''),
                'AUTO_n_peaks_exported': len({(x['peak_start0'], x['peak_end0'])
                                               for x in peaks.get(cid, [])}),
                'AUTO_n_neighbor_CDS': len(ngenes.get(cid, [])),
                '_copies': copies.get(cid, []), '_hyp': hyp.get(cid, []),
                '_members': members.get(cid, []), '_run_path': R,
                'AUTO_n_tRNA': (trna.get(cid, {}) or {}).get('n_tRNA_overlapping', '0'),
                'AUTO_tRNA_names': (trna.get(cid, {}) or {}).get('tRNA_names', ''),
                'AUTO_n_CRISPR_arrays': (crispr.get(cid, {}) or {}).get(
                    'n_CRISPR_arrays_overlapping', '0'),
                '_trna': trna.get(cid, {}), '_crispr': crispr.get(cid, {}),
                '_peaks': peaks.get(cid, []), '_ngenes': ngenes.get(cid, []),
                '_model_row': m})
        tracks_index[lab] = str(R / 'scores')

    # ---------------------------------------------------------- sampling
    chosen = defaultdict(list)
    if a.selection:
        by_cid = {(c['run'], c['candidate_id']): c for c in allc}
        missing = []
        for r in tsv(a.selection):
            k = (r['run'], r['candidate_id'])
            c = by_cid.get(k)
            if c is None:
                missing.append(k)
                continue
            strata = ['random'] if r.get('in_random_sample', '').lower() == 'true' else []
            strata += [f'enriched:{x}' for x in r['rules'].split(',')
                       if x and not x.startswith('random')]
            chosen[c['uid']] = strata
        if missing:
            raise ValueError(f'selection lists candidates absent from the runs: {missing[:5]}'
                             f' ({len(missing)} total). Never silently drop a selected row.')
        book = [c for c in allc if c['uid'] in chosen]
        n_rand = sum(1 for c in book if 'random' in chosen[c['uid']])
        print(f'browse book from frozen selection {a.selection}: {len(book)} candidates '
              f'({n_rand} random, {len(book) - n_rand} enriched-only)')
    else:
        chosen = internal_sampling(a, allc)
        book = [c for c in allc if c['uid'] in chosen]
        n_rand = sum(1 for c in book if 'random' in chosen[c['uid']])
        print(f'browse book: {len(book)} candidates ({n_rand} random, '
              f'{len(book) - n_rand} enriched)')
    # ---------------------------------------------------------- tables
    cols = [k for k in allc[0] if not k.startswith('_')]
    with (a.out / 'tables/all_candidates.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, cols + ['in_browse_book', 'strata'], delimiter='\t',
                           extrasaction='ignore')
        w.writeheader()
        for c in allc:
            w.writerow({**{k: c[k] for k in cols},
                        'in_browse_book': c['uid'] in chosen,
                        'strata': ','.join(chosen.get(c['uid'], []))})
    review_cols = ['uid', 'round', 'route', 'candidate_id', 'strata', 'AUTO_score',
                   'AUTO_explainer_status', 'AUTO_architecture', 'AUTO_n_copies',
                   'AUTO_probe_verdict', 'HUMAN_sequence_evidence',
                   'HUMAN_explanation_coverage', 'HUMAN_biological_context', 'HUMAN_notes']
    with (a.out / 'review/browse_book_review.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, review_cols, delimiter='\t')
        w.writeheader()
        for c in book:
            w.writerow({**{k: c.get(k, '') for k in review_cols},
                        'strata': ','.join(chosen[c['uid']]),
                        'HUMAN_sequence_evidence': '', 'HUMAN_explanation_coverage': '',
                        'HUMAN_biological_context': '', 'HUMAN_notes': ''})
    if a.existing_review and a.existing_review.exists():
        import shutil
        shutil.copy(a.existing_review, a.out / 'review' / a.existing_review.name)

    (a.out / 'provenance/track_archive_index.json').write_text(json.dumps({
        'note': 'full per-base GPU tracks are NOT inside this package; they are kept '
                'in the companion archive at the paths below',
        'paths': tracks_index}, indent=2) + '\n')
    # ---------------------------------------------------------- per-candidate pages
    import numpy as np
    genomes, gffs, trackcache = {}, {}, {}

    def genome_of(run_path):
        if run_path not in genomes:
            genomes[run_path] = read_fasta(run_path / 'input.fa')
        return genomes[run_path]

    def gff_of(run_path):
        if run_path not in gffs:
            d = defaultdict(list)
            f = run_path / 'cds/genes.gff3'
            if f.exists():
                for line in f.open():
                    if line.startswith('#'):
                        continue
                    t = line.rstrip('\n').split('\t')
                    if len(t) > 8 and t[2] == 'CDS':
                        m = re.search(r'ID=([^;]+)', t[8])
                        d[t[0]].append((int(t[3]) - 1, int(t[4]), t[6],
                                        m.group(1) if m else 'CDS'))
            gffs[run_path] = d
        return gffs[run_path]

    def track_for(run_path, contig, lo, hi):
        """Orientation-averaged track over [lo,hi) assembled from window NPZs."""
        key = str(run_path)
        if key not in trackcache:
            idx = defaultdict(list)
            for rj in sorted(run_path.glob('scores/**/run.json')):
                meta = json.loads(rj.read_text())
                for rec in meta['records']:
                    idx[rec['id']].append((rj.parent / rec['track'],
                                           rec['scoring_orientation']))
            wins = {w['window_id']: w for w in tsv(run_path / 'prepared/windows.tsv')}
            trackcache[key] = (idx, wins)
        idx, wins = trackcache[key]
        acc = np.full(hi - lo, np.nan)
        for wid, w in wins.items():
            if w['contig'] != contig:
                continue
            ws, we = int(w['start0']), int(w['end0'])
            if we <= lo or ws >= hi:
                continue
            vals = []
            for path, _o in idx.get(wid, []):
                try:
                    vals.append(np.load(path)['above_uniform'].astype(float))
                except Exception:
                    pass
            if not vals:
                continue
            v = np.mean(vals, axis=0)
            s0, e0 = max(ws, lo), min(we, hi)
            seg = v[s0 - ws:e0 - ws]
            tgt = acc[s0 - lo:e0 - lo]
            acc[s0 - lo:e0 - lo] = np.where(np.isnan(tgt), seg, (tgt + seg) / 2)
        return acc

    manifest_rows = {}
    for src in (Path.home() / 'repeat_db_v1/manifests/manifest.tsv',
                Path.home() / 'repeat_db_v1/manifests/download_ledger_round2.tsv'):
        for r in tsv(src):
            manifest_rows[r['accession_version']] = r

    for c in book:
        R = c['_run_path']
        g = genome_of(R)[c['contig']]
        cs, ce = c['start0'], c['end0']
        pad = max(1500, (ce - cs) * 2)
        lo, hi = max(0, cs - pad), min(len(g), ce + pad)
        cps = [(int(x['start0']), int(x['end0'])) for x in c['_copies']]
        genes = [x for x in gff_of(R).get(c['contig'], []) if x[0] < hi and lo < x[1]]
        tr = track_for(R, c['contig'], lo, hi)
        vals = [0.0 if np.isnan(x) else float(x) for x in tr]
        pk_spans = sorted({(int(x['peak_start0']), int(x['peak_end0']))
                           for x in c['_peaks']})
        fig1 = svg_track(vals, cand=(cs, ce), copies=cps, peaks=pk_spans, x0=lo,
                         title=f"detector track {c['candidate_id']}")
        zlo, zhi = max(0, cs - 300), min(len(g), ce + 300)
        tz = track_for(R, c['contig'], zlo, zhi)
        vz = [0.0 if np.isnan(x) else float(x) for x in tz]
        fig2 = svg_track(vz, cand=(cs, ce), copies=cps, peaks=pk_spans, x0=zlo,
                         height=96, title='zoom')
        fig3 = svg_genes(genes, lo, hi, cand=(cs, ce))
        (a.out / f'figures/{c["uid"]}_track.svg').write_text(fig1)
        (a.out / f'figures/{c["uid"]}_genes.svg').write_text(fig3)
        seq = g[cs:ce]
        wrap = lambda t: '\n'.join(t[i:i + 70] for i in range(0, len(t), 70))
        (a.out / f'sequences/{c["uid"]}.fa').write_text(
            f">{c['candidate_id']}|{c['round']}|{c['route']}\n" + wrap(seq) + '\n')
        ctx_pad = 500
        xs, xe = max(0, cs - ctx_pad), min(len(g), ce + ctx_pad)
        ctx = g[xs:xe]
        (a.out / f'sequences/{c["uid"]}_context.fa').write_text(
            f">{c['contig']}:{xs}-{xe}|context|candidate_at_{cs}-{ce}\n"
            + wrap(ctx) + '\n')
        mrow = manifest_rows.get(c['contig'], {})
        tax = mrow.get('taxonomy', mrow.get('family', ''))
        ncds = len(genes)
        ncop = len(cps)
        summary = (f"候选区间 {ce - cs:,} nt；" +
                   (f"发现 {ncop} 个相似 copy；" if ncop else "当前方法未找到 copy；") +
                   f"邻域含 {ncds} 个预测 CDS。")
        # ---- copies grouped under the model that proposed them -------------
        hyp_by_id = {h['model_id']: h for h in c['_hyp']}
        cop_by_model = defaultdict(list)
        for x in c['_copies']:
            cop_by_model[x.get('model_id', '')].append(x)
        primary_id = (c['_model_row'] or {}).get('model_id', '')
        order = sorted(cop_by_model, key=lambda m: (m != primary_id, m))
        blocks = []
        for mid in order:
            h = hyp_by_id.get(mid, {})
            cl = sorted(cop_by_model[mid], key=lambda x: int(x['start0']))
            is_p = mid == primary_id
            rows = ''.join(
                f"<tr><td>{i}</td><td><code>{x['contig']}:{x['start0']}-{x['end0']}</code></td>"
                f"<td>{int(x['end0']) - int(x['start0'])}</td>"
                f"<td>{x.get('mismatches_to_consensus','')}</td>"
                f"<td>{x.get('consensus_identity','')}</td>"
                f"<td>{x.get('loo_consensus_identity','')}</td>"
                f"<td>{x.get('overlaps_candidate','')}</td>"
                f"<td>{x.get('contained_in_candidate','')}</td>"
                f"<td class='seq'>{ESC(x.get('sequence',''))}</td></tr>"
                for i, x in enumerate(cl, 1))
            spac = h.get('median_start_to_start', '') or '—'
            blocks.append(f'''<div class="card">
<h3>{'主展示模型 · ' if is_p else '替代假设 · '}<code>{ESC(mid)}</code>
{'<span class="tag">primary</span>' if is_p else ''}
<span class="tag">{ESC(h.get('status','') or '—')}</span></h3>
<div class="grid">
<div class="kv"><b>结构</b>{ESC(h.get('architecture','') or '—')}</div>
<div class="kv"><b>来源</b>{ESC(h.get('model_source','') or '—')}</div>
<div class="kv"><b>consensus</b><code>{ESC(h.get('consensus','') or '—')}</code></div>
<div class="kv"><b>单元长度</b>{ESC(h.get('unit_len','') or '—')}</div>
<div class="kv"><b>copy 数</b>{ESC(h.get('n_copies','') or str(len(cl)))}</div>
<div class="kv"><b>落在候选内</b>{ESC(h.get('n_copies_contained_in_candidate','') or '—')}</div>
<div class="kv"><b>中位 identity</b>{ESC(h.get('median_consensus_identity','') or '—')}</div>
<div class="kv"><b>中位留一法 identity</b>{ESC(h.get('median_loo_consensus_identity','') or '—')}</div>
<div class="kv"><b>起点间距中位数</b>{ESC(str(spac))}</div>
<div class="kv"><b>间距 CV</b>{ESC(h.get('spacing_cv','') or '—')}</div>
<div class="kv"><b>周期</b>{ESC(h.get('period','') or '—')}</div>
<div class="kv"><b>相位</b>{ESC(h.get('phase_status','') or '—')}</div>
</div>
<table><thead><tr><th>#</th><th>坐标</th><th>长度</th><th>mismatch</th>
<th>identity</th><th>留一法 identity</th><th>与候选相交</th><th>完全在候选内</th>
<th>序列（完整）</th></tr></thead><tbody>{rows or '<tr><td colspan=9>该模型无 copy 记录</td></tr>'}</tbody></table>
</div>''')
        models_html = ''.join(blocks) or '<p class="mut">当前方法未给出任何重复模型。</p>'
        no_cop = [h for h in c['_hyp'] if h['model_id'] not in cop_by_model]
        if no_cop:
            models_html += ('<p class="mut">另有 ' + str(len(no_cop)) +
                            ' 个假设没有 copy 记录：' +
                            ', '.join(f'<code>{ESC(h["model_id"])}</code>' for h in no_cop[:12]) +
                            ('…' if len(no_cop) > 12 else '') + '</p>')

        # ---- peaks actually called by the frozen detector -------------------
        pk_u = {}
        for x in c['_peaks']:
            k = (int(x['peak_start0']), int(x['peak_end0']))
            pk_u.setdefault(k, []).append(x)
        rows_pk = ''.join(
            f"<tr><td>{i}</td><td><code>{k[0]:,}-{k[1]:,}</code></td><td>{k[1]-k[0]}</td>"
            f"<td>{int(v[0]['summit0']):,}</td>"
            f"<td>{max(float(y['peak_mean_above_uniform']) for y in v):.4f}</td>"
            f"<td>{'、'.join(sorted({y['orientation'] for y in v}))}</td>"
            f"<td>{v[0]['contained_in_candidate']}</td>"
            f"<td class='mut'>{ESC('、'.join(sorted({y['window_id'] for y in v})))}</td></tr>"
            for i, (k, v) in enumerate(sorted(pk_u.items()), 1))
        gaps_pk = [b[0] - a_[1] for a_, b in zip(sorted(pk_u)[:-1], sorted(pk_u)[1:])]

        # ---- existing annotation in the neighbourhood ----------------------
        rows_gn = ''.join(
            f"<tr><td><code>{ESC(gn['gene_id'])}</code></td>"
            f"<td><code>{int(gn['start0']):,}-{int(gn['end0']):,}</code></td>"
            f"<td>{ESC(gn['strand'])}</td><td>{ESC(gn['relation'])}</td>"
            f"<td>{ESC(gn.get('distance_nt','') or '0')}</td>"
            f"<td>{ESC(gn.get('overlap_nt','') or '0')}</td>"
            f"<td>{ESC(gn.get('product',''))}</td>"
            f"<td>{ESC(gn.get('protein_available',''))}</td></tr>"
            for gn in sorted(c['_ngenes'], key=lambda g: int(g['start0'])))
        rows_hyp = ''
        mrow_m = c['_model_row'] or {}
        tr_, cr_ = c['_trna'] or {}, c['_crispr'] or {}
        n_tr = int(tr_.get('n_tRNA_overlapping', 0) or 0)
        n_cr = int(cr_.get('n_CRISPR_arrays_overlapping', 0) or 0)
        ann_bits = []
        if n_tr:
            ann_bits.append(
                f'<div class="card"><b>tRNA</b>：与 {n_tr} 个注释 tRNA 相交，覆盖候选的 '
                f"{float(tr_.get('tRNA_coverage_fraction', 0) or 0):.0%}。"
                f"<br>名称：<code>{ESC(tr_.get('tRNA_names', ''))}</code>"
                f"<br>坐标：<code>{ESC(tr_.get('tRNA_coords', ''))}</code>"
                '<br><span class="mut">tRNA 是高度受约束的结构性基因，不是序列重复。'
                'detector 在这里打分高，并不意味着存在 repeat。</span></div>')
        if n_cr:
            ann_bits.append(
                f'<div class="card"><b>CRISPR</b>：位于 pilercr 注释的 {n_cr} 个阵列内，'
                f"覆盖候选的 {float(cr_.get('CRISPR_coverage_fraction', 0) or 0):.0%}。"
                f"<br>拷贝数 <code>{ESC(cr_.get('CRISPR_array_copies', ''))}</code>，"
                f"repeat <code>{ESC(cr_.get('CRISPR_repeat_len', ''))}</code> nt，"
                f"spacer <code>{ESC(cr_.get('CRISPR_spacer_len', ''))}</code> nt"
                f"<br>阵列 consensus：<code>{ESC(cr_.get('CRISPR_consensus', ''))}</code>"
                '<br><span class="mut">explainer 的 12-nt consensus 通常只是完整 CRISPR '
                'repeat 的一个片段，因此它描述了真实阵列，但不等于其真正的重复单元。</span></div>')
        if cr_.get('CRISPR_repeat_without_called_array') == 'True':
            ann_bits.append(
                '<div class="warn"><b>携带 CRISPR repeat 序列，但不在任何 pilercr 阵列内。</b>'
                f"<br>匹配到的阵列 consensus：<code>{ESC(cr_.get('model_consensus_in_CRISPR_repeat', ''))}</code>"
                '<br>可能是退化或残缺的阵列，专用工具未能识别；值得单独核查。</div>')
        if not ann_bits:
            ann_bits.append('<p class="mut">未与任何注释 tRNA 或 CRISPR 阵列相交。</p>')
        ann_html = ''.join(ann_bits)
        rows_win = ''.join(
            f"<tr><td><code>{ESC(m.get('window_id',''))}</code></td>"
            f"<td>{m.get('orientation','')}</td><td>{m.get('score','')}</td></tr>"
            for m in c['_members'][:40])
        unexplained = ''
        if c['AUTO_explainer_status'] == 'no_model' or not cps:
            unexplained = ('<div class="warn"><b>当前方法未解释该候选。</b>'
                           'explainer 未给出 repeat 模型。下方提供原始候选序列与上下文，'
                           'probe 结论仅为补充证据：probe 阳性不等于可信 repeat，'
                           'probe 阴性不等于该区间没有重复（短单元灵敏度仅约 0.22）。</div>')
        body = f"""<p><a href="../index.html">← 候选索引</a> · <a href="../START_HERE.html">首页</a></p>
<h1>{ESC(c['candidate_id'])}</h1>
<p><span class="tag">{c['round']}</span><span class="tag">{c['route']}</span>
<span class="tag">{ESC(','.join(chosen[c['uid']]))}</span></p>
<div class="warn">计算候选，尚未人工确认。本页所有自动标签以 AUTO_ 标注。</div>
<p>{summary}</p>
<div class="card"><div class="grid">
<div class="kv"><b>accession.version</b>{ESC(c['contig'])}</div>
<div class="kv"><b>分类</b>{ESC(tax or '—')}</div>
<div class="kv"><b>坐标 (0-based, half-open)</b>{cs:,}–{ce:,}</div>
<div class="kv"><b>AUTO detector 分数</b>{c['AUTO_score']}</div>
<div class="kv"><b>AUTO peaks</b>{ESC(str(c['AUTO_n_peaks']))}</div>
<div class="kv"><b>AUTO explainer 状态</b>{ESC(c['AUTO_explainer_status'])}</div>
<div class="kv"><b>AUTO 结构</b>{ESC(c['AUTO_architecture'] or '—')}</div>
<div class="kv"><b>AUTO copy 数</b>{ESC(str(c['AUTO_n_copies']))}</div>
<div class="kv"><b>AUTO 假设数</b>{c['AUTO_n_hypotheses']}</div>
<div class="kv"><b>AUTO probe</b>{ESC(c['AUTO_probe_verdict'] or '未运行')}</div>
</div></div>
{unexplained}
<h2>detector 信号</h2>
<p class="mut">纵轴为 masked predictability（nats above uniform），上限 ln(4)=1.386。
双方向分数在同一窗口内平均；跨窗口仅在重叠处平均，未做跨窗口任意平滑。</p>
{fig1}
<h3>候选范围放大</h3>{fig2}
<h2>基因邻域</h2>{fig3}
<h2>detector peak（冻结阈值 τ=0.80 实际调用结果）</h2>
<p class="mut">这些就是 detector 自己调用的 peak，按窗口与方向原样复算，再映射到基因组坐标；
同一 peak 可能由多个窗口或两个方向各报一次，表中已按坐标去重并记录来源。
peak 的均值上限为 ln(4)=1.386。</p>
<table><thead><tr><th>#</th><th>坐标</th><th>长度</th><th>summit</th>
<th>峰均值 (nats)</th><th>方向</th><th>完全在候选内</th><th>来源窗口</th></tr></thead>
<tbody>{rows_pk or '<tr><td colspan=8>无</td></tr>'}</tbody></table>
<p class="mut">共 {len(pk_u)} 个去重 peak{f'；相邻 peak 边缘间隔：{", ".join(str(x) for x in gaps_pk[:20])}{"…" if len(gaps_pk) > 20 else ""} nt' if gaps_pk else ''}。
这里报告的是<b>边缘间隔</b>，不是起点到起点的 pitch，两者不可直接与注释周期比较。</p>

<h2>repeat 模型（按模型分组展示 copies）</h2>
<p class="mut">每个模型单独成块，copies 列在提出它们的模型之下，<b>序列完整未截断</b>。
主展示模型只是展示选择，不代表它比替代假设更正确；全部假设都保留。
未知基因保留未知，不根据邻近位置推断功能。</p>
{models_html}

<h2>已有注释：非编码 RNA 与 CRISPR</h2>
<p class="mut">pyrodigal 只预测蛋白编码基因，结构性非编码 RNA 与 CRISPR 阵列在其中<b>不可见</b>。
这里补上 aragorn（tRNA）与 pilercr（CRISPR 阵列）的注释。
这是<b>注释</b>，不是判断：它说明 detector 打在了什么已知结构上。</p>
{ann_html}

<h2>已有注释（邻域内预测 CDS）</h2>
<p class="mut">基因为 pyrodigal 从头预测，<b>没有功能注释</b>——所有 product 都是
"hypothetical protein"。因此下表只能说明位置与朝向关系，不含任何功能判断。
AUTO CDS 重叠数 {ESC(str(mrow_m.get('n_CDS_overlapping','') or '—'))}，
重叠比例 {ESC(str(mrow_m.get('CDS_overlap_fraction','') or '—'))}。</p>
<table><thead><tr><th>gene_id</th><th>坐标</th><th>链</th><th>relation</th>
<th>距离 nt</th><th>重叠 nt</th><th>product</th><th>有蛋白序列</th></tr></thead>
<tbody>{rows_gn or '<tr><td colspan=8>邻域内无预测 CDS</td></tr>'}</tbody></table>

<h3>原始窗口成员</h3>
<table><thead><tr><th>window_id</th><th>方向</th><th>分数</th></tr></thead>
<tbody>{rows_win or '<tr><td colspan=3>—</td></tr>'}</tbody></table>

<h2>序列</h2>
<h3>候选区间（{ce - cs:,} nt，完整）</h3>
<div class="seq">{ESC(seq)}</div>
<h3>上下文（候选 ±{ctx_pad} nt，{xe - xs:,} nt；候选区间以 [ ] 标出）</h3>
<div class="seq">{ESC(ctx[:cs - xs])}<b>[{ESC(ctx[cs - xs:ce - xs])}]</b>{ESC(ctx[ce - xs:])}</div>
<p><a href="../sequences/{c['uid']}.fa">下载候选 FASTA</a> ·
<a href="../sequences/{c['uid']}_context.fa">下载上下文 FASTA</a> ·
<a href="../figures/{c['uid']}_track.svg">track SVG</a> ·
<a href="../figures/{c['uid']}_genes.svg">邻域 SVG</a> ·
<a href="../tables/all_candidates.tsv">总表 TSV</a> ·
<a href="../review/browse_book_review.tsv">审查表 TSV</a></p>
<h2>人工审查</h2>
<p class="mut">请在 <code>review/browse_book_review.tsv</code> 中填写；此处仅为字段说明。</p>
<table><tbody>
<tr><th>HUMAN_sequence_evidence</th><td>clear_repeat / suspected_repeat / no_repeat_found / cannot_judge</td></tr>
<tr><th>HUMAN_explanation_coverage</th><td>model_adequate / partial / boundary_or_phase_problem / structure_unsupported</td></tr>
<tr><th>HUMAN_biological_context</th><td>known_common_structure / worth_followup / function_unknown</td></tr>
<tr><th>HUMAN_notes</th><td>自由文本</td></tr>
</tbody></table>
<p class="mut">局限：环状基因组原点跨越未检测；VIPR 已知漏检；probe 对短单元灵敏度低。</p>"""
        (a.out / f'candidates/{c["uid"]}.html').write_text(page(c['candidate_id'], body))

    # ---------------------------------------------------------- index
    rows_idx = ''.join(
        f"<tr><td><a href='candidates/{c['uid']}.html'>{ESC(c['candidate_id'])}</a></td>"
        f"<td>{c['round']}</td><td>{c['route']}</td><td>{ESC(c['contig'])}</td>"
        f"<td>{c['AUTO_score']}</td><td>{ESC(c['AUTO_explainer_status'])}</td>"
        f"<td>{ESC(c['AUTO_architecture'] or '')}</td><td>{ESC(str(c['AUTO_n_copies']))}</td>"
        f"<td>{ESC(','.join(chosen.get(c['uid'], [])) or '')}</td></tr>"
        for c in sorted(book, key=lambda x: -x['AUTO_score']))
    idx_body = f"""<p><a href="START_HERE.html">← 首页</a></p>
<h1>浏览册候选索引（{len(book)}）</h1>
<div class="warn">计算候选，尚未人工确认。</div>
<p><input id="q" placeholder="搜索 accession / 结构 / 状态" size="34">
<select id="r"><option value="">全部轮次</option><option>round-1</option><option>round-2</option><option>round-3</option></select>
<select id="t"><option value="">全部路线</option><option>host</option><option>phage</option></select>
<select id="s"><option value="">全部状态</option><option>candidate_associated</option>
<option>flank_only</option><option>no_model</option></select></p>
<table id="tb"><thead><tr><th>候选</th><th>轮次</th><th>路线</th><th>基因组</th>
<th>分数</th><th>explainer</th><th>结构</th><th>copies</th><th>抽样</th></tr></thead>
<tbody>{rows_idx}</tbody></table>
<p class="mut">全部 {len(allc)} 个候选见 <a href="tables/all_candidates.tsv">tables/all_candidates.tsv</a>。</p>
<script>
const q=document.getElementById('q'),r=document.getElementById('r'),t=document.getElementById('t'),s=document.getElementById('s');
function f(){{const Q=q.value.toLowerCase();document.querySelectorAll('#tb tbody tr').forEach(tr=>{{
const x=tr.innerText.toLowerCase();const c=tr.children;
tr.style.display=(x.includes(Q)&&(!r.value||c[1].innerText==r.value)&&(!t.value||c[2].innerText==t.value)&&(!s.value||c[5].innerText==s.value))?'':'none';}});}}
[q,r,t,s].forEach(e=>e.addEventListener('input',f));
document.querySelectorAll('#tb th').forEach((th,i)=>th.addEventListener('click',()=>{{
const tb=document.querySelector('#tb tbody');const rows=[...tb.rows];const n=i==4||i==7;
const asc=th.dataset.a!=='1';th.dataset.a=asc?'1':'0';
rows.sort((x,y)=>{{const A=x.cells[i].innerText,B=y.cells[i].innerText;
return (n?(parseFloat(A)||0)-(parseFloat(B)||0):A.localeCompare(B))*(asc?1:-1);}});
rows.forEach(z=>tb.appendChild(z));}}));
</script>"""
    (a.out / 'index.html').write_text(page('候选索引', idx_body))
    print(f'wrote {len(book)} candidate pages + index')
    return a, runs, allc, book, chosen


if __name__ == '__main__':
    main()
