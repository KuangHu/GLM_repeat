#!/usr/bin/env python3
"""Validate the offline review package before it is handed to a reviewer.

Checks that would otherwise only surface while someone is reading the package:

  sequences    every candidate FASTA matches the genome at its stated coordinates
  consistency  copy counts agree between the copies table, the candidate page and
               the figure; browse-book membership agrees between the selection
               manifest, the review sheet and the full table
  offline      no absolute cluster path, no network reference, no missing local
               target among the links and figures the pages reference
  figures      every referenced SVG exists and carries real plotted geometry
  labels       AUTO_/HUMAN_ discipline: HUMAN_ columns blank, AUTO_ labels present
  checksums    MANIFEST.sha256 verifies against the files on disk

Exit status is non-zero if any check fails, so this cannot pass silently.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

FAIL, WARN = [], []


def fail(msg):
    FAIL.append(msg)
    print(f'  FAIL  {msg}')


def warn(msg):
    WARN.append(msg)
    print(f'  WARN  {msg}')


def ok(msg):
    print(f'  ok    {msg}')


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
    if not p.exists():
        return []
    with p.open() as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--package', type=Path, required=True)
    ap.add_argument('--run', action='append', required=True, metavar='LABEL=DIR')
    ap.add_argument('--selection', type=Path, required=True)
    a = ap.parse_args()
    P = a.package
    runs = dict(s.partition('=')[::2] for s in a.run)

    print('\n[1] candidate sequences vs genome coordinates')
    genomes = {lab: read_fasta(Path(d) / 'input.fa') for lab, d in runs.items()}
    n_seq = n_bad = 0
    for fa in sorted((P / 'sequences').glob('*.fa')):
        uid = fa.stem
        lab = uid.split('__')[0]
        rec = read_fasta(fa)
        if not rec:
            fail(f'{fa.name}: empty FASTA')
            continue
        header = next(iter(rec))
        seq = rec[header]
        m = re.search(r'([A-Za-z0-9_.]+):(\d+)-(\d+)', header)
        if not m:
            fail(f'{fa.name}: header carries no parsable coordinates: {header}')
            continue
        contig, s0, e0 = m.group(1), int(m.group(2)), int(m.group(3))
        g = genomes.get(lab, {}).get(contig)
        n_seq += 1
        if g is None:
            fail(f'{fa.name}: contig {contig} absent from {lab} input.fa')
            n_bad += 1
        elif g[s0:e0] != seq:
            fail(f'{fa.name}: sequence != genome[{contig}:{s0}-{e0}] '
                 f'(len {len(seq)} vs {e0-s0})')
            n_bad += 1
    ok(f'{n_seq} candidate FASTAs checked against the genome, {n_bad} mismatched')

    print('\n[2] copy counts: copies table vs page vs figure')
    # The page renders EVERY copy row the explainer emitted for the candidate
    # (all hypotheses), while AUTO_n_copies is the primary model's own count.
    # These are different numbers by design, so they are checked separately.
    all_copies, primary_copies = {}, {}
    for lab, d in runs.items():
        pm = {m['candidate_id']: m['model_id']
              for m in tsv(Path(d) / 'explain_full/models.tsv')
              if m.get('primary', '').lower() in ('true', '1', 'yes')}
        a_by, p_by = defaultdict(int), defaultdict(int)
        for c in tsv(Path(d) / 'explain_full/copies.tsv'):
            a_by[c['candidate_id']] += 1
            if pm.get(c['candidate_id']) == c['model_id']:
                p_by[c['candidate_id']] += 1
        all_copies[lab], primary_copies[lab] = a_by, p_by
    book = tsv(P / 'review/browse_book_review.tsv')
    n_chk = n_mis = 0
    for r in book:
        uid = r['uid']
        lab = uid.split('__')[0]
        cid = r['candidate_id']
        pg = P / f'candidates/{uid}.html'
        if not pg.exists():
            fail(f'{uid}: candidate page missing')
            continue
        html = pg.read_text()
        declared = r.get('AUTO_n_copies', '') or '0'
        table_rows = len(re.findall(
            r'<tr><td>\d+</td><td><code>[A-Za-z0-9_.]+:\d+-\d+</code></td>', html))
        n_all_cop = all_copies[lab].get(cid, 0)
        n_prim = primary_copies[lab].get(cid, 0)
        n_chk += 1
        if declared.isdigit() and n_prim and int(declared) != n_prim:
            fail(f'{uid}: AUTO_n_copies={declared} but the primary model has '
                 f'{n_prim} copy rows')
            n_mis += 1
        if table_rows != n_all_cop:
            fail(f'{uid}: page renders {table_rows} copy rows but the explainer '
                 f'emitted {n_all_cop} for this candidate')
            n_mis += 1
        fig = P / f'figures/{uid}_track.svg'
        if n_all_cop and fig.exists():
            bands = len(re.findall(r'class="copy"', fig.read_text()))
            if bands == 0 and 'copy' not in fig.read_text():
                warn(f'{uid}: {n_all_cop} copies but the track SVG draws no copy band')
    ok(f'{n_chk} browse-book pages cross-checked, {n_mis} inconsistent')

    print('\n[3] browse-book membership agrees across manifest / sheet / table')
    sel = tsv(a.selection)
    sel_keys = {(r['run'], r['candidate_id']) for r in sel}
    sheet_keys = {(r['uid'].split('__')[0], r['candidate_id']) for r in book}
    full = tsv(P / 'tables/all_candidates.tsv')
    tbl_keys = {(r['run'], r['candidate_id']) for r in full
                if r['in_browse_book'].lower() == 'true'}
    if sel_keys != sheet_keys:
        fail(f'selection manifest vs review sheet differ by '
             f'{len(sel_keys ^ sheet_keys)} candidates')
    if sel_keys != tbl_keys:
        fail(f'selection manifest vs full table in_browse_book differ by '
             f'{len(sel_keys ^ tbl_keys)} candidates')
    if sel_keys == sheet_keys == tbl_keys:
        ok(f'{len(sel_keys)} browse-book candidates agree across all three')
    n_rand_sel = sum(1 for r in sel if r['in_random_sample'].lower() == 'true')
    n_rand_sheet = sum(1 for r in book if 'random' in r['strata'].split(','))
    if n_rand_sel != n_rand_sheet:
        fail(f'random stratum size differs: manifest {n_rand_sel}, sheet {n_rand_sheet}')
    else:
        ok(f'random stratum consistently {n_rand_sel} candidates')
    prov = json.loads((P / 'manifests/selection_provenance.json').read_text())
    listed = sum(len(v) for v in prov['random_sampling_lists'].values())
    if listed != n_rand_sel:
        fail(f'saved sampling list has {listed} ids but {n_rand_sel} rows are flagged random')
    else:
        ok(f'saved sampling list matches the random flag ({listed} ids)')

    print('\n[4] offline integrity: links, absolute paths, network references')
    pages = sorted(P.rglob('*.html'))
    abs_path = net = broken = 0
    for pg in pages:
        t = pg.read_text()
        for m in re.findall(r'(?:href|src)="([^"]+)"', t):
            if m.startswith(('http://', 'https://', '//')):
                net += 1
                fail(f'{pg.relative_to(P)}: network reference {m}')
            elif m.startswith('#') or m.startswith('data:'):
                continue
            elif m.startswith('/'):
                abs_path += 1
                fail(f'{pg.relative_to(P)}: absolute path {m}')
            else:
                tgt = (pg.parent / m).resolve()
                if not tgt.exists():
                    broken += 1
                    fail(f'{pg.relative_to(P)}: broken link {m}')
        for m in re.findall(r'/global/(?:scratch|home)/[^\s"\'<>]+', t):
            abs_path += 1
            fail(f'{pg.relative_to(P)}: cluster path leaked: {m[:80]}')
    ok(f'{len(pages)} pages scanned: {net} network refs, {abs_path} absolute paths, '
       f'{broken} broken links')

    print('\n[5] figures carry real geometry')
    empty = 0
    figs = sorted((P / 'figures').glob('*.svg'))
    for f in figs:
        t = f.read_text()
        pts = re.findall(r'points="([^"]*)"', t) + re.findall(r'\sd="([^"]*)"', t)
        shapes = len(re.findall(r'<(?:rect|polyline|path|polygon|circle|line)\b', t))
        coords = sum(len(p.split()) for p in pts)
        if shapes < 2 or (coords < 10 and '<rect' not in t):
            empty += 1
            fail(f'figures/{f.name}: no plotted geometry ({shapes} shapes, {coords} points)')
    ok(f'{len(figs)} SVGs checked, {empty} empty')
    referenced = set()
    for pg in pages:
        for m in re.findall(r'(?:href|src)="\.\./figures/([^"]+)"', pg.read_text()):
            referenced.add(m)
    missing = referenced - {f.name for f in figs}
    if missing:
        fail(f'{len(missing)} referenced figures do not exist: {sorted(missing)[:3]}')
    else:
        ok(f'{len(referenced)} referenced figures all present')

    print('\n[6] AUTO_/HUMAN_ discipline')
    human_cols = [c for c in (book[0] if book else {}) if c.startswith('HUMAN_')]
    nonblank = [(r['uid'], c) for r in book for c in human_cols if r[c].strip()]
    if not human_cols:
        fail('review sheet has no HUMAN_ columns')
    elif nonblank:
        fail(f'{len(nonblank)} HUMAN_ cells are pre-filled, e.g. {nonblank[:3]}')
    else:
        ok(f'{len(human_cols)} HUMAN_ columns present and blank in all {len(book)} rows')
    auto_cols = [c for c in (full[0] if full else {}) if c.startswith('AUTO_')]
    ok(f'{len(auto_cols)} AUTO_ columns in the full table')
    no_notice = [pg.name for pg in (P / 'candidates').glob('*.html')
                 if '尚未人工确认' not in pg.read_text()]
    if no_notice:
        fail(f'{len(no_notice)} candidate pages lack the unconfirmed notice')
    else:
        ok('every candidate page carries the "尚未人工确认" notice')

    print('\n[7] peaks, context sequences and model grouping on the pages')
    peaks_by = {}
    for lab, d in runs.items():
        by = defaultdict(set)
        for r in tsv(Path(d) / 'calls/peaks.tsv'):
            by[r['candidate_id']].add((r['peak_start0'], r['peak_end0']))
        peaks_by[lab] = by
    n_pk = n_ctx = n_trunc = n_grp = 0
    for r in book:
        uid, cid = r['uid'], r['candidate_id']
        lab = uid.split('__')[0]
        html = (P / f'candidates/{uid}.html').read_text()
        want = len(peaks_by[lab].get(cid, ()))
        got = len(re.findall(
            r'<tr><td>\d+</td><td><code>[\d,]+-[\d,]+</code></td>', html))
        if want != got:
            fail(f'{uid}: {want} deduplicated peaks exported but the page lists {got}')
            n_pk += 1
        ctx = P / f'sequences/{uid}_context.fa'
        if not ctx.exists():
            fail(f'{uid}: context FASTA missing')
            n_ctx += 1
        else:
            rec = read_fasta(ctx)
            hdr = next(iter(rec))
            m = re.match(r'([A-Za-z0-9_.]+):(\d+)-(\d+)', hdr)
            g = genomes[lab].get(m.group(1), '')
            if g[int(m.group(2)):int(m.group(3))] != rec[hdr]:
                fail(f'{uid}: context FASTA does not match the genome')
                n_ctx += 1
        if re.search(r'候选区间[^<]*完整[^<]*</h3>\s*<div class="seq">[^<]*…', html):
            fail(f'{uid}: candidate sequence is still truncated on the page')
            n_trunc += 1
        mids = {x['model_id'] for x in tsv(Path(runs[lab]) / 'explain_full/copies.tsv')
                if x['candidate_id'] == cid}
        blocks = html.count('主展示模型 ·') + html.count('替代假设 ·')
        if blocks != len(mids):
            fail(f'{uid}: {len(mids)} models have copies but the page renders '
                 f'{blocks} model blocks')
            n_grp += 1
    ok(f'{len(book)} pages: peaks {len(book)-n_pk} consistent, context seqs '
       f'{len(book)-n_ctx} verified, {len(book)-n_trunc} untruncated, '
       f'{len(book)-n_grp} correctly grouped by model')

    print('\n[8] MANIFEST.sha256')
    man = P / 'MANIFEST.sha256'
    listed_files, bad = 0, 0
    recorded = {}
    for line in man.read_text().splitlines():
        h, _, rel = line.partition('  ')
        recorded[rel] = h
        listed_files += 1
    on_disk = {str(f.relative_to(P)) for f in P.rglob('*')
               if f.is_file() and f.name != 'MANIFEST.sha256'}
    if set(recorded) != on_disk:
        fail(f'manifest lists {len(recorded)} files, disk has {len(on_disk)}; '
             f'symmetric difference {len(set(recorded) ^ on_disk)}')
    for rel, h in recorded.items():
        f = P / rel
        if not f.exists():
            bad += 1
            continue
        d = hashlib.sha256()
        with f.open('rb') as fh:
            for b in iter(lambda: fh.read(1 << 20), b''):
                d.update(b)
        if d.hexdigest() != h:
            bad += 1
            fail(f'checksum mismatch: {rel}')
    if bad == 0 and set(recorded) == on_disk:
        ok(f'all {listed_files} files verify against MANIFEST.sha256')

    print(f'\n{"=" * 64}')
    print(f'FAIL {len(FAIL)}   WARN {len(WARN)}')
    if FAIL:
        print('\nPackage is NOT ready. Failures:')
        for m in FAIL[:40]:
            print(f'  - {m}')
    sys.exit(1 if FAIL else 0)


if __name__ == '__main__':
    main()
