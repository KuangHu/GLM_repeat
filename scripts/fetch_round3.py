#!/usr/bin/env python3
"""Download round-3 genomes and record a QC ledger.

Every sequence is checked against the length recorded in the selection manifest
before it is accepted, and its sha256 is written to the ledger, so the scanned
bytes can always be tied back to the frozen selection. Sequences that fail are
logged and excluded rather than silently truncated or retried into the panel.

The GenBank/FASTA cross-check from earlier rounds taught that an over-eager
sanitiser can delete IUPAC ambiguity letters; here only whitespace is stripped
and any non-ACGT content is measured, not removed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import time
import urllib.parse
import urllib.request
from pathlib import Path

E = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'


def efetch(accs):
    q = urllib.parse.urlencode({'db': 'nuccore', 'id': ','.join(accs),
                                'rettype': 'fasta', 'retmode': 'text',
                                'tool': 'GLM_repeat_round3'})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(f'{E}/efetch.fcgi?{q}', timeout=300) as fh:
                return fh.read().decode()
        except Exception:
            if attempt == 4:
                raise
            time.sleep(3 * (attempt + 1))


def parse(text):
    out, name, buf = {}, None, []
    for line in text.splitlines():
        if line.startswith('>'):
            if name:
                out[name] = ''.join(buf)
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(''.join(line.split()))     # strip whitespace ONLY
    if name:
        out[name] = ''.join(buf)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest', action='append', required=True, metavar='ROUTE=FILE')
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--batch', type=int, default=50)
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)
    ledger, failures = [], []

    for spec in a.manifest:
        route, _, f = spec.partition('=')
        rows = list(csv.DictReader(Path(f).open(), delimiter='\t'))
        want = {r['accession_version']: int(r['length_nt']) for r in rows}
        meta = {r['accession_version']: r for r in rows}
        got = {}
        accs = list(want)
        for i in range(0, len(accs), a.batch):
            chunk = accs[i:i + a.batch]
            got.update(parse(efetch(chunk)))
            print(f'  {route}: {len(got)}/{len(accs)}', flush=True)
            time.sleep(0.4)
        out = a.outdir / f'{route}.fa'
        with out.open('w') as fh:
            for acc in accs:
                seq = got.get(acc, '').upper()
                if not seq:
                    failures.append((route, acc, 'not returned by efetch'))
                    continue
                if len(seq) != want[acc]:
                    failures.append((route, acc,
                                     f'length {len(seq)} != manifest {want[acc]}'))
                    continue
                n_bad = sum(1 for c in seq if c not in 'ACGT')
                fh.write(f'>{acc}\n')
                for j in range(0, len(seq), 70):
                    fh.write(seq[j:j + 70] + '\n')
                m = meta[acc]
                ledger.append({
                    'accession_version': acc, 'route': route,
                    'length_nt': len(seq),
                    'gc_fraction': round((seq.count('G') + seq.count('C')) / len(seq), 5),
                    'non_acgt_fraction': round(n_bad / len(seq), 8),
                    'sequence_sha256': hashlib.sha256(seq.encode()).hexdigest(),
                    'length_matches_manifest': True,
                    'organism': m.get('organism', m.get('virus_name', '')),
                    'family': m.get('family', ''), 'phylum': m.get('phylum', ''),
                    'taxid': m.get('taxid', m.get('virus_taxid', '')),
                })
        print(f'{route}: wrote {out} ({sum(1 for l in ledger if l["route"]==route)} records)')

    with (a.outdir / 'download_ledger_round3.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, list(ledger[0]), delimiter='\t')
        w.writeheader()
        w.writerows(ledger)
    if failures:
        with (a.outdir / 'download_failures_round3.log').open('w') as fh:
            for r, acc, why in failures:
                fh.write(f'{r}\t{acc}\t{why}\n')
    tot = {}
    for l in ledger:
        tot[l['route']] = tot.get(l['route'], 0) + l['length_nt']
    for r, n in sorted(tot.items()):
        print(f'{r}: {n:,} nt ({n/1e6:.2f} Mb)')
    print(f'TOTAL {sum(tot.values())/1e6:.2f} Mb; {len(failures)} failures')


if __name__ == '__main__':
    main()
