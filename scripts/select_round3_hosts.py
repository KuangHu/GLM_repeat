#!/usr/bin/env python3
"""Select round-3 host genomes from phyla no previous round has scanned.

Every accession's lineage is VERIFIED against the phylum it was selected for
before it is accepted. An earlier round used four taxids without verification and
three were wrong (an ASFV eukaryotic virus, a mislabelled family and a bamboo), so
an unverified identifier is never allowed into a manifest here.

Selection uses only taxonomy, assembly completeness and sequence length. It never
uses detector scores, repeat content or neighbouring gene content, so the round-3
panel cannot be tuned toward what the detector already likes.
"""
from __future__ import annotations

import argparse
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

E = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
TOOL = 'GLM_repeat_round3'


def get(endpoint, **params):
    params.setdefault('tool', TOOL)
    url = f'{E}/{endpoint}.fcgi?' + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                return fh.read().decode()
        except Exception as exc:
            if attempt == 4:
                raise
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f'unreachable: {endpoint}')


def esearch(db, term, retmax=40):
    txt = get('esearch', db=db, term=term, retmax=retmax, retmode='json')
    return json.loads(txt)['esearchresult'].get('idlist', [])


def esummary(db, ids):
    if not ids:
        return {}
    txt = get('esummary', db=db, id=','.join(ids), retmode='json')
    d = json.loads(txt).get('result', {})
    return {k: v for k, v in d.items() if k != 'uids'}


def lineage(taxid):
    xml = get('efetch', db='taxonomy', id=str(taxid), retmode='xml')
    names = re.findall(r'<Lineage>(.*?)</Lineage>', xml, re.S)
    sci = re.findall(r'<ScientificName>(.*?)</ScientificName>', xml)
    return (names[0] if names else ''), (sci[0] if sci else '')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--phyla', required=True, help='comma-separated target phyla')
    ap.add_argument('--exclude', required=True, help='comma-separated accessions already used')
    ap.add_argument('--min-nt', type=int, default=1_000_000)
    ap.add_argument('--max-nt', type=int, default=4_500_000)
    ap.add_argument('--target-nt', type=int, default=20_000_000)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    used = {x.strip().split('.')[0] for x in a.exclude.split(',') if x.strip()}
    chosen, total, rejected = [], 0, []

    for phylum in [p.strip() for p in a.phyla.split(',') if p.strip()]:
        if total >= a.target_nt:
            break
        term = (f'"{phylum}"[Organism] AND srcdb_refseq[PROP] AND '
                f'"complete genome"[Title] AND {a.min_nt}:{a.max_nt}[SLEN] '
                f'AND biomol_genomic[PROP]')
        ids = esearch('nuccore', term, retmax=30)
        if not ids:
            rejected.append((phylum, 'no esearch hit'))
            continue
        summ = esummary('nuccore', ids[:20])
        picked = None
        # Curated RefSeq reference records (NC_) are preferred; NZ_ complete
        # chromosomes are accepted as a fallback because several phyla have no
        # NC_ record in the requested length band.
        ordered = sorted(summ.items(),
                         key=lambda kv: (not kv[1].get('accessionversion', '').startswith('NC_'),
                                         kv[1].get('accessionversion', '')))
        for uid, rec in ordered:
            acc = rec.get('accessionversion', '')
            if not acc.split('_')[0] in ('NC', 'NZ') or acc.split('.')[0] in used:
                continue
            slen = int(rec.get('slen', 0) or 0)
            if not (a.min_nt <= slen <= a.max_nt):
                continue
            taxid = rec.get('taxid')
            if not taxid:
                continue
            lin, sci = lineage(taxid)
            # MANDATORY: the organism's lineage must actually contain the phylum
            if phylum.lower() not in lin.lower():
                rejected.append((phylum, f'{acc} lineage mismatch: {sci} / {lin[:70]}'))
                continue
            picked = {'phylum': phylum, 'accession_version': acc, 'taxid': taxid,
                      'organism': sci, 'length_nt': slen, 'lineage': lin,
                      'title': rec.get('title', ''),
                      'selection_reason': f'phylum {phylum} not scanned in any prior round',
                      'lineage_verified': True}
            break
            time.sleep(0.4)
        if picked:
            chosen.append(picked)
            total += picked['length_nt']
            used.add(picked['accession_version'].split('.')[0])
            print(f"  + {picked['accession_version']:14s} {picked['length_nt']:>9,} nt  "
                  f"{picked['phylum']:22s} {picked['organism'][:44]}")
        else:
            rejected.append((phylum, 'no acceptable record'))
        time.sleep(0.5)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    cols = ['accession_version', 'taxid', 'organism', 'phylum', 'length_nt',
            'lineage_verified', 'selection_reason', 'lineage', 'title']
    with a.out.open('w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for c in chosen:
            fh.write('\t'.join(str(c[k]).replace('\t', ' ') for k in cols) + '\n')
    print(f'\n{len(chosen)} host genomes, {total:,} nt ({total/1e6:.2f} Mb) -> {a.out}')
    if rejected:
        print(f'{len(rejected)} phyla produced no accepted record:')
        for p, why in rejected:
            print(f'  - {p}: {why}')


if __name__ == '__main__':
    main()
