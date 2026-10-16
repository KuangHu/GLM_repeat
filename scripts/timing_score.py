#!/usr/bin/env python3
"""Timed wrapper around reference_score.score_sequence (scoring logic unchanged).

Records model-load, per-orientation inference (CUDA-synchronized) and write time
separately, with per-orientation peak GPU memory. Tracks are written in the same
NPZ layout as reference_score.py.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reference_score as rs  # noqa: E402


def sync(device):
    if device.startswith('cuda'):
        torch.cuda.synchronize()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--snapshot', type=Path, required=True)
    p.add_argument('--fasta', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--device', default='cuda:0')
    p.add_argument('--mode', choices=['exact', 'random_partition'], default='exact')
    p.add_argument('--groups', type=int, default=20)
    p.add_argument('--repeats', type=int, default=2)
    p.add_argument('--batch-size', type=int, default=1)
    p.add_argument('--seed', type=int, default=20261014)
    p.add_argument('--both-strands', action='store_true')
    p.add_argument('--records', nargs='*', help='only score these FASTA ids')
    p.add_argument('--warmup-nt', type=int, default=3000,
                   help='length of real-sequence prefix used for warm-up passes')
    p.add_argument('--warmup-variants', type=int, default=8)
    p.add_argument('--attention-mask-policy', choices=rs.ATTENTION_MASK_POLICIES,
                   default=rs.DEFAULT_ATTENTION_MASK_POLICY)
    p.add_argument('--trust-reviewed-code', action='store_true')
    a = p.parse_args()
    if not a.trust_reviewed_code:
        p.error('Review the downloaded model code, then pass --trust-reviewed-code.')
    if a.out.exists() and any(a.out.iterdir()):
        p.error('Output directory is nonempty; use a new directory.')
    a.out.mkdir(parents=True, exist_ok=True)

    records = [(n, s) for n, s in rs.read_fasta(a.fasta) if not a.records or n in a.records]
    if not records:
        raise ValueError('No FASTA records selected.')

    sync(a.device)
    t0 = time.monotonic()
    model, tokenizer = rs.load_model(a.snapshot.resolve(), a.device)
    sync(a.device)
    load_s = time.monotonic() - t0

    # Warm-up: a few real-length forward passes, discarded.
    warm = records[0][1][:a.warmup_nt]
    ids, attn, positions, _ = rs.tokenize_dna(tokenizer, warm)
    with torch.inference_mode():
        for _ in range(max(1, a.warmup_variants // a.batch_size)):
            x = ids.repeat(a.batch_size, 1).to(a.device)
            x[:, positions[0]] = tokenizer.mask_token_id
            model(input_ids=x, return_dict=True, attention_mask=rs.resolve_attention_mask(
                attn.repeat(a.batch_size, 1).to(a.device), a.attention_mask_policy))
    sync(a.device)

    manifest = json.loads((a.snapshot / 'snapshot_manifest.json').read_text())
    meta = {'model': manifest['model_id'], 'revision': manifest['revision'],
            'mode': a.mode, 'groups': a.groups, 'repeats': a.repeats, 'seed': a.seed,
            'batch_size': a.batch_size, 'both_strands': a.both_strands,
            'attention_mask_policy': a.attention_mask_policy,
            'code_sha256': {f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                            for f in (Path(__file__).resolve(), Path(rs.__file__).resolve())},
            'device_name': torch.cuda.get_device_name(0) if a.device.startswith('cuda') else 'cpu',
            'dtype': str(next(model.parameters()).dtype),
            'model_load_seconds': load_s,
            'score_interpretation': 'masked predictability, not repeat probability',
            'records': []}
    for name, sequence in records:
        orientations = [('+', sequence)]
        if a.both_strands:
            orientations.append(('-', rs.reverse_complement(sequence)))
        for orientation, oriented in orientations:
            if a.device.startswith('cuda'):
                torch.cuda.reset_peak_memory_stats()
            sync(a.device)
            t0 = time.monotonic()
            result = rs.score_sequence(model, tokenizer, oriented, mode=a.mode,
                                       groups=a.groups, repeats=a.repeats, seed=a.seed,
                                       batch_size=a.batch_size, device=a.device,
                                       attention_mask_policy=a.attention_mask_policy)
            sync(a.device)
            infer_s = time.monotonic() - t0
            peak_mb = (torch.cuda.max_memory_allocated() / 2**20
                       if a.device.startswith('cuda') else None)
            t0 = time.monotonic()
            arrays = {k: v[::-1] if orientation == '-' else v for k, v in result.items()
                      if isinstance(v, np.ndarray)}
            key = hashlib.sha256((name + '\0' + orientation).encode()).hexdigest()[:24]
            outfile = a.out / (key + '.npz')
            np.savez_compressed(outfile, **arrays)
            write_s = time.monotonic() - t0
            meta['records'].append({
                'id': name, 'length': len(sequence),
                'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest(),
                'scoring_orientation': orientation, 'track': outfile.name,
                'coordinates': 'input_fasta_record_0_based',
                'n_variants': result['n_variants'],
                'inference_seconds': infer_s, 'write_seconds': write_s,
                'variants_per_second': result['n_variants'] / infer_s,
                'peak_gpu_mb': peak_mb})
            print(f"{name} {orientation} {len(sequence)}nt {result['n_variants']} variants "
                  f"{infer_s:.1f}s peak {peak_mb:.0f}MB", flush=True)
            (a.out / 'run.json').write_text(json.dumps(meta, indent=2) + '\n')
    print(a.out / 'run.json')


if __name__ == '__main__':
    main()
