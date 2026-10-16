#!/usr/bin/env python3
"""Reference gLM2 masked scoring of pre-windowed DNA FASTA records.

This is a scoring primitive, NOT a trained repeat detector or a whole-genome
caller. Exact: one hidden nucleotide per variant. Random partition: an explicitly
approximate multi-mask statistic, with every nucleotide scored in every repeat.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import time
from pathlib import Path
from typing import Iterator
import numpy as np
import torch

DNA = 'acgt'


def normalize_dna(sequence: str) -> str:
    sequence = sequence.strip().lower()
    if not sequence:
        raise ValueError('Empty DNA sequence.')
    if set(sequence) - set(DNA):
        raise ValueError('Expected A/C/G/T only. Split at ambiguous bases upstream; do not delete them.')
    return sequence


def reverse_complement(sequence: str) -> str:
    return normalize_dna(sequence).translate(str.maketrans('acgt', 'tgca'))[::-1]


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    name, parts = None, []
    seen = set()
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield name, normalize_dna(''.join(parts))
                fields = line[1:].split()
                if not fields or fields[0] in seen:
                    raise ValueError('Empty or duplicated FASTA identifier.')
                name, parts = fields[0], []
                seen.add(name)
            elif name is None:
                raise ValueError('FASTA sequence before its first header.')
            else:
                parts.append(line)
    if name is not None:
        yield name, normalize_dna(''.join(parts))


def mask_plan(length: int, mode: str, groups: int = 20,
              repeats: int = 2, seed: int = 0) -> list[np.ndarray]:
    if length <= 0 or groups <= 0 or repeats <= 0:
        raise ValueError('length, groups and repeats must be positive.')
    if mode == 'exact':
        return [np.array([i], dtype=np.int64) for i in range(length)]
    if mode != 'random_partition':
        raise ValueError(f'Unknown mask mode: {mode}')
    rng = np.random.default_rng(seed)
    plan = []
    for _ in range(repeats):
        plan.extend(np.array_split(rng.permutation(length), min(groups, length)))
    return plan


def tokenize_dna(tokenizer, sequence: str, max_tokens: int = 4096):
    sequence = normalize_dna(sequence)
    enc = tokenizer('<+>' + sequence, add_special_tokens=False,
                    truncation=False, return_tensors='pt')
    ids = enc['input_ids']
    if ids.ndim != 2 or ids.shape[0] != 1 or ids.shape[1] > max_tokens:
        raise ValueError('Unexpected batch shape or context overflow; never silently truncate.')
    base_ids = [tokenizer.convert_tokens_to_ids(b) for b in DNA]
    if len(set(base_ids)) != 4 or tokenizer.unk_token_id in base_ids or None in base_ids:
        raise ValueError('Tokenizer does not expose four distinct lowercase nucleotide tokens.')
    positions = [i for i, token in enumerate(ids[0].tolist()) if token in base_ids]
    expected = [base_ids[DNA.index(b)] for b in sequence]
    if len(positions) != len(sequence) or ids[0, positions].tolist() != expected:
        raise ValueError('Nucleotide/token mapping failed; abort rather than shift coordinates.')
    mask_id = tokenizer.mask_token_id
    if mask_id is None or mask_id == tokenizer.unk_token_id or mask_id in base_ids:
        raise ValueError('No valid, distinct mask token.')
    attention = enc.get('attention_mask', torch.ones_like(ids))
    if attention.shape != ids.shape:
        raise ValueError('Invalid attention mask.')
    return ids, attention, np.asarray(positions, dtype=np.int64), base_ids


ATTENTION_MASK_POLICIES = ('auto', 'always')
DEFAULT_ATTENTION_MASK_POLICY = 'always'


def resolve_attention_mask(attention: torch.Tensor, policy: str):
    """Attention mask actually passed to the model.

    'always' passes the tokenizer mask unchanged. 'auto' omits it only when every
    entry is 1 (unpadded input), which lets SDPA use fused kernels; any mask that
    contains a 0 is preserved. MLM <mask> placement is unaffected.
    """
    if policy not in ATTENTION_MASK_POLICIES:
        raise ValueError(f'Unknown attention_mask_policy: {policy}')
    if policy == 'auto' and bool((attention == 1).all()):
        return None
    return attention


@torch.inference_mode()
def score_sequence(model, tokenizer, sequence: str, *, mode: str = 'exact',
                   groups: int = 20, repeats: int = 2, seed: int = 0,
                   batch_size: int = 1, device: str = 'cpu', max_tokens: int = 4096,
                   attention_mask_policy: str = DEFAULT_ATTENTION_MASK_POLICY):
    if batch_size <= 0:
        raise ValueError('batch_size must be positive.')
    if attention_mask_policy not in ATTENTION_MASK_POLICIES:
        raise ValueError(f'Unknown attention_mask_policy: {attention_mask_policy}')
    sequence = normalize_dna(sequence)
    original, attention, positions, base_ids = tokenize_dna(tokenizer, sequence, max_tokens)
    plan = mask_plan(len(sequence), mode, groups, repeats, seed)
    sums = np.zeros(len(sequence), dtype=np.float64)
    counts = np.zeros(len(sequence), dtype=np.int64)
    truth_index = np.array([DNA.index(b) for b in sequence])
    model.eval()
    for offset in range(0, len(plan), batch_size):
        batch_plan = plan[offset:offset + batch_size]
        ids = original.repeat(len(batch_plan), 1).to(device)
        attn = resolve_attention_mask(attention.repeat(len(batch_plan), 1).to(device),
                                      attention_mask_policy)
        for row, nucleotide_positions in enumerate(batch_plan):
            token_positions = torch.as_tensor(positions[nucleotide_positions], device=device)
            ids[row, token_positions] = tokenizer.mask_token_id
        out = model(input_ids=ids, attention_mask=attn, return_dict=True)
        logits = getattr(out, 'logits', None)
        if logits is None or logits.shape[:2] != ids.shape:
            raise ValueError('Need pretrained MLM logits, not encoder embeddings.')
        for row, nucleotide_positions in enumerate(batch_plan):
            token_positions = torch.as_tensor(positions[nucleotide_positions], device=device)
            selected = logits[row, token_positions][:, base_ids].float()
            log_p4 = torch.log_softmax(selected, dim=-1)
            target = torch.as_tensor(truth_index[nucleotide_positions], device=device)
            values = log_p4.gather(1, target[:, None]).squeeze(1).cpu().numpy()
            if not np.isfinite(values).all():
                raise FloatingPointError('Non-finite masked probabilities.')
            sums[nucleotide_positions] += values
            counts[nucleotide_positions] += 1
    expected_count = 1 if mode == 'exact' else repeats
    if not np.all(counts == expected_count):
        raise RuntimeError('Mask coverage is incomplete or duplicated unexpectedly.')
    logp = sums / counts
    return {
        'log_p4': logp.astype(np.float32),
        'above_uniform': (logp + math.log(4)).astype(np.float32),
        'n_evaluations': counts,
        'n_variants': len(plan),
    }


def load_model(snapshot: Path, device: str):
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    config_path = snapshot / 'config.json'
    if not (snapshot / 'snapshot_manifest.json').exists() or not config_path.exists():
        raise ValueError('Use a local snapshot created by download_model.py.')
    config = json.loads(config_path.read_text())
    if 'AutoModelForMaskedLM' not in config.get('auto_map', {}):
        raise ValueError('Snapshot lacks the AutoModelForMaskedLM mapping; inspect upstream config.')
    manifest = json.loads((snapshot / 'snapshot_manifest.json').read_text())
    for relative, expected in manifest.get('sha256', {}).items():
        path = (snapshot / relative).resolve()
        if not path.is_relative_to(snapshot.resolve()) or not path.is_file():
            raise ValueError('Snapshot file missing or path escaped snapshot directory.')
        digest = hashlib.sha256()
        with path.open('rb') as fh:
            for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b''):
                digest.update(chunk)
        if digest.hexdigest() != expected:
            raise ValueError(f'Snapshot integrity check failed: {relative}')
    dtype = torch.float32
    if device.startswith('cuda'):
        if not torch.cuda.is_available():
            raise RuntimeError('CUDA unavailable. Request a GPU compute node.')
        if not torch.cuda.is_bf16_supported():
            raise RuntimeError('This setup requests BF16; validate a different dtype explicitly.')
        dtype = torch.bfloat16
    tokenizer = AutoTokenizer.from_pretrained(str(snapshot), trust_remote_code=True,
                                              local_files_only=True)
    model, info = AutoModelForMaskedLM.from_pretrained(
        str(snapshot), trust_remote_code=True, local_files_only=True,
        use_safetensors=True, torch_dtype=dtype, output_loading_info=True)
    if info.get('missing_keys') or info.get('mismatched_keys') or info.get('error_msgs'):
        raise RuntimeError(f'Incomplete pretrained checkpoint load: {info}')
    if info.get('unexpected_keys'):
        raise RuntimeError(f'Unexpected checkpoint keys; review before scoring: {info}')
    model.requires_grad_(False)
    return model.to(device).eval(), tokenizer


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
    p.add_argument('--attention-mask-policy', choices=ATTENTION_MASK_POLICIES,
                   default=DEFAULT_ATTENTION_MASK_POLICY)
    p.add_argument('--trust-reviewed-code', action='store_true')
    a = p.parse_args()
    if not a.trust_reviewed_code:
        p.error('Review the downloaded model code, then pass --trust-reviewed-code.')
    if a.out.exists() and any(a.out.iterdir()):
        p.error('Output directory is nonempty; use a new directory.')
    a.out.mkdir(parents=True, exist_ok=True)
    model, tokenizer = load_model(a.snapshot.resolve(), a.device)
    manifest = json.loads((a.snapshot / 'snapshot_manifest.json').read_text())
    meta = {'model': manifest['model_id'], 'revision': manifest['revision'],
            'mode': a.mode, 'groups': a.groups, 'repeats': a.repeats,
            'seed': a.seed, 'both_strands': a.both_strands,
            'attention_mask_policy': a.attention_mask_policy,
            'dtype': str(next(model.parameters()).dtype),
            'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'score_interpretation': 'masked predictability, not repeat probability',
            'records': []}
    for name, sequence in read_fasta(a.fasta):
        orientations = [('+', sequence)]
        if a.both_strands:
            orientations.append(('-', reverse_complement(sequence)))
        for orientation, oriented_sequence in orientations:
            started = time.monotonic()
            result = score_sequence(model, tokenizer, oriented_sequence, mode=a.mode,
                groups=a.groups, repeats=a.repeats, seed=a.seed,
                batch_size=a.batch_size, device=a.device,
                attention_mask_policy=a.attention_mask_policy)
            # Both tracks are stored in coordinates of the input FASTA record.
            arrays = {k: v[::-1] if orientation == '-' else v for k, v in result.items()
                      if isinstance(v, np.ndarray)}
            key = hashlib.sha256((name + '\0' + orientation).encode()).hexdigest()[:24]
            outfile = a.out / (key + '.npz')
            np.savez_compressed(outfile, **arrays)
            meta['records'].append({'id': name, 'length': len(sequence),
                'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest(),
                'scoring_orientation': orientation, 'track': outfile.name,
                'coordinates': 'input_fasta_record_0_based',
                'n_variants': result['n_variants'], 'seconds': time.monotonic()-started})
    if not meta['records']:
        raise ValueError('FASTA contained no sequences.')
    (a.out / 'run.json').write_text(json.dumps(meta, indent=2) + '\n')
    print(a.out / 'run.json')


if __name__ == '__main__':
    main()
