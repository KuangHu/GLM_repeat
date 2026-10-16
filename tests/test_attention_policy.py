"""Attention-mask policy: auto drops only an all-ones mask; padding is preserved."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import torch

path = Path(__file__).parents[1] / 'scripts' / 'reference_score.py'
spec = importlib.util.spec_from_file_location('reference_score', path)
r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r)

VOCAB = {'<+>': 1, 'a': 4, 'c': 7, 'g': 8, 't': 11}
PAD, MASK = 2, 12


class Tokenizer:
    mask_token_id = MASK
    unk_token_id = 0

    def __init__(self, pad_tail=0):
        self.pad_tail = pad_tail

    def convert_tokens_to_ids(self, token):
        return VOCAB.get(token, 0)

    def __call__(self, sequence, **kwargs):
        ids = [1] + [VOCAB[c] for c in sequence[3:]] + [PAD] * self.pad_tail
        attn = [1] * (len(ids) - self.pad_tail) + [0] * self.pad_tail
        return {'input_ids': torch.tensor([ids]), 'attention_mask': torch.tensor([attn])}


class RecordingModel(torch.nn.Module):
    """Logits depend on whether a mask was passed, so a dropped meaningful mask
    would change the scores, not just the call log."""

    def __init__(self):
        super().__init__()
        self.masks = []

    def forward(self, input_ids, attention_mask, return_dict):
        self.masks.append(None if attention_mask is None else attention_mask.clone())
        logits = torch.zeros((*input_ids.shape, 13))
        logits[..., VOCAB['a']] = 1.0
        if attention_mask is not None and (attention_mask == 0).any():
            logits[..., VOCAB['a']] = 2.0  # padding-aware path
        return SimpleNamespace(logits=logits)


def test_resolve_all_ones():
    ones = torch.ones(3, 9, dtype=torch.long)
    assert r.resolve_attention_mask(ones, 'auto') is None
    assert r.resolve_attention_mask(ones, 'always') is ones


def test_resolve_preserves_padding():
    m = torch.ones(2, 9, dtype=torch.long)
    m[1, -2:] = 0
    assert r.resolve_attention_mask(m, 'auto') is m
    assert r.resolve_attention_mask(m, 'always') is m


def test_resolve_rejects_unknown_policy():
    with pytest.raises(ValueError):
        r.resolve_attention_mask(torch.ones(1, 3), 'never')
    with pytest.raises(ValueError):
        r.score_sequence(RecordingModel(), Tokenizer(), 'acgt', attention_mask_policy='x')


def test_default_policy_is_conservative_until_validated():
    assert r.DEFAULT_ATTENTION_MASK_POLICY in r.ATTENTION_MASK_POLICIES


@pytest.mark.parametrize('policy', ['auto', 'always'])
def test_padded_input_keeps_mask(policy):
    model = RecordingModel()
    out = r.score_sequence(model, Tokenizer(pad_tail=2), 'acgta', batch_size=2,
                           attention_mask_policy=policy)
    assert all(m is not None and (m[:, -2:] == 0).all() for m in model.masks)
    ref = r.score_sequence(RecordingModel(), Tokenizer(pad_tail=2), 'acgta',
                           batch_size=1, attention_mask_policy='always')
    assert np.array_equal(out['log_p4'], ref['log_p4'])


def test_unpadded_auto_drops_mask_always_keeps_it():
    auto, always = RecordingModel(), RecordingModel()
    a = r.score_sequence(auto, Tokenizer(), 'acgtac', batch_size=4, attention_mask_policy='auto')
    b = r.score_sequence(always, Tokenizer(), 'acgtac', batch_size=4,
                         attention_mask_policy='always')
    assert auto.masks and all(m is None for m in auto.masks)
    assert always.masks and all(m is not None and (m == 1).all() for m in always.masks)
    assert np.array_equal(a['log_p4'], b['log_p4'])


@pytest.mark.parametrize('policy', ['auto', 'always'])
def test_mlm_masking_unchanged_by_policy(policy):
    seen = []

    class Spy(RecordingModel):
        def forward(self, input_ids, attention_mask, return_dict):
            seen.extend(input_ids.clone())
            return super().forward(input_ids, attention_mask, return_dict)

    r.score_sequence(Spy(), Tokenizer(), 'acgtacg', batch_size=3, attention_mask_policy=policy)
    assert len(seen) == 7
    for i, row in enumerate(seen):
        assert (row == MASK).sum() == 1 and row[i + 1] == MASK
