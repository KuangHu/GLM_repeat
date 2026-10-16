"""CPU-only engineering tests with dummy models, NOT gLM2 biological validation."""
import importlib.util
import math
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
import torch

path = Path(__file__).parents[1] / 'scripts' / 'reference_score.py'
spec = importlib.util.spec_from_file_location('reference_score', path)
r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r)


class Tokenizer:
    # Deliberately not the real gLM2 IDs: the adapter must look IDs up.
    vocab = {'<+>': 1, 'a': 4, 'c': 7, 'g': 8, 't': 11}
    mask_token_id = 12
    unk_token_id = 0

    def convert_tokens_to_ids(self, token):
        return self.vocab.get(token, 0)

    def __call__(self, sequence, **kwargs):
        assert sequence.startswith('<+>')
        assert kwargs['truncation'] is False
        assert kwargs['add_special_tokens'] is False
        ids = torch.tensor([[1] + [self.vocab[c] for c in sequence[3:]]])
        return {'input_ids': ids, 'attention_mask': torch.ones_like(ids)}


class UniformModel(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.calls = []

    def forward(self, input_ids, attention_mask, return_dict):
        assert return_dict
        assert (attention_mask == 1).all()
        self.calls.extend(input_ids.clone())
        logits = torch.zeros((*input_ids.shape, 13))
        logits[..., 0] = 100  # Non-DNA logit must not affect the four-base softmax.
        return SimpleNamespace(logits=logits)


def test_normalization():
    assert r.normalize_dna('ACGT') == 'acgt'


@pytest.mark.parametrize('s', ['', 'acngt', 'acugt', 'ac-gt', 'ac gt'])
def test_ambiguous_and_empty_rejected(s):
    with pytest.raises(ValueError):
        r.normalize_dna(s)


def test_reverse_complement():
    assert r.reverse_complement('AACG') == 'cgtt'
    assert r.reverse_complement(r.reverse_complement('AACG')) == 'aacg'


def test_exact_plan():
    assert [p.tolist() for p in r.mask_plan(4, 'exact')] == [[0], [1], [2], [3]]


def test_random_plan_coverage_and_determinism():
    a = r.mask_plan(99, 'random_partition', groups=20, repeats=2, seed=21)
    b = r.mask_plan(99, 'random_partition', groups=20, repeats=2, seed=21)
    assert len(a) == 40
    assert all(np.array_equal(x, y) for x, y in zip(a, b))
    assert np.array_equal(np.bincount(np.concatenate(a)), np.full(99, 2))


def test_short_random_plan():
    a = r.mask_plan(3, 'random_partition', groups=20, repeats=2)
    assert len(a) == 6
    assert np.array_equal(np.bincount(np.concatenate(a)), np.full(3, 2))


def test_mapping():
    ids, _, pos, bases = r.tokenize_dna(Tokenizer(), 'ACGT')
    assert pos.tolist() == [1, 2, 3, 4]
    assert ids[0, pos].tolist() == bases


def test_overflow_rejected():
    with pytest.raises(ValueError):
        r.tokenize_dna(Tokenizer(), 'ACGT', max_tokens=4)


def test_uniform_scores_and_actual_masking():
    model = UniformModel()
    out = r.score_sequence(model, Tokenizer(), 'acgtacg', batch_size=3)
    assert np.allclose(out['log_p4'], -math.log(4), atol=1e-6)
    assert np.allclose(out['above_uniform'], 0, atol=1e-6)
    assert (out['n_evaluations'] == 1).all()
    assert out['n_variants'] == 7
    assert len(model.calls) == 7
    for i, row in enumerate(model.calls):
        assert row[0] == 1
        assert (row == Tokenizer.mask_token_id).sum() == 1
        assert row[i + 1] == Tokenizer.mask_token_id


def test_random_scores():
    model = UniformModel()
    out = r.score_sequence(model, Tokenizer(), 'acgt' * 25,
        mode='random_partition', groups=20, repeats=2, batch_size=3)
    assert out['n_variants'] == 40
    assert (out['n_evaluations'] == 2).all()
    assert np.allclose(out['above_uniform'], 0, atol=1e-6)


def test_batch_invariance():
    a = r.score_sequence(UniformModel(), Tokenizer(), 'acgtacgt', batch_size=1)
    b = r.score_sequence(UniformModel(), Tokenizer(), 'acgtacgt', batch_size=5)
    assert np.array_equal(a['log_p4'], b['log_p4'])


def test_fasta(tmp_path):
    p = tmp_path / 't.fa'
    p.write_text('>one description\nACGT\nAC\n>two\nTTGC\n')
    assert list(r.read_fasta(p)) == [('one', 'acgtac'), ('two', 'ttgc')]
    p.write_text('>same\nAC\n>same\nGT\n')
    with pytest.raises(ValueError):
        list(r.read_fasta(p))
