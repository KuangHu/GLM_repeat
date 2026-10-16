#!/usr/bin/env python3
"""Download and pin a gLM2 snapshot without executing its remote Python code."""
import argparse
import hashlib
import json
import re
from pathlib import Path


def main():
    from huggingface_hub import HfApi, snapshot_download
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model', default='tattabio/gLM2_650M')
    p.add_argument('--revision', default='main')
    p.add_argument('--out', required=True)
    a = p.parse_args()
    out = Path(a.out).resolve()
    if (out / 'snapshot_manifest.json').exists():
        raise SystemExit('Existing pinned snapshot: use a different output directory.')
    info = HfApi().model_info(a.model, revision=a.revision)
    sha = info.sha
    if not sha or not re.fullmatch(r'[0-9a-f]{40}', sha):
        raise SystemExit('Cannot resolve a full commit SHA; refusing an unpinned download.')
    snapshot_download(repo_id=a.model, revision=sha, local_dir=str(out),
                      allow_patterns=['*.json', '*.py', '*.safetensors', 'README.md', 'LICENSE*'])
    hashes = {}
    for path in sorted(out.rglob('*')):
        if not path.is_file() or '.cache' in path.relative_to(out).parts:
            continue
        h = hashlib.sha256()
        with path.open('rb') as fh:
            for block in iter(lambda: fh.read(8 * 1024 * 1024), b''):
                h.update(block)
        hashes[str(path.relative_to(out))] = h.hexdigest()
    (out / 'snapshot_manifest.json').write_text(json.dumps({
        'model_id': a.model, 'revision': sha, 'sha256': hashes,
        'code_review_required': True,
    }, indent=2) + '\n')
    print(out)
    print('Review downloaded Python files before enabling trust_remote_code.')


if __name__ == '__main__':
    main()
