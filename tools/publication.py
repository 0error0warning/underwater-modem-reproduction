"""Resolve historical hashes through an explicit, path-sanitized public export."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def verified_source_hash(path, expected):
    path = Path(path).resolve()
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual == expected:
        return expected
    manifest = json.loads((ROOT / 'publication_manifest.json').read_text(encoding='utf-8'))
    item = manifest['files'].get(path.relative_to(ROOT).as_posix(), {})
    if item.get('source_sha256') == expected and item.get('sha256') == actual:
        return expected
    raise AssertionError(f'Unexplained source change: {path.relative_to(ROOT)}')

