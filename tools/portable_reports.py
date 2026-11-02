"""Convert generated local report links and JSON paths to repository paths."""
import json
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def convert(value):
    if isinstance(value, dict):
        return {convert(key): convert(item) for key, item in value.items()}
    if isinstance(value, list):
        return [convert(item) for item in value]
    if isinstance(value, str):
        normalized = value.replace('\\', '/')
        prefix = ROOT.as_posix() + '/'
        if normalized.startswith(prefix):
            return normalized[len(prefix):]
    return value

def main():
    changed = 0
    pattern = re.compile(re.escape(ROOT.as_posix()) + r'/([^)>\n]+)')
    for path in (ROOT/'results').rglob('*'):
        if path.suffix not in {'.md', '.json'}:
            continue
        old = path.read_text(encoding='utf-8-sig')
        if path.suffix == '.json':
            obj = json.loads(old); normalized = convert(obj)
            if obj == normalized:
                continue
            new = json.dumps(normalized, ensure_ascii=False, indent=2)+'\n'
        else:
            new = pattern.sub(lambda m: os.path.relpath(ROOT/m.group(1), path.parent).replace('\\','/'), old)
        if new != old:
            path.write_text(new, encoding='utf-8')
            changed += 1
    print(f'Converted {changed} report/metadata files to repository-relative paths.')

if __name__ == '__main__':
    main()
