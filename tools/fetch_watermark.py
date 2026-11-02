"""Fetch selected files from the official Watermark zip using HTTP ranges."""
import argparse
import hashlib
import io
import json
import pathlib
import urllib.request
import zipfile

URL = 'https://www.ffi.no/en/research/watermark/_/attachment/download/cc46828f-0078-46a1-bf77-8cee63b50056%3Af3b02654a06b2405fa042a5d0a05ce1b77305e40/WatermarkV1.zip'
ROOT = pathlib.Path(__file__).resolve().parents[1]


class RemoteFile(io.RawIOBase):
    def __init__(self):
        with urllib.request.urlopen(urllib.request.Request(URL, method='HEAD'), timeout=30) as response:
            self.size = int(response.headers['Content-Length'])
        self.position = 0

    def seekable(self):
        return True

    def seek(self, offset, whence=0):
        self.position = offset if whence == 0 else (self.position if whence == 1 else self.size) + offset
        return self.position

    def tell(self):
        return self.position

    def read(self, size=-1):
        size = min(self.size-self.position, size if size >= 0 else self.size)
        if size <= 0:
            return b''
        start = self.position
        request = urllib.request.Request(URL, headers={'Range': f'bytes={start}-{start+size-1}'})
        with urllib.request.urlopen(request, timeout=120) as response:
            if response.status != 206:
                raise RuntimeError(f'Server ignored range: {response.status}')
            actual = response.headers.get('Content-Range', '')
            if not actual.startswith(f'bytes {start}-{start+size-1}/'):
                raise RuntimeError(f'Unexpected range: {actual}')
            data = response.read()
        if len(data) != size:
            raise RuntimeError(f'Truncated response: {len(data)} != {size}')
        self.position += size
        return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fetch', nargs='*', help='Case-insensitive name substrings; omit for inventory only')
    args = parser.parse_args()
    destination = ROOT/'vendor'
    destination.mkdir(exist_ok=True)
    with zipfile.ZipFile(RemoteFile()) as archive:
        inventory = [{'name': entry.filename, 'bytes': entry.file_size, 'compressed_bytes': entry.compress_size} for entry in archive.infolist() if not entry.is_dir()]
        (ROOT/'references'/'watermark_archive_inventory.json').write_text(json.dumps({'source': URL, 'files': inventory}, indent=2), encoding='utf-8')
        if args.fetch is None:
            print(json.dumps(inventory, indent=2))
            return
        records = []
        for entry in archive.infolist():
            if entry.is_dir() or not any(pattern.lower() in entry.filename.lower() for pattern in args.fetch):
                continue
            target = (destination/entry.filename).resolve()
            if not target.is_relative_to(destination.resolve()):
                raise RuntimeError('Unsafe archive path')
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists() and target.stat().st_size == entry.file_size:
                print(f'Exists: {entry.filename}', flush=True)
                data = target.read_bytes()
            else:
                print(f'Fetching {entry.filename} ({entry.compress_size} compressed bytes)', flush=True)
                data = archive.read(entry)
                target.write_bytes(data)
            records.append({'name': entry.filename, 'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)})
        manifest = ROOT/'references'/'watermark_download_manifest.json'
        old = json.loads(manifest.read_text()) if manifest.exists() else {'source': URL, 'files': []}
        merged = {x['name']: x for x in old['files']}
        merged.update({x['name']: x for x in records})
        manifest.write_text(json.dumps({'source': URL, 'files': list(merged.values())}, indent=2), encoding='utf-8')
        print(f'Verified {len(records)} files (ZIP CRC and SHA256 recorded).', flush=True)


if __name__ == '__main__':
    main()
