"""Download the published fields, verify SHA256, and restore repository paths.

Python 3.9+ standard library only. Different existing files are never overwritten.
"""
import hashlib
import json
from pathlib import Path
import shutil
import urllib.request
import zipfile


ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    manifest = json.loads((ROOT / 'results/ni_al_lookup/release_manifest.json').read_text())
    archive = manifest['archive']
    entries = manifest['files']
    targets = {}
    # Validate paths and existing files before downloading or extracting anything.
    for entry in entries:
        target = (ROOT / entry['path']).resolve()
        target.relative_to((ROOT / 'results').resolve())
        if target.exists() and sha256(target) != entry['sha256']:
            raise SystemExit(f'Different local file; keep it or move it aside before retrying: {target}')
        targets[entry['path']] = target
    if all(path.exists() for path in targets.values()):
        print('All published files are already present and verified.')
        return
    cache = ROOT / 'downloads' / archive['name']
    cache.parent.mkdir(exist_ok=True)
    if not cache.exists() or sha256(cache) != archive['sha256']:
        partial = cache.with_suffix('.zip.part')
        print(f"Downloading {archive['size_bytes']:,} bytes from {archive['url']}", flush=True)
        request = urllib.request.Request(archive['url'], headers={'User-Agent': 'MOOSE-Calc-results'})
        with urllib.request.urlopen(request, timeout=90) as response, partial.open('wb') as output:
            shutil.copyfileobj(response, output)
        if sha256(partial) != archive['sha256']:
            raise SystemExit(f'Archive SHA256 mismatch; incomplete file retained at {partial}')
        partial.replace(cache)
    with zipfile.ZipFile(cache) as source:
        if sorted(source.namelist()) != sorted(targets):
            raise SystemExit('Archive member list differs from the publication manifest.')
        for entry in entries:
            target = targets[entry['path']]
            if target.exists():
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            temporary = target.with_name(target.name + '.download-part')
            with source.open(entry['path']) as incoming, temporary.open('wb') as output:
                shutil.copyfileobj(incoming, output)
            if temporary.stat().st_size != entry['size_bytes'] or sha256(temporary) != entry['sha256']:
                raise SystemExit(f'Extracted file checksum mismatch: {temporary}')
            temporary.replace(target)
            print(f"Restored {entry['path']}", flush=True)
    print('Published results restored and verified.')


if __name__ == '__main__':
    main()
