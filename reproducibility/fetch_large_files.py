#!/usr/bin/env python3
"""Download the paper's frozen derived inputs without authentication and verify SHA-256."""
import hashlib, json, os, subprocess, tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()
def main():
    items = json.loads((ROOT / 'reproducibility/public_inputs.json').read_text())
    for item in items:
        target = ROOT / item['path']
        if target.exists() and digest(target) == item['sha256']:
            print('VERIFIED', item['path'], flush=True)
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        fd, name = tempfile.mkstemp(prefix='.download-', dir=target.parent)
        os.close(fd)
        temporary = Path(name)
        try:
            subprocess.run(['curl', '--fail', '--location', '--retry', '3', '--silent', '--show-error', '--output', name, item['url']], check=True)
            if temporary.stat().st_size != item['bytes'] or digest(temporary) != item['sha256']:
                raise RuntimeError('Input checksum/size mismatch: ' + item['path'])
            temporary.replace(target)
            print('DOWNLOADED AND VERIFIED', item['path'], flush=True)
        finally:
            temporary.unlink(missing_ok=True)
    print('All public inputs verified.')
if __name__ == '__main__':
    main()
