#!/usr/bin/env python3
"""Fetch the selected runtime or upstream inputs and verify recorded SHA256s.

No sky realizations are downloaded. A system Python 3 is sufficient for this
bootstrap step; the numerical commands use the separately installed runtime.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import tempfile
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def safe_child(root, relative):
    target = (root / relative).resolve()
    target.relative_to(root.resolve())
    return target


def verify(path, record):
    if not path.is_file() or path.stat().st_size != record['bytes']:
        raise ValueError(f'Missing file or size mismatch: {path}')
    if digest(path) != record['sha256']:
        raise ValueError(f'SHA256 mismatch: {path}')


def acquire(record, target, args):
    print(f"{target.relative_to(ROOT)} ({record['bytes']:,} bytes)", flush=True)
    if args.plan:
        return
    if target.exists():
        verify(target, record)
        print('  verified existing file', flush=True)
        return
    if args.verify_only:
        raise FileNotFoundError(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    local = args.local_assets / record['filename'] if args.local_assets else None
    with tempfile.NamedTemporaryFile(dir=target.parent, prefix='.download-', delete=False) as stream:
        temporary = Path(stream.name)
        try:
            if local is not None:
                with local.open('rb') as source:
                    shutil.copyfileobj(source, stream, 8 * 1024 * 1024)
            else:
                request = Request(record['url'], headers={'User-Agent': 'cmb-sky-sim/1.1.0'})
                with urlopen(request, timeout=180) as source:
                    shutil.copyfileobj(source, stream, 8 * 1024 * 1024)
            stream.flush()
            verify(temporary, record)
            if target.exists():
                raise FileExistsError(f'Destination appeared during download: {target}')
            temporary.replace(target)
        finally:
            temporary.unlink(missing_ok=True)
    print('  SHA256 verified', flush=True)


def install_wheels(archive, args):
    manifest = json.loads((ROOT / 'cmb_sehgal2010/evidence/offline_wheels_manifest.json').read_text())
    expected = {record['path']: record for record in manifest['files']}
    destination = ROOT / 'cmb_sehgal2010/wheels'
    if args.plan:
        return
    if args.verify_only:
        for name, record in expected.items():
            verify(safe_child(destination, name), record)
        return
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as source:
        members = source.getmembers()
        if sorted(m.name for m in members) != sorted('wheels/' + name for name in expected):
            raise ValueError('Wheel archive inventory differs from the recorded inventory')
        if not all(m.isfile() for m in members):
            raise ValueError('Wheel archive must contain ordinary files only')
        with tempfile.TemporaryDirectory(dir=destination, prefix='.extract-') as temporary:
            stage = Path(temporary)
            for member in members:
                name = Path(member.name).name
                target = safe_child(stage, name)
                with source.extractfile(member) as src, target.open('wb') as out:
                    shutil.copyfileobj(src, out)
                verify(target, expected[name])
            for name, record in expected.items():
                target = safe_child(destination, name)
                if target.exists():
                    verify(target, record)
                else:
                    (stage / name).replace(target)
    print('Seven CMB wheel files verified', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', action='store_true', help='Base Python/MPI runtime for noise and observation tools')
    parser.add_argument('--cmb', action='store_true', help='Base runtime and the seven CMB dependency wheels')
    parser.add_argument('--templates', action='store_true', help='Seven PySM input templates from their provider')
    parser.add_argument('--theory-sources', action='store_true', help='Reference documents/code checked by the theory builder')
    parser.add_argument('--local-assets', type=Path, help='Use matching filenames in this directory instead of network downloads')
    parser.add_argument('--plan', action='store_true', help='List destinations and sizes without reading or writing downloads')
    parser.add_argument('--verify-only', action='store_true', help='Check installed download files without downloading')
    args = parser.parse_args()
    if not any([args.runtime, args.cmb, args.templates, args.theory_sources]):
        parser.error('select --runtime, --cmb, --templates and/or --theory-sources')
    if args.plan and args.verify_only:
        parser.error('--plan and --verify-only are mutually exclusive')
    assets = {r['filename']: r for r in json.loads((ROOT / 'release_assets.json').read_text())['files']}
    if args.runtime or args.cmb:
        record = assets['ali_sky_runtime_1.0.tar.gz']
        acquire(record, ROOT / 'runtime.tar.gz', args)
        if args.cmb:
            target = ROOT / 'cmb_sehgal2010/runtime_base.tar.gz'
            if args.plan:
                print('CMB base runtime will also be placed at', target.relative_to(ROOT))
            elif target.exists() or args.verify_only:
                verify(target, record)
            else:
                # Copy the verified archive; the two installed runtime prefixes are separate.
                shutil.copy2(ROOT / 'runtime.tar.gz', target)
                verify(target, record)
    if args.cmb:
        record = assets['cmb_wheels_1.0.tar.gz']
        archive = ROOT / '.downloads' / record['filename']
        acquire(record, archive, args)
        install_wheels(archive, args)
    if args.templates:
        manifest = json.loads((ROOT / 'runtime_evidence/template_manifest.json').read_text())
        for record in manifest['files']:
            row = {**record, 'filename': record['logical_name']}
            acquire(row, safe_child(ROOT / 'template_cache', record['logical_name']), args)
    if args.theory_sources:
        source_root = ROOT / 'cmb_sehgal2010/theory/source_evidence'
        for record in json.loads((source_root / 'source_manifest.json').read_text()):
            acquire(record, safe_child(source_root, record['filename']), args)


if __name__ == '__main__':
    main()
