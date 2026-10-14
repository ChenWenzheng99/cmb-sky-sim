#!/usr/bin/env python3
"""Verify 43 vendored sources, 41 installed Python/data files and the ducc route.

This read-only check imports and inspects the installed library; it does not
generate maps, run a harmonic transform, compile Fortran, or create a ready
marker. The two .f90 files are preserved source-only archive members.
"""
import argparse
import hashlib
import importlib.metadata
import inspect
import json
import os
from pathlib import Path
import platform
import sys
import tempfile


SOURCE_ONLY_F90 = {'lenspyx/fortran/remapping.f90', 'lenspyx/shts/shts.f90'}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def inside(path, parent):
    try:
        Path(path).resolve().relative_to(Path(parent).resolve())
        return True
    except ValueError:
        return False


def function_record(function, package_root):
    filename = Path(inspect.getsourcefile(function)).resolve()
    source_lines, first_line = inspect.getsourcelines(function)
    return {'relative_path': filename.relative_to(package_root).as_posix(),
            'first_line': first_line, 'last_line': first_line + len(source_lines) - 1,
            'function': function.__qualname__, 'source_sha256': sha256(filename)}


def verify(args):
    root = args.root.resolve()
    vendor = (args.vendor_root or root / 'vendor/lenspyx-2.0.0').resolve()
    manifest_path = (args.manifest or root / 'evidence/lenspyx_vendor_manifest.json').resolve()
    runtime = (args.runtime_root or root / 'runtime').resolve()
    manifest = json.loads(manifest_path.read_text())
    records = manifest['files']
    paths = [item['path'] for item in records]
    require(len(records) == 43 and len(set(paths)) == 43, 'Expected 43 distinct retained source paths')
    require({path for path in paths if path.endswith('.f90')} == SOURCE_ONLY_F90,
            'Unexpected source-only Fortran inventory')

    import lenspyx
    import lenspyx.lensing as lensing
    from lenspyx.remapping import deflection_028, deflection_029, utils_geom
    import ducc0
    import numpy

    installed = Path(lenspyx.__file__).resolve().parent
    require(inside(installed, runtime), 'lenspyx import is outside the dedicated runtime')
    require(inside(ducc0.__file__, runtime) and inside(numpy.__file__, runtime),
            'A numerical dependency was imported outside the dedicated runtime')
    versions = {name: importlib.metadata.version(name) for name in ['lenspyx', 'ducc0', 'numpy']}
    require(versions == {'lenspyx':'2.0.0', 'ducc0':'0.38.0', 'numpy':'1.26.4'},
            f'Unexpected execution-library versions: {versions}')

    file_checks = []
    for record in records:
        relative = Path(record['path'])
        require(not relative.is_absolute() and '..' not in relative.parts and relative.parts[0] == 'lenspyx',
                f'Invalid manifest relative path: {relative}')
        vendored_file = vendor / relative
        require(vendored_file.is_file(), f'Missing vendored source: {relative}')
        vendor_sha = sha256(vendored_file)
        require(vendored_file.stat().st_size == record['bytes'] and vendor_sha == record['sha256'],
                f'Vendored source differs from original manifest: {relative}')
        installed_file = installed / Path(*relative.parts[1:])
        source_only = relative.as_posix() in SOURCE_ONLY_F90
        installed_sha = None
        if source_only:
            # Source retention is mandatory; installation/compilation is not.
            if installed_file.is_file():
                installed_sha = sha256(installed_file)
                require(installed_sha == record['sha256'], f'Installed optional F90 source differs: {relative}')
        else:
            require(relative.suffix in {'.py', '.dat', '.ini'}, f'Unexpected installed-file type: {relative}')
            require(installed_file.is_file(), f'Missing installed Python/data file: {relative}')
            installed_sha = sha256(installed_file)
            require(installed_file.stat().st_size == record['bytes'] and installed_sha == record['sha256'],
                    f'Installed Python/data file differs from original manifest: {relative}')
        file_checks.append({'relative_path': relative.as_posix(), 'bytes': record['bytes'],
                            'expected_sha256': record['sha256'], 'vendored_sha256': vendor_sha,
                            'classification': 'source_only_Fortran' if source_only else 'installed_Python_or_data',
                            'installed_present': installed_file.is_file(), 'installed_sha256': installed_sha,
                            'vendored_exact': True, 'installed_exact_when_required': True})

    require(lenspyx.alm2lenmap is lensing.alm2lenmap, 'Unexpected top-level remapping entry point')
    require(lensing.deflection is deflection_029.deflection, 'alm2lenmap does not use deflection_029')
    require(deflection_029.deflection.__base__ is deflection_028.deflection,
            'Unexpected inherited pointing implementation')
    require(deflection_029.synthesis_general is ducc0.sht.experimental.synthesis_general,
            'General synthesis does not resolve to ducc')
    require(utils_geom.synthesis is ducc0.sht.experimental.synthesis,
            'Geometry synthesis does not resolve to ducc')
    require(deflection_028.get_deflected_angles is ducc0.misc.get_deflected_angles,
            'Pointing calculation does not resolve to ducc')
    flags = {'deflection028_HAS_DUCCPOINTING': bool(deflection_028.HAS_DUCCPOINTING),
             'deflection028_HAS_DUCCROTATE': bool(deflection_028.HAS_DUCCROTATE),
             'deflection028_HAS_DUCCGRADONLY': bool(deflection_028.HAS_DUCCGRADONLY),
             'deflection029_HAS_DUCCROTATE': bool(deflection_029.HAS_DUCCROTATE),
             'deflection029_HAS_DUCCGRADONLY': bool(deflection_029.HAS_DUCCGRADONLY),
             'deflection028_HAS_FORTRAN': bool(deflection_028.HAS_FORTRAN),
             'deflection029_HAS_FORTRAN': bool(deflection_029.HAS_FORTRAN)}
    require(all(value for name, value in flags.items() if 'DUCC' in name), 'Missing ducc route capability')
    require(not flags['deflection028_HAS_FORTRAN'] and not flags['deflection029_HAS_FORTRAN'],
            'Unexpected optional compiled Fortran backend in the pure-Python wheel')
    compiled_files = sorted(path.relative_to(installed).as_posix() for path in installed.rglob('*')
                            if path.is_file() and (path.suffix in {'.so', '.pyd', '.dylib'}))
    require(not compiled_files, f'Unexpected compiled extension under lenspyx: {compiled_files}')
    require('lenspyx.fortran.remapping' not in sys.modules and 'lenspyx.shts.shts' not in sys.modules,
            'Legacy Fortran execution route was loaded')
    route_functions = [lensing.alm2lenmap, deflection_029.deflection.gclm2lenmap,
                       deflection_028.deflection._build_d1, deflection_028.deflection._build_angles,
                       utils_geom.Geom.synthesis]
    return {
        'schema': 'lenspyx_vendored_and_execution_sources_v1',
        'status': 'PASS_VENDORED43_INSTALLED41_DUCC_EXECUTION_ROUTE',
        'release_root': str(root), 'runtime_prefix': str(runtime), 'installed_lenspyx': str(installed),
        'vendored_root': str(vendor), 'manifest_sha256': sha256(manifest_path),
        'helper_sha256': sha256(__file__), 'versions': versions, 'python': platform.python_version(),
        'vendored_source_files_checked': 43, 'installed_python_data_files_checked': 41,
        'source_only_Fortran_files_checked': 2, 'source_only_Fortran_paths': sorted(SOURCE_ONLY_F90),
        'file_checks': file_checks, 'backend_flags': flags,
        'compiled_lenspyx_extensions_found': compiled_files,
        'active_source_route': [function_record(function, installed) for function in route_functions],
        'route': ['lenspyx.alm2lenmap', 'lensing.deflection_029.deflection',
                  'ducc0.sht.experimental.synthesis_general', 'ducc0.misc.get_deflected_angles',
                  'ducc0.misc.lensing_rotate'],
        'conclusion': 'The retained alm2lenmap forward T/QU path uses ducc synthesis, pointing and polarization rotation. The two original F90 files are preserved in vendor, are not compiled into this wheel and are not used by this execution route.',
        'scope': 'Read-only original-file hashes, dedicated-prefix imports and active-function/capability inspection; no native map generation, SHT, compilation, environment mutation or ready-marker write.',
        'limitations': ['The check does not claim legacy APIs that directly require optional Fortran extensions are available.',
                        'A small-array or native remapping numerical validation is separate from this source/backend check.']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--vendor-root', type=Path)
    parser.add_argument('--manifest', type=Path)
    parser.add_argument('--runtime-root', type=Path)
    parser.add_argument('--report', type=Path, help='Optional explicit evidence JSON; never writes a ready marker')
    args = parser.parse_args()
    try:
        result = verify(args)
        exit_code = 0
    except Exception as error:
        result = {'schema':'lenspyx_vendored_and_execution_sources_v1', 'status':'FAIL',
                  'error_type':type(error).__name__, 'error':str(error), 'helper_sha256':sha256(__file__)}
        exit_code = 1
    encoded = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile('w', dir=args.report.parent, prefix='.lenspyx_check_', delete=False) as handle:
            handle.write(encoded)
            temporary = handle.name
        os.replace(temporary, args.report)
    print(encoded, end='')
    return exit_code


if __name__ == '__main__':
    raise SystemExit(main())
