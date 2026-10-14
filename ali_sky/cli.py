"""One-command noise campaigns and forward observation synthesis."""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

from .config import load_config, frequency_tag
from .io import atomic_json


def under_mpi():
    return any(name in os.environ for name in ('OMPI_COMM_WORLD_SIZE', 'PMI_SIZE', 'PMIX_RANK', 'MPI_LOCALNRANKS'))


def noise_mpi(config):
    """Strided (band,ID) task allocation; seeds never use rank or task order."""
    from mpi4py import MPI
    from .noise import generate_noisefile
    comm = MPI.COMM_WORLD
    rank, size = comm.Get_rank(), comm.Get_size()
    tasks = [(band, rid) for band in range(len(config['channels'])) for rid in config['ids']]
    results, errors = [], []
    for task in tasks[rank::size]:
        try:
            result = generate_noisefile(config, task)
            results.append({'frequency_GHz': config['channels'][task[0]]['frequency_GHz'],
                            'id': task[1], 'output': result['output'], 'reused': result.get('reused', False)})
            print(f'rank {rank}: {result["output"]["path"]}', flush=True)
        except Exception as exc:
            errors.append({'rank': rank, 'task': list(task), 'error': f'{type(exc).__name__}: {exc}'})
    all_results = comm.gather(results, root=0)
    all_errors = comm.gather(errors, root=0)
    failed = None
    if rank == 0:
        flat_results = [x for group in all_results for x in group]
        flat_errors = [x for group in all_errors for x in group]
        atomic_json(Path(config['output_dir']) / 'noise_campaign.json',
                    {'status': 'FAILED' if flat_errors else 'COMPLETE', 'configuration': config,
                     'MPI_ranks': size, 'tasks': sorted(flat_results, key=lambda x: (x['frequency_GHz'], x['id'])),
                     'errors': flat_errors})
        failed = bool(flat_errors)
        if failed:
            print(json.dumps(flat_errors, indent=2), file=sys.stderr)
    failed = comm.bcast(failed, root=0)
    return 1 if failed else 0


def launch_noise(config):
    if under_mpi():
        return noise_mpi(config)
    ranks = config['mpi_ranks']
    mpiexec = os.environ.get('ALI_SKY_MPIEXEC', str(Path(sys.executable).parent / 'mpiexec'))
    if not Path(mpiexec).exists():
        raise FileNotFoundError(f'Matching MPI launcher not found: {mpiexec}; use bundled runtime/setup_env.sh')
    out = Path(config['output_dir'])
    out.mkdir(parents=True, exist_ok=True)
    task_file = out / 'noise_configuration.json'
    atomic_json(task_file, config)
    command = [mpiexec, '-n', str(ranks), sys.executable, '-m', 'ali_sky.cli', '_noise-worker', '--config', str(task_file)]
    print(f'Launching {ranks} MPI ranks for {len(config["channels"])} frequency bands and {len(config["ids"])} IDs', flush=True)
    return subprocess.run(command, check=False).returncode


def run_observe(config):
    if under_mpi():
        raise ValueError('run_observe.sh is a serial coordinator; automatic multiband noise launches its own MPI workers')
    from .observation import observe_one, preflight
    # Resolve sky/mask inputs before any noise generation; missing explicit archive IDs are errors.
    preflight(config)
    if config['add_noise'] and config['custom_noise'] and not config['noise_dir']:
        print('NOTICE: custom noise parameters were supplied without a noise directory. '
              f'Multifrequency noise will be generated in parallel using {config["mpi_ranks"]} MPI ranks '
              '(default 5); the sky beam is not applied to noise.', flush=True)
        generated = Path(config['output_dir']) / 'generated_noise'
        noise_config = {k: config[k] for k in ('nside', 'seed', 'instrument', 'scale_factor', 'coordinate', 'overwrite', 'mpi_ranks')}
        noise_config.update({'channels': [{k: v for k, v in c.items() if k != 'beam_transfer_path'} for c in config['channels']],
                             'ids': list(dict.fromkeys(config['noise_ids'])), 'output_dir': str(generated), 'mode': 'noise'})
        code = launch_noise(noise_config)
        if code:
            return code
        config = {**config, 'noise_dir': str(generated), 'noise_coordinate': config['coordinate'],
                  'noise_pattern': '{instrument}/{frequency}/noise_{id:05d}.fits'}
    results = []
    for channel in config['channels']:
        for rid, noise_id in zip(config['ids'], config['noise_ids']):
            result = observe_one(config, channel, rid, noise_id)
            results.append({'frequency_GHz': channel['frequency_GHz'], 'sky_id': rid,
                            'noise_id': noise_id, 'output': result['output']})
    atomic_json(Path(config['output_dir']) / 'observation_campaign.json',
                {'status': 'COMPLETE', 'configuration': config, 'products': results})
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ('noise', 'observe', 'validate-config', '_noise-worker'):
        p = sub.add_parser(name, help=argparse.SUPPRESS if name.startswith('_') else None)
        p.add_argument('--config', required=True, help='JSON configuration')
        if name == 'validate-config':
            p.add_argument('--mode', choices=('noise', 'observe'), default='observe')
        else:
            p.add_argument('--n-realizations', type=int)
            p.add_argument('--start-id', type=int)
            p.add_argument('--seed', type=int)
            p.add_argument('--output-dir')
            p.add_argument('--mpi-ranks', type=int)
    args = parser.parse_args(argv)
    try:
        if args.command == '_noise-worker':
            return noise_mpi(json.loads(Path(args.config).read_text()))
        mode = args.mode if args.command == 'validate-config' else args.command
        overrides = {name: getattr(args, name, None) for name in ('n_realizations', 'start_id', 'seed', 'output_dir', 'mpi_ranks')}
        config = load_config(args.config, mode, overrides)
        if args.command == 'validate-config':
            print(json.dumps(config, indent=2))
            return 0
        return launch_noise(config) if mode == 'noise' else run_observe(config)
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
