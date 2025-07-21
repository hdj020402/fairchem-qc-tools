#!/usr/bin/env python

import sys
import os
import argparse
import json
import subprocess

def str2bool(v: str):
    if v is None:
        return None
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    if v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    raise argparse.ArgumentTypeError('Boolean value expected.')

def build_run_parser(subparsers: argparse._SubParsersAction):
    parser_run: argparse.ArgumentParser = subparsers.add_parser(
        'run',
        help='Submit structure and task config to /run-web',
        description='Submit structure and task config to /run-web',
        formatter_class=argparse.RawTextHelpFormatter
        )
    parser_run.epilog = """
Example:
python web.py run \\
    --structure.file_path=test/test.xyz \\
    --task.optimization.run \\
    --task.optimization.fmax=0.02 \\
    --task.spe_calculation.run
"""

    service_group = parser_run.add_argument_group('Service Options')
    structure_group = parser_run.add_argument_group('Structure Options')
    task_group = parser_run.add_argument_group('Task Options')

    # Service section
    service_group.add_argument('--service.host', type=str, help='Service host address')
    service_group.add_argument('--service.port', type=int, help='Service port number')

    # Structure section
    structure_group.add_argument('--structure.file_path', type=str, required=True,
                            help='Path to input structure file (.xyz or .gjf)')
    structure_group.add_argument('--structure.charge', type=int, help='Charge of the molecule')
    structure_group.add_argument('--structure.spin', type=int, help='Spin multiplicity of the molecule')

    # Task - Optimization
    task_group.add_argument('--task.optimization.run', nargs='?', const=True, default=None,
                            type=str2bool, help='Whether to run optimization task')
    task_group.add_argument('--task.optimization.traj', nargs='?', const=True, default=None,
                            type=str2bool, help='Whether to save the optimization trajectory')
    task_group.add_argument('--task.optimization.optimized_structure', nargs='?', const=True, default=None,
                            type=str2bool, help='Whether to save the optimized structure')
    task_group.add_argument('--task.optimization.fmax', type=float,
                            help='Force convergence threshold (eV/Å)')
    task_group.add_argument('--task.optimization.steps', type=int,
                            help='Maximum number of optimization steps')

    # Task - SPE Calculation
    task_group.add_argument('--task.spe_calculation.run', nargs='?', const=True, default=None,
                            type=str2bool, help='Whether to run single-point energy calculation')

    # Task - Vibration Calculation
    task_group.add_argument('--task.vib_calculation.run', nargs='?', const=True, default=None,
                            type=str2bool, help='Whether to run vibrational analysis')

    # Task - Force Calculation
    task_group.add_argument('--task.force_calculation.run', nargs='?', const=True, default=None,
                            type=str2bool, help='Whether to run force calculation')

    return parser_run

def build_log_parser(subparsers: argparse._SubParsersAction):
    parser_log: argparse.ArgumentParser = subparsers.add_parser(
        'log',
        help='Get log from /logs/{task_id}',
        description='Get log from /logs/{task_id}',
        formatter_class=argparse.RawTextHelpFormatter
        )
    parser_log.epilog = """
Example:
python web.py log \\
    --task_id abc123 \\
    --service.host 0.0.0.0 \\
    --service.port 8000
"""

    log_group = parser_log.add_argument_group('Log Options')
    service_group = parser_log.add_argument_group('Service Options')

    log_group.add_argument('--task_id', type=str, required=True, help='Task ID to fetch log')
    service_group.add_argument('--service.host', type=str, default='localhost',
                            help='Service host address')
    service_group.add_argument('--service.port', type=int, default=8000,
                            help='Service port number')
    return parser_log

def build_download_parser(subparsers: argparse._SubParsersAction):
    parser_download: argparse.ArgumentParser = subparsers.add_parser(
        'download',
        help='Download task folder from /download/{task_id}',
        description='Download task folder from /download/{task_id}',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser_download.epilog = """
Example:
  python web.py download \\
    --task_id abc123 \\
    --archive_format zip \\
    --output abc123.zip
"""

    download_group = parser_download.add_argument_group('Download Options')
    service_group = parser_download.add_argument_group('Service Options')

    download_group.add_argument('--task_id', type=str, required=True,
                                help='Task ID to download')
    download_group.add_argument('--archive_format', type=str, default='zip',
                                choices=['zip', 'tar.gz'],
                                help='Archive format (default: zip)')
    download_group.add_argument('--output', type=str, default=None,
                                help='Output file path (default: {task_id}.{format})')
    service_group.add_argument('--service.host', type=str, default='localhost',
                               help='Service host address')
    service_group.add_argument('--service.port', type=int, default=8000,
                               help='Service port number')

    return parser_download

def set_nested_value(d: dict, key: str, value):
    keys = key.split('.')
    for k in keys[:-1]:
        d = d.setdefault(k, {})
    d[keys[-1]] = value

def parse_args_into_config():
    parser = argparse.ArgumentParser(description="QC Task CLI Tool")
    subparsers = parser.add_subparsers(dest='command', required=True)

    run_parser = build_run_parser(subparsers)
    log_parser = build_log_parser(subparsers)
    download_parser = build_download_parser(subparsers)

    args = parser.parse_args()

    if args.command == 'run':
        config = {}
        for key, value in vars(args).items():
            if key not in ['command'] and value is not None:
                set_nested_value(config, key, value)
        return 'run', config
    elif args.command == 'log':
        return 'log', args
    elif args.command == 'download':
        return 'download', args
    else:
        raise ValueError("Unknown command")

def main():
    command, data = parse_args_into_config()

    if command == 'run':
        config = data

        structure_file_path = config.get('structure', {}).get('file_path')
        if not structure_file_path or not os.path.isfile(structure_file_path):
            print("Error: structure.file_path is required and must be a valid file.")
            sys.exit(1)

        host = config.get('service', {}).get('host', 'localhost')
        port = config.get('service', {}).get('port', 8000)
        url = f"http://{host}:{port}/run-web"

        task_config_json = json.dumps(config, ensure_ascii=False)

        curl_cmd = [
            'curl', '-X', 'POST',
            url,
            '-H', 'Content-Type: multipart/form-data',
            '-F', f'structure_file=@{structure_file_path}',
            '-F', f'task_config={task_config_json}'
        ]

    elif command == 'log':
        host = getattr(data, 'service.host', 'localhost')
        port = getattr(data, 'service.port', 8000)
        task_id = data.task_id
        url = f"http://{host}:{port}/logs/{task_id}"

        curl_cmd = [
            'curl', '-X', 'GET',
            url
        ]

    elif command == 'download':
        host = getattr(data, 'service.host', 'localhost')
        port = getattr(data, 'service.port', 8000)
        task_id = data.task_id
        archive_format = data.archive_format
        output = data.output or f"{task_id}.{archive_format}"

        url = f"http://{host}:{port}/download/{task_id}?archive_format={archive_format}"

        curl_cmd = [
            'curl', '-X', 'GET',
            url,
            '--output', output
        ]

    print("Executing curl command:")
    print(' '.join(curl_cmd))
    print()

    try:
        result = subprocess.check_output(curl_cmd, stderr=subprocess.STDOUT)
        print("Response Body:")
        print(result.decode('utf-8'))
    except subprocess.CalledProcessError as e:
        print("Error executing curl:")
        print(e.output.decode('utf-8'))
        sys.exit(1)

if __name__ == '__main__':
    main()
