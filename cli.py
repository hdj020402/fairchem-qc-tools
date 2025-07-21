import yaml, argparse, os
import requests
from core.utils import str2bool

def parse_args_into_config(config: dict):
    parser = argparse.ArgumentParser(description="Override config.yml parameters via command line.")

    # Service section
    parser.add_argument('--service.host', type=str, help='')
    parser.add_argument('--service.port', type=int, help='')

    # Structure section
    parser.add_argument('--structure.file_path', type=str,
                        help='Path to input structure file (.xyz or .gjf)')
    parser.add_argument('--structure.charge', type=int,
                        help='Charge of the molecule')
    parser.add_argument('--structure.spin', type=int,
                        help='Spin multiplicity of the molecule')

    # Task - Optimization
    parser.add_argument('--task.optimization.run', nargs='?', const=True, default=None,
                        type=str2bool, help='Whether to run optimization task')
    parser.add_argument('--task.optimization.traj', nargs='?', const=True, default=None,
                        type=str2bool, help='Whether to save the optimization trajectory')
    parser.add_argument('--task.optimization.optimized_structure', nargs='?', const=True, default=None,
                        type=str2bool, help='Whether to save the optimized structure')
    parser.add_argument('--task.optimization.fmax', type=float,
                        help='Force convergence threshold (eV/Å)')
    parser.add_argument('--task.optimization.steps', type=int,
                        help='Maximum number of optimization steps')

    # Task - SPE Calculation
    parser.add_argument('--task.spe_calculation.run', nargs='?', const=True, default=None,
                        type=str2bool, help='Whether to run single-point energy calculation')

    # Task - Vibration Calculation
    parser.add_argument('--task.vib_calculation.run', nargs='?', const=True, default=None,
                        type=str2bool, help='Whether to run vibrational analysis')

    # Task - Force Calculation
    parser.add_argument('--task.force_calculation.run', nargs='?', const=True, default=None,
                            type=str2bool, help='Whether to run force calculation')

    args = parser.parse_args()

    # Apply overrides using dot notation
    def set_nested_value(d: dict, key: str, value):
        keys = key.split('.')
        for k in keys[:-1]:
            d = d.setdefault(k, {})
        d[keys[-1]] = value

    for key, value in vars(args).items():
        if value is not None:
            set_nested_value(config, key, value)

    return config

def main():
    TASK_CONFIG_PATH = "task_config.yml"
    with open(TASK_CONFIG_PATH, "r", encoding="utf-8") as f:
        task_config: dict[str, dict] = yaml.safe_load(f)

    task_config = parse_args_into_config(task_config)
    structure_file_path = task_config.get("structure", {}).get("file_path")
    if structure_file_path and os.path.isfile(structure_file_path):
        task_config['structure']['file_path'] = os.path.abspath(structure_file_path)

    host = task_config.get("service", {}).get("host", "localhost")
    port = task_config.get("service", {}).get("port", 8000)
    url = f"http://{host}:{port}/run-cli"

    response = requests.post(url, json=task_config)

    print("Status Code:", response.status_code)
    print("Response Body:", response.json())

if __name__ == "__main__":
    main()
