import logging, argparse, re, os, time
import string, secrets
from ase import Atoms
from ase.io import read

def setup_logger(logger_name: str, log_file: str, level=logging.INFO) -> logging.Logger:
    logger = logging.getLogger(logger_name)
    logger.setLevel(level)

    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    file_handler = logging.FileHandler(log_file, mode='a')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.propagate = False

    return logger

def gen_atoms(file_path: str, charge: int=0, spin: int=1):
    atoms = read(file_path)
    if file_path.endswith('.xyz'):
        ...
    elif file_path.endswith('.gjf'):
        with open(file_path) as f:
            content = f.read()
        match = re.search(r'#.*\n\n*.*\n\n*(\d+)\s+(\d+)', content)
        if match:
            charge = int(match.group(1))
            spin = int(match.group(2))
        else:
            raise ValueError("Cannot parse charge and spin from gjf file")
    else:
        raise ValueError("Unsupported file format. Only .xyz and .gjf files are allowed.")

    atoms.info['charge'] = charge
    atoms.info['spin'] = spin

    return atoms

def gen_xyz(atoms: Atoms):
    atoms_list = atoms.get_chemical_symbols()
    positions = atoms.positions
    xyz = [f'{len(atoms_list)}', '']
    for atom, pos in zip(atoms_list, positions):
        xyz.append(f'{atom} {pos[0]: .8f} {pos[1]: .8f} {pos[2]: .8f}')
    return '\n'.join(xyz)

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

def parse_args_into_config(config: dict):
    parser = argparse.ArgumentParser(description="Override config.yml parameters via command line.")

    # General section
    parser.add_argument('--general.log_file', type=str,
                        help='Path to the log file. Default: <file_path>.log')

    # Model section
    parser.add_argument('--model.size', type=str, choices=['small', 'medium', 'large'],
                        help='Size of the model to use. Choices: small, medium, large')
    parser.add_argument('--model.device', type=str, choices=['cpu', 'cuda'],
                        help='Device to run the model on. Choices: cpu, cuda')

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
    parser.add_argument('--task.optimization.traj_file', type=str,
                        help='Path to save optimization trajectory (.traj)')
    parser.add_argument('--task.optimization.optimized_atoms_file', type=str,
                        help='Path to save optimized atoms')
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
    parser.add_argument('--task.vib_calculation.name', type=str,
                        help='Base name for vibration output files')

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

def generate_uid(length: int=6):
    characters = string.ascii_letters + string.digits
    return ''.join(secrets.choice(characters) for _ in range(length))

def build_task_log_info(filepath: str):
    filename = os.path.splitext(os.path.basename(filepath))[0]
    safe_filename = re.sub(r'[^\w\-]', '_', filename)
    timestamp = time.strftime('%Y%m%d_%H%M%S')
    uid = generate_uid()
    task_name = f'{safe_filename}_{timestamp}_{uid}'

    log_dir = os.path.join(os.path.dirname(filepath), "logs")
    os.makedirs(log_dir, exist_ok=True)

    log_file = os.path.join(log_dir, f"{task_name}.log")

    return log_file, task_name
