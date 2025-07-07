import argparse, logging, yaml, os, re
from fairchem.core.units import mlip_unit
from fairchem.core import FAIRChemCalculator
from ase import Atoms
from ase.optimize import LBFGS
from ase.vibrations import Vibrations
from ase.io import read, write

MODEL = {
    'small': 'model/UMA/uma-s-1p1',
    'medium': 'model/UMA/uma-m-1p1',
    'large': ''
    }

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

def parse_args_into_config(param: dict):
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
            set_nested_value(param, key, value)

    return param

def setup_logger(logger_name: str, log_file: str, level=logging.INFO) -> logging.Logger:
    logger = logging.getLogger(logger_name)
    logger.setLevel(level)

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

class QCCalculator:
    def __init__(
        self,
        atoms: Atoms,
        calc: FAIRChemCalculator,
        logger: logging.Logger,
        ):
        self.atoms = atoms
        self.atoms.calc = calc
        self.logger = logger
        self.log_file = logger.handlers[0].baseFilename

    def optimization(
        self,
        traj_file: str | None,
        optimized_atoms_file: str | None,
        fmax: float=0.05,
        steps: int=100_000_000
        ):
        self.logger.info('Start optimization...')
        optimizer = LBFGS(self.atoms, trajectory=traj_file, logfile=self.log_file)
        optimizer.run(fmax=fmax, steps=steps)
        optimized_energy = self.atoms.get_potential_energy()
        self.logger.info(f'Optimized energy: {optimized_energy} eV')
        self.logger.info(f'Standard orientation: \n{gen_xyz(self.atoms)}')

        if optimized_atoms_file:
            write(optimized_atoms_file, self.atoms)
            self.logger.info(f'Optimized structure has been saved to {optimized_atoms_file}')

        self.logger.info('End optimization...')

        return optimized_energy

    def spe_calculation(
        self,
        ):
        self.logger.info('Start SPE calculation...')
        energy = self.atoms.get_potential_energy()
        self.logger.info(f'Energy: {energy} eV')
        self.logger.info('End SPE calculation...')

        return energy

    def vib_calculation(
        self,
        name: str
        ):
        self.logger.info('Start vib calculation...')
        vib = Vibrations(atoms=self.atoms, name=name)
        vib.run()
        vib.summary(log=self.log_file) # 打印频率总结
        vib.write_jmol() # 写入振动模式文件，可在Jmol中查看
        vib.clean()
        self.logger.info('End vib calculation...')


def main(
    param: dict[str, dict]
    ):
    if param['general']['log_file']:
        log_file = param['general']['log_file']
    else:
        name, suffix = os.path.splitext(param['structure']['file_path'])
        log_file = name + '.log'
    logger = setup_logger('QC', log_file)

    try:
        logger.info('Loading calculator...')
        if not MODEL[param['model']['size']]:
            raise ValueError('Model path is empty or invalid.')
        predictor = mlip_unit.load_predict_unit(
            path=MODEL[param['model']['size']],
            device=param['model']['device'])
        calc = FAIRChemCalculator(predictor, task_name='omol')
        logger.info('Done!')

        atoms = gen_atoms(
            param['structure']['file_path'],
            param['structure']['charge'],
            param['structure']['spin'])

        QCC = QCCalculator(atoms=atoms, calc=calc, logger=logger)
        logger.info(f'file path: {param["structure"]["file_path"]}')
        logger.info(f'Input orientation: \n{gen_xyz(atoms)}')
        logger.info(f'charge: {atoms.info["charge"]}')
        logger.info(f'spin: {atoms.info["spin"]}')
        logger.info(f'model: {MODEL[param["model"]["size"]]}')
        for task, args in param['task'].items():
            args: dict
            if args['run']:
                kwargs = {k: v for k, v in args.items() if k != 'run'}
                getattr(QCC, task)(**kwargs)

    except Exception as e:
        logger.error("Error occurred during execution:", exc_info=True)
        raise
    finally:
        if param['model'].get('device') == 'cuda':
            import torch
            torch.cuda.empty_cache()

if __name__ == '__main__':
    # TODO: 支持批量计算，类似 `qg09 -a`
    with open('config.yml', 'r', encoding='utf-8') as f:
        param: dict = yaml.full_load(f)
    param = parse_args_into_config(param)
    main(param)

