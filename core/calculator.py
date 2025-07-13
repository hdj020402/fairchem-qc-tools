import logging, os
from fairchem.core import FAIRChemCalculator
from ase import Atoms
from ase.optimize import LBFGS
from ase.vibrations import Vibrations
from ase.io import write
from core.utils import gen_xyz, gen_forces, calc_forces_stats

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
        self.log_dir = os.path.dirname(self.log_file)

    def optimization(
        self,
        traj: bool,
        optimized_structure: bool,
        fmax: float=0.05,
        steps: int=100_000_000
        ):
        self.logger.info('Start optimization...')
        if traj:
            traj_file = os.path.join(self.log_dir, 'optimization.traj')
        else:
            traj_file = None
        optimizer = LBFGS(self.atoms, trajectory=traj_file, logfile=self.log_file)
        optimizer.run(fmax=fmax, steps=steps)
        optimized_energy = self.atoms.get_potential_energy()
        self.logger.info(f'Optimized energy: {optimized_energy} eV')
        self.logger.info(f'Standard orientation: \n{gen_xyz(self.atoms)}')

        if optimized_structure:
            optimized_atoms_file = os.path.join(self.log_dir, 'optimized_structure.xyz')
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

    def vib_calculation(self):
        self.logger.info('Start vib calculation...')
        name = os.path.join(self.log_dir, 'vib')
        vib = Vibrations(atoms=self.atoms, name=name)
        vib.run()
        vib.summary(log=self.log_file) # 打印频率总结
        vib.write_jmol() # 写入振动模式文件，可在Jmol中查看
        # vib.clean()
        self.logger.info('End vib calculation...')

    def force_calculation(self):
        self.logger.info('Start force calculation...')
        force = self.atoms.get_forces()
        self.logger.info(f'Forces: \n{gen_forces(self.atoms)}')
        max_force_per_atom, rms_force_per_atom, max_force_component, rms_force_component = calc_forces_stats(force)
        self.logger.info(f'Max force per atom: {max_force_per_atom:.6f}')
        self.logger.info(f'RMS force per atom: {rms_force_per_atom:.6f}')
        self.logger.info(f'Max force component: {max_force_component:.6f}')
        self.logger.info(f'RMS force component: {rms_force_component:.6f}')
        self.logger.info('End force calculation...')
        return force, max_force_per_atom, rms_force_per_atom, max_force_component, rms_force_component
