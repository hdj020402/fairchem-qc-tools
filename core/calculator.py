import logging
from fairchem.core import FAIRChemCalculator
from ase import Atoms
from ase.optimize import LBFGS
from ase.vibrations import Vibrations
from ase.io import write
from core.utils import gen_xyz

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
