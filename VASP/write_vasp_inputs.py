import ase
from pathlib import Path

def write_vasp_inputs(atoms: ase.Atoms,
                      save_path: Path,
                      NSW: int = 200
    ) -> None:
    r'''
    Write vasp inputs which have same settings with acpy.

    Args:
        atoms: ase.Atoms
        save_path, Path

    Returns:
        None
    '''
    import os
    from ase.calculators.vasp import Vasp
    from ..utils.globals import MAGMOMS

    os.environ['VASP_PP_PATH'] = '/home/shared/programs/vasp/vasp_pp'
    
    magmom_temp = [MAGMOMS.get(atom.symbol, 0.0) for atom in atoms]

    vasp_calculator = Vasp(encut=400, 
                           xc = 'rpbe', 
                           algo='fast', 
                           prec='Normal', 
                           ismear=1, 
                           sigma=0.1,
                           ediff=1e-5, 
                           ispin=2,
                           kpts=[4,4,1],
                           ibrion=2,
                           potim=0.1,
                           ediffg=-0.05,
                           isif=0, 
                           nsw=200,
                           nwrite=1,
                           lcharg=False,
                           lwave=False,
                           lvtot=False,
                           istart=0, 
                           magmom=magmom_temp, 
                           npar=2,
                           lreal='auto',
                           setups={"base":"recommended"})
    
    # write vasp input file
    vasp_calculator.directory = save_path
    vasp_calculator.write_input(atoms)

