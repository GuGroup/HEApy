import ase
from pathlib import Path

def write_vasp_inputs(atoms: ase.Atoms,
                      save_path: Path,
                      NSW: int = 200
                      calculation_type: str = "slab",
    ) -> None:
    r'''
    Write vasp inputs in pre-determined settings.

    Args:
        atoms: ase.Atoms
        save_path: Path
        calculation_type: str 

    Returns:
        None
    '''
    import os
    from ase.calculators.vasp import Vasp
    from ..utils.globals import MAGMOMS

    assert calculation_type in ["slab", "bulk"], f"Choose calculation type between 'slab' or 'bulk'."

    os.environ['VASP_PP_PATH'] = '/home/shared/programs/vasp/vasp_pp'
    
    magmom_temp = [MAGMOMS.get(atom.symbol, 0.0) for atom in atoms]
    
    if calculation_type == bulk:
        vasp_calculator = Vasp(encut=400, 
                               xc = 'rpbe', 
                               algo='fast', 
                               prec='Normal', 
                               ismear=1, 
                               sigma=0.1,
                               ediff=1e-5, 
                               ispin=2,
                               kpts=[4,4,4],
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

    elif calculation_type == "slab":
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

