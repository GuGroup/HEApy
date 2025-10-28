# ytk/ktg
# ytk/ktg

import ase
import numpy as np

from pathlib import Path
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.build import fcc111, add_adsorbate

from ..utils.logging import setup_logger

def generate_random_slab(constituent_atoms: list[str],
                         calculation_path: Path,
                         size: tuple | list,
                         constraint: bool | None,
                         adsorbate: Path | None,
                         adsorbate_positions: tuple | None,
                         count: int = 1,
    ) -> None:
    r'''
    ===========================================================================
    Generate random high entropy alloy slab (+ adsorbate)
    Half bottom layers are fixed if constraint

    Args:
        constituent_atoms: Atom symbols in high-entropy alloy, ["Co", "Ni", "Cu", "Ru"]
        calculation_path: where the calculations would be saved
        size: size of slab, (2, 2, 4) or (3, 3, 6)
        constraint: whether fix atoms on low half
        adsorbate: the path of adsorbate
        adsorption_positions: adsorption (x, y) positions, height is set on 2
        count: the number of total calculation

    Returns:
        None

    ===========================================================================
    '''
    from .dataset_utils import _set_constraint, get_binding_atom
    from ..utils.globals import LATTICE_CONSTANTS
    from ..VASP.write_vasp_inputs import write_vasp_inputs

    np.random.seed(1308)
    logger = setup_logger("HEA")
    
    # Set lattice constant
    lattice = 0
    for constituent_atom in constituent_atoms:
        lattice += LATTICE_CONSTANTS[constituent_atom]
    lattice = lattice / len(constituent_atoms)
    
    # Set calculation dict for how many calculations have done
    calculation_count = {}
    for constituent_atom in constituent_atoms:
        calculation_count[constituent_atom] = 0

    
    for _ in range(count):
        # Make default slab
        slab = fcc111("Cu", 
                       size=size, 
                       a=lattice, 
                       vacuum=10)
        
        # Convert atoms.symbol -> atomic number
        constituent_atom_numbers = []
        for atom_symbol in constituent_atoms:
            atom_temp = ase.Atoms(atom_symbol)
            constituent_atom_numbers.append(atom_temp.get_atomic_numbers()[0])
        
        # Choose random constituent atoms sequence with np.random
        random_atoms_sequence = np.random.choice(constituent_atom_numbers, 
                                                 len(slab), 
                                                 replace=True).tolist()
        slab.set_atomic_numbers(random_atoms_sequence)
        
        binding_atom_symbol = get_binding_atom(slab)
        vasp_path = calculation_path / f"{binding_atom_symbol}_{calculation_count[binding_atom_symbol]}"
        calculation_count[binding_atom_symbol] += 1

        if constraint:
            slab = _set_constraint(slab)

        if adsorbate:
            adsorbate_atoms = read(adsorbate, index="-1")
            adslab = slab.copy()
            add_adsorbate(adslab, 
                          adsorbate_atoms, 
                          height=2, 
                          position=adsorbate_positions
                          )
            slab_save_path = vasp_path / "slab"
            adslab_save_path = vasp_path / "adslab"

            slab_save_path.mkdir(exist_ok=True, parents=True)
            adslab_save_path.mkdir(exist_ok=True, parents=True)

            write_vasp_inputs(slab, slab_save_path)
            write_vasp_inputs(adslab, adslab_save_path)
            logger.info(f"{vasp_path} calculation prepared!")

        if not adsorbate:
            slab_save_path = vasp_path / "slab"
            slab_save_path.mkdir(exist_ok=True, parents=True)
            write_vasp_inputs(slab, slab_save_path)
            logger.info(f"{vasp_path} calculation prepared!")
    

def generate_random_bulk(constituent_atoms: list[str],
                         size: tuple | list = (2, 2, 2),
    ) -> ase.Atoms:
    from ase.build import bulk
    
    from ..utils.globals import LATTICE_CONSTANTS
    r'''
    ===========================================================================
    Generate random high entropy alloy bulk

    Args:
        constituent_atoms: Atom symbols in high-entropy alloy, ["Co", "Ni", "Cu", "Ru"]
        size: repeat of slab, (2, 2, 2)

    Output:
        ase.Atoms

    ===========================================================================
    '''
    
    lattice = 0
    for constituent_atom in constituent_atoms:
        lattice += LATTICE_CONSTANTS[constituent_atom]
    lattice = lattice / len(constituent_atoms)
    
    constituent_atom_numbers = []
    for atom_symbol in constituent_atoms:
        atom_temp = ase.Atoms(atom_symbol)
        constituent_atom_numbers.append(atom_temp.get_atomic_numbers()[0])
    
    atoms = bulk("Cu",
                 "fcc",
                 a=lattice,
                 cubic=True)

    atoms = atoms.repeat(size)
    
    random_atom_sequence = np.random.choice(constituent_atom_numbers,
                                            len(atoms),
                                            replace=True).tolist()

    atoms.set_atomic_numbers(random_atom_sequence)

    return atoms

def generate_adslab_zpe(contcar_path: Path,
                        poscar_save_path: Path
    ) -> None:
    from ase.constraints import FixAtoms
    from ..VASP.write_vasp_inputs import write_vasp_inputs
    
    logger = setup_logger("HEA")

    adslab = read(contcar_path, index="-1")
    adslab.set_constraint()
    fix = [idx for idx, atom in enumerate(adslab) if atom.symbol not in ["N", "O"]]
    adslab.set_constraint(FixAtoms(fix))
    
    poscar_save_path.mkdir(exist_ok=True, parents=True)
    write_vasp_inputs(adslab, poscar_save_path, calculation_type="adslab_zpe")

    logger.info(f"generate_zpe_adslab on {poscar_save_path} done!")




