# ytk/ktg
# ytk/ktg

import ase
import numpy as np

from pathlib import Path
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.build import fcc111, add_adsorbate

def generate_random_slab(constituent_atoms: list[str],
                         size: tuple | list,
                         constraint: bool | None,
                         adsorbate: Path | None,
                         adsorbate_positions: tuple | None,
    ) -> ase.Atoms:
    r'''
    ===========================================================================
    Generate random high entropy alloy slab (+ adsorbate)
    Half bottom layers are fixed if constraint

    Args:
        constituent_atoms: Atom symbols in high-entropy alloy, ["Co", "Ni", "Cu", "Ru"]
        size: size of slab, (2, 2, 4) or (3, 3, 6)
        constraint: whether fix atoms on low half
        adsorbate: the path of adsorbate
        adsorption_positions: adsorption (x, y) positions, height is set on 2

    Output:
        ase.Atoms

    ===========================================================================
    '''
    from .dataset_utils import _set_constraint
    from ..utils.globals import LATTICE_CONSTANTS
    
    lattice = 0
    for constituent_atom in constituent_atoms:
        lattice += LATTICE_CONSTANTS[constituent_atom]
    lattice = lattice / len(constituent_atoms)

    atoms = fcc111("Cu", 
                   size=size, 
                   a=lattice, 
                   vacuum=10)

    constituent_atom_numbers = []
    for atom_symbol in constituent_atoms:
        atom_temp = ase.Atoms(atom_symbol)
        constituent_atom_numbers.append(atom_temp.get_atomic_numbers()[0])
    
    random_atoms_sequence = np.random.choice(constituent_atom_numbers, 
                                             len(atoms), 
                                             replace=True).tolist()
    
    atoms.set_atomic_numbers(random_atoms_sequence)
    
    if constraint:
        atoms = _set_constraint(atoms)

    if adsorbate:
        adsorbate_atoms = read(adsorbate, index="-1")
        adslab = atoms.copy()
        add_adsorbate(adslab, 
                      adsorbate_atoms, 
                      height=2, 
                      position=adsorbate_positions
                      )
        return atoms, adslab

    return atoms

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


