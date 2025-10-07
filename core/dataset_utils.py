import ase
from pathlib import Path

def check_validity(atoms: ase.Atoms
    ) -> bool:
    import numpy as np
    import math
    r'''
    Check whether ase.Atoms can be a data of dataset.
    If there is a value of NaN or inf,
    or if any atom have force more than 5 eV / A
    it fail.

    Args:
        atoms: ase.Atoms

    Returns:
        bool

    Example:
        valid = check_validity(atoms)

        if valid:
            ...
    '''
    forces = atoms.calc.get_forces()
    #forces = torch.tensor(forces, dtype=torch.float32)
    energy = atoms.calc.get_potential_energy()
    if True in np.isnan(forces) | True in np.isinf(forces):
        return False
    if math.isnan(energy) | math.isinf(energy):
        return False
    for force in forces:
        if np.linalg.norm(force) > 5:
            return False

    return True

