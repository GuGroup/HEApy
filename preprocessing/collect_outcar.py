"""
Collect OUTCAR and save as train.xyz / test.xyz / val.xyz
"""

from __future__ import annotations

import logging
from pathlib import Path

import re
import os
import math
import ase
import numpy as np
from ase.io import read, write
from tqdm import tqdm

logger = logging.getLogger("HEA")

def get_paths(calculation_path: Path
    ) -> list:
    r'''
    Get paths of relaxation results.

    Args:
        path: Path
        (The folder where relaxation results are have place)

    Returns:
        calculation_paths: list
        (List of each OUTCAR file path)

    '''
    calculation_paths = []

    for folder in calculation_path.glob("*"):
        if folder.is_dir():
            slab = folder / "slab" / "OUTCAR"
            adslab = folder / "adslab" / "OUTCAR"
            if slab.is_file() & adslab.is_file():
                calculation_paths.append([slab, adslab])

    logger.info(f"{len(calculation_paths)} calculations collected")
    
    return calculation_paths

def split_paths(calculation_paths: list
    ) -> tuple[list, list, list]:
    r'''
    Split the calculation paths into train/test/val
    Now split ratio is 90/5/5.

    Args:
        calculation_paths: list

    Returns:
        ttv_paths

    '''
    np.random.seed(1307)
    
    ttv_paths = {"train": [],
                 "test": [],
                 "val": [],
                 }

    for path in tqdm(calculation_paths):
        i = np.random.randint(100)

        if i <= 89:
            ttv_paths["train"].append(path)

        if i >= 90 and i<= 94:
            ttv_paths["test"].append(path)

        if i >= 95:
            ttv_paths["val"].append(path)

    return ttv_paths

def check_oszicar(path: Path
    ) -> list[bool]:
    r'''
    Check osizicar file in calculation folder, finding frame idxs which are not converged on ediff 60.

    Args:
        path: Path, path of calculation folder

    Returns:
        valid: list[bool], True for ediff converged fram and False for unconverged Frame ex) [T, T, T, F, F, ..]
    '''
    with open(path / 'OSZICAR') as f:
        oszicar = f.read()
    p = re.compile(r'^(.*)\n.*F=', re.MULTILINE)
    lines = p.findall(oszicar)
    ediffs = [float(line.split()[3]) for line in lines]
    valid = [True if math.fabs(ediff) < 10E-6 else False for ediff in ediffs]
    return valid

def check_validity(atoms: ase.Atoms,
                   force_threshold: float=5.0,
    ) -> bool:
    r'''
    Check validity of ase.Atoms for training
    If inf or NaN value in force or energy, 
    or force is over force_threshold eV/A,
    return False.
    (else True.)

    Args:
        atoms: ase.Atoms

    Output:
        bool

    '''
    # Force test
    forces = atoms.calc.get_forces()
    energy = atoms.calc.get_potential_energy()
    if True in np.isnan(forces) | True in np.isinf(forces):
        return False
    if math.isnan(energy) | math.isinf(energy):
        return False
    for force in forces:
        if np.linalg.norm(force) > force_threshold:
            return False

    return True

def collect_outcar(calculation_path: Path=Path("/home/ytk/MH/jobs"),
                   save_path: Path=Path("./"),
                   force_threshold: float=5.0
    ) -> None:
    import json 
    r'''
    Collect OUTCAR and save into
    train.xyz / test.xyz / val.xyz

    Args:
        calculation_path: Path
        save_path: Path

    Returns:
        None

    '''

    try:
        save_path.mkdir(parents=True, exist_ok=True)
    except:
        pass

    paths = get_paths(calculation_path)
    ttv_paths = split_paths(paths)

    save_log = {"train": {},
                "test": {},
                "val": {},
                }
    
    for key in save_log.keys():
        ttv_save_path = save_path / key
        ttv_save_path.mkdir(exist_ok=True)
        dataset = []
        idx = 0
        for path in tqdm(ttv_paths[key]):
            slab_path = path[0]
            adslab_path = path[1]
            try:
                slab = read(slab_path, index=":")
                adslab = read(adslab_path, index=":")
                save_log[key][str(slab_path)] = []
                save_log[key][str(adslab_path)] = []
                E_struc_converged_list_slab = check_oszicar(slab_path.parent)
                E_struc_converged_list_adslab = check_oszicar(adslab_path.parent)
            except:
                logger.debug(f"Failed to read {str(slab_path)} or {str(adslab_path)}")
                continue
            
            for atoms_idx, (atoms, E_struc_converged) in enumerate(zip(slab, E_struc_converged_list_slab)):
                validity = check_validity(atoms, force_threshold=force_threshold)
                if E_struc_converged is False:
                    logger.debug(f"Self Consistent Loop failed on {str(path)}, {atoms_idx}")
                if validity is False:
                    logger.debug(f"Validity failed on {str(path)}, {atoms_idx}")
                else:
                    dataset.append(atoms)
                    save_log[key][str(slab_path)].append(atoms_idx)
                    idx += 1

            for atoms_idx, (atoms, E_struc_converged) in enumerate(zip(adslab, E_struc_converged_list_adslab)):
                validity = check_validity(atoms, force_threshold=force_threshold)
                if E_struc_converged is False:
                    logger.debug(f"Self Consistent Loop failed on {str(path)}, {atoms_idx}")
                if validity is False:
                    logger.debug(f"Validity failed on {str(path)}, {atoms_idx}")
                else:
                    dataset.append(atoms)
                    save_log[key][str(adslab_path)].append(atoms_idx)
                    idx += 1

        write(ttv_save_path / (key + ".xyz"), dataset)
    
    with open(save_path / "split_log.json", "w") as f:
        json.dump(save_log, f)

    logger.info(f"Done! Datasets are saved on {str(save_path)}")

if __name__ == "__main__":
    collect_outcar()
