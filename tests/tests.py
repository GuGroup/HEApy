'''
Two test method, singlepoint / full relxation

'''
from __future__ import annotations

import logging
from pathlib import Path

import torch
import json
import ase
import os
from ase.io import read, write
from tqdm import tqdm
from fairchem.core.units.mlip_unit import load_predict_unit
from fairchem.core import pretrained_mlip, FAIRChemCalculator

from HEApy.tests.test_utils import get_gibbs_free_energy
from ..utils.logging import setup_logger

logger = setup_logger("HEA")

def singlepoint_prediction(ckpt_path: Path,
                           dataset_path: Path,
                           result_path: Path,
    ) -> dict:
    from ase.optimize import BFGS
    from datetime import datetime
    
    predictor = load_predict_unit(ckpt_path, device="cuda")
    calc = FAIRChemCalculator(predictor, task_name="oc20")

    dataset = read(dataset_path, index=":")

    result = {"target": {"total_energy": [],
                         "force_x": [],
                         "force_y": [],
                         "force_z": []},

              "prediction": {"total_energy": [],
                             "force_x": [],
                             "force_y": [],
                             "force_z": []}}

    for atoms in tqdm(dataset, desc="Evaluating structures"):
        target_energy = atoms.calc.get_potential_energy()
        target_forces = atoms.calc.get_forces()
        target_forces_x = target_forces[:, 0].tolist()
        target_forces_y = target_forces[:, 1].tolist()
        target_forces_z = target_forces[:, 2].tolist()

        atoms_copy = atoms.copy()
        atoms_copy.calc = calc
        
        opt = BFGS(atoms_copy)
        opt.run(fmax=0.05, steps=0)

        prediction_energy = atoms_copy.calc.get_potential_energy()
        prediction_forces = atoms_copy.calc.get_forces()
        prediction_forces_x = prediction_forces[:, 0].tolist()
        prediction_forces_y = prediction_forces[:, 1].tolist()
        prediction_forces_z = prediction_forces[:, 2].tolist()

        result["target"]["total_energy"].append(target_energy)
        result["target"]["force_x"].append(target_forces_x)
        result["target"]["force_y"].append(target_forces_y)
        result["target"]["force_z"].append(target_forces_z)

        result["prediction"]["total_energy"].append(prediction_energy)
        result["prediction"]["force_x"].append(prediction_forces_x)
        result["prediction"]["force_y"].append(prediction_forces_y)
        result["prediction"]["force_z"].append(prediction_forces_z)
   
    
    result_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(result, result_path)
    
    logger.info("Finished singlepoint prediction")

    return result


def relaxation_prediction(ckpt_path: Path,
                          calculation_path: Path,
                          file_name: str = "uma_prediction.traj",
                          steps: int=400
    ) -> dict:
    from ase.optimize import BFGS
    
    predictor = load_predict_unit(ckpt_path, device="cuda")
    calc = FAIRChemCalculator(predictor, task_name="oc20")
    poscar = read(calculation_path / "POSCAR")
    poscar.calc = calc
    opt = BFGS(poscar, trajectory=str(calculation_path / file_name))
    opt.run(fmax=0.05, steps=steps)  
    

    logger.info(f"Finished relaxation prediction on {calculation_path}")

