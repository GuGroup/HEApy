# ytk
# ytk

import json
import torch
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from tqdm import tqdm
from ase.io import read

from ..utils.logging import setup_logger
from ..core.dataset_utils import get_binding_atom
from ..tests.test_utils import get_gibbs_free_energy

logger = setup_logger("HEA")

def assemble_test_result(
    test_paths: Path,
    save_path: Path
    ) -> None:

    target = []
    prediction = []

    # Add data points
    for test_path in test_paths.glob("*"):
        if test_path.is_dir():
            slab_target = read(test_path / "slab" / "OUTCAR", index="-1")
            adslab_target = read(test_path / "adslab" / "OUTCAR", index="-1")
            
            binding_site_atom = get_binding_atom(read(test_path / "slab" / "POSCAR"))

            target_gibbs_free_energy = get_gibbs_free_energy(slab_target.calc.get_potential_energy(),
                                                             adslab_target.calc.get_potential_energy(),
                                                             hea_type="AlFeCoNiCu",
                                                             binding_site_atom=binding_site_atom
                                                             )

            slab_prediction = read(test_path / "slab" / "uma_prediction.traj", index="-1")
            adslab_prediction = read(test_path / "adslab" / "uma_prediction.traj", index="-1")

            prediction_gibbs_free_energy = get_gibbs_free_energy(slab_prediction.calc.get_potential_energy(),
                                                                 adslab_prediction.calc.get_potential_energy(),
                                                                 hea_type="AlFeCoNiCu",
                                                                 binding_site_atom=binding_site_atom
                                                                 )

            target.append(target_gibbs_free_energy)
            prediction.append(prediction_gibbs_free_energy)

    result = torch.tensor([target, prediction])
    torch.save(result, save_path)

def assemble_gibbs_free_energy_result(
    calculation_paths: Path,
    save_path: Path,
    file_name: str = "OUTCAR"
    ) -> None:

    result = {"total": [], "Al": [], "Fe": [], "Co": [], "Ni": [], "Cu": []}

    for calculation_path in tqdm(calculation_paths.glob("*")):
        if calculation_path.is_dir():
            slab = read(calculation_path / "slab" / file_name, index="-1")
            adslab = read(calculation_path / "adslab" / file_name, index="-1")
            
            binding_site_atom = get_binding_atom(read(calculation_path / "slab" / "POSCAR"))

            gibbs_free_energy = get_gibbs_free_energy(
                    slab.calc.get_potential_energy(),
                    adslab.calc.get_potential_energy(),
                    hea_type="AlFeCoNiCu",
                    binding_site_atom=binding_site_atom
            )
            
            result[binding_site_atom].append(gibbs_free_energy)
            result["total"].append(gibbs_free_energy)

    json.dump(result, open(save_path, "w"))

    logger.info(f"Assemble result from {calculation_paths} done!")

    


