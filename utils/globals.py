# ytk
# ytk

from pathlib import Path

UMA_CKPT_PATH = Path("/home/ytk/HAE/DFTs/CoNiCuRu/dataset_224/uma_fine_tuning/202510-2318-0020-853f/checkpoints/final/inference_ckpt.pt")

COLORMAP = {"total": "grey",
            "Co": "lightcoral",
            "Ni": "lime",
            "Cu": "darkgoldenrod",
            "Ru": "teal",
            }


LATTICE_CONSTANTS = {
                     "Al": 4.049,
                     "Au": 4.189,
                     "Co": 3.549,
                     "Cr": 3.633,
                     "Cu": 3.676,
                     "Fe": 3.672,
                     "Hf": 4.504,
                     "Ir": 3.887,
                     "Lu": 4.814,
                     "Mn": 3.518,
                     "Mo": 4.010,
                     "Nb": 4.259,
                     "Ni": 3.552,
                     "Os": 3.873,
                     "Pd": 3.974,
                     "Pt": 3.990,
                     "Re": 3.926,
                     "Rh": 3.847,
                     "Ru": 3.816,
                     "Sc": 4.643,
                     "Ta": 4.227,
                     "Tc": 3.867,
                     "Ti": 4.136,
                     "V": 3.820,
                     "W": 4.033,
                     "Zr": 4.549}

MAGMOMS = {"Al": 0.0,
           "Au": 0.0,
           "Co": 2.0,
           "Cr": 3.0,
           "Cu": 1.0,
           "Fe": 3.0,
           "Hf": 0.0,
           "Ir": 0.0,
           "Lu": 0.0,
           "Mn": 3.0,
           "Mo": 0.0,
           "Nb": 1.0,
           "Ni": 1.0,
           "Os": 0.0,
           "Pd": 0.0,
           "Pt": 0.0,
           "Re": 0.0,
           "Rh": 0.0,
           "Ru": 0.0,
           "Sc": 1.0,
           "Ta": 0.0,
           "Tc": 2.0,
           "Ti": 0.0,
           "V": 1.0,
           "W": 0.0,
           "Zr": 0.0}

POTENTIAL_ENERGY = {
        "H2": -6.97736696,
        "H2O": -14.15696269,
        "HNO3": -27.693,
        }

ZPE = {
       "CoNiCuRu": {
           "Co": 0.1220,
           "Ni": 0.1309,
           "Cu": 0.1306,
           "Ru": 0.1939,
        },
       "H2": 0.277,
       "H2O": 0.586,
       "HNO3": 0.703
       }

# on 300K, only vib!
TS = {
      "CoNiCuRu":{
          "Co": 0.177,
          "Ni": 0.258,
          "Cu": 0.125,
          "Ru": 0.099,
      },
      "H2": 0.407,
      "H2O": 0.578,
      "HNO3": 0.781
      }

