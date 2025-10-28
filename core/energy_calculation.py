# ytk
# ytk

import re
from pathlib import Path

from ..utils.logging import setup_logger

logger = setup_logger("HEA")

def calculate_zpe(calculation_path: Path,
    ) -> float:
    r"""
    calculate zpe based on freqs_THz
    grep THz from OUTCAR, pass i values.
    """
    freqs_THz = []

    with open(calculation_path / "OUTCAR", "r") as f:
        for line in f:
            match = re.search(r'([\d.]+)\s+THz', line)
            if match:
                freq = float(match.group(1))
                if "f/i" not in line and freq > 0:
                    freqs_THz.append(freq)

    if not freqs_THz:
        logging.warning(f"No freqency information in {calcualtion_path}. Check it")
        return 0.0

    freqs_eV = [f * 0.0041357 for f in freqs_THz]
    zpe = 0.5 * sum(f * 0.0041357 for f in freqs_THz)
    logger.info(f"{calculation_path} on ZPE: {zpe:.5f} eV")

    return zpe

def calculate_entropy(
        calculation_path: Path,
        atoms_type: str = "adslab"
    ) -> float:
    r"""
    Calculate entropy on 300K with freqs_THz
    """
    from ase.io import read
    from ase.thermochemistry import IdealGasThermo

    
    freqs_THz = []

    with open(calculation_path / "OUTCAR", "r") as f:
        for line in f:
            match = re.search(r'([\d.]+)\s+THz', line)
            if match:
                freq = float(match.group(1))
                if "f/i" not in line and freq > 0:
                    freqs_THz.append(freq)
    
    contcar = read(calculation_path / "CONTCAR", index="-1")
    freqs_eV = [f * 0.0041357 for f in freqs_THz]
    
    #put 2 on symmetrynumber if you calculate new gas
    thermo = IdealGasThermo(
            vib_energies=freqs_eV,
            potentialenergy=0.0,
            atoms=contcar,
            geometry="linear",
            symmetrynumber=1,
            spin=0
            )
    T = 300
    P = 101325
    S_vib = thermo.get_entropy(T, P)
    logger.info(f"{calculation_path} on Vibrational_entropy: {S_vib:.5f} eV on {T}K, {P} pressure")

    return S_vib

