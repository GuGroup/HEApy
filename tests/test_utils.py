# mjs0000
# ytk

import numpy as np

def analyze_singlepoint_prediction_result(result: dict) -> None:
    """
    From the `result` dictionary, which is a return of `singlepoint_prediction`,
    Print energy/atom MAE, energ/atom RMSE, force MAE, force RMSE
    """
    target_e = np.array(result["target"]["energy_per_atom"])
    pred_e = np.array(result["prediction"]["energy_per_atom"])
    
    # Energy statistics
    energy_mae = np.mean(np.abs(target_e - pred_e))
    energy_rmse = np.sqrt(np.mean((target_e - pred_e)**2))
    
    # Force statistics
    target_fx = np.concatenate(result["target"]["force_x"])
    target_fy = np.concatenate(result["target"]["force_y"])
    target_fz = np.concatenate(result["target"]["force_z"])
    pred_fx = np.concatenate(result["prediction"]["force_x"])
    pred_fy = np.concatenate(result["prediction"]["force_y"])
    pred_fz = np.concatenate(result["prediction"]["force_z"])
    
    force_mae = np.mean([
        np.abs(target_fx - pred_fx),
        np.abs(target_fy - pred_fy),
        np.abs(target_fz - pred_fz)
    ])  
    
    force_rmse = np.sqrt(np.mean([
        (target_fx - pred_fx)**2,
        (target_fy - pred_fy)**2,
        (target_fz - pred_fz)**2
    ])) 
    
    print("\n" + "="*50)
    print("EVALUATION METRICS")
    print("="*50)
    print(f"Energy/Atom MAE: {energy_mae*1000:.3f} meV/atom")
    print(f"Energy/Atom RMSE: {energy_rmse*1000:.3f} meV/atom")
    print(f"Force MAE: {force_mae*1000:.3f} meV/Å")
    print(f"Force RMSE: {force_rmse*1000:.3f} meV/Å")
    print("="*50 + "\n")

def get_gibbs_free_energy(energy_slab: float,
                          energy_adslab: float,
                          hea_type: str,
                          binding_site_atom: str
    ) -> float:
    from ..utils.globals import POTENTIAL_ENERGY, ZPE, TS
    r'''
    ===========================================================================
    Calculate adsorption energy. (It might be called as relative energy..)
    Formula is,
    E_ads = E_adslab - E_slab - E_NO
    
    Args:
        adslab: ase.Atoms
        slab: ase.Atoms

    Returns:
        labeled_site: dict
    ===========================================================================
    '''
    adsorption_energy = energy_adslab + 2 * POTENTIAL_ENERGY["H2O"] - (energy_slab + 1.5 * POTENTIAL_ENERGY["H2"] + POTENTIAL_ENERGY["HNO3"])
    zero_point_energy = ZPE[hea_type][binding_site_atom] + 2 * ZPE["H2O"] - (1.5 * ZPE["H2"] + ZPE["HNO3"])
    temperature_entropy = TS[hea_type][binding_site_atom] + 2 * TS["H2O"] - (1.5 * TS["H2"] + TS["HNO3"]) 
    
    gibbs_free_energy = adsorption_energy + zero_point_energy - temperature_entropy - 1.12
    
    return gibbs_free_energy
