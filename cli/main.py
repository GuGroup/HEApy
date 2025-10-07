# ytk
# ytk

from pathlib import Path

def run_VASP_preparation():
    from ..utils.logging import setup_logger
    from ..core.generate_slab import generate_random_slab
    from ..VASP.write_vasp_inputs import write_vasp_inputs

    logger = setup_logger("HEA")
    calculation_path = Path("/home/ytk/hea_test")

    for i in range(2):
        vasp_path = calculation_path / f"calculation_{i}"
        vasp_path.mkdir(exist_ok=True, parents=True)
        adsorbate_path = Path(__file__).parents[1] / "VASP" / "NO_CONTCAR"
        
        slab, adslab = generate_random_slab(constituent_atoms=["Co", "Ni", "Cu", "Ru"],
                                            size=(3, 3, 6),
                                            constraint=True,
                                            adsorbate=adsorbate_path,
                                            adsorbate_positions=(0, 0)
                                            )
        
        slab_save_path = vasp_path / "slab"
        adslab_save_path = vasp_path / "adslab"
        
        slab_save_path.mkdir(exist_ok=True, parents=True)
        adslab_save_path.mkdir(exist_ok=True, parents=True)
        
        write_vasp_inputs(slab, slab_save_path)
        write_vasp_inputs(adslab, adslab_save_path)
        logger.info(f"calculation_{i} saved in {calculation_path}")

if __name__ == "__main__":
    run_VASP_preparation()


