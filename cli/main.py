# ytk
# ytk

from pathlib import Path

def run_assemble_result():
    from ..draw.assemble_result import assemble_gibbs_free_energy_result
    
    dft_paths = Path("/home/ytk/HAE/DFTs/CoNiCuRu/slab_224")
    dft_save_path = Path("/home/ytk/HAE/DFTs/CoNiCuRu/results/dft_result.json")

    assemble_gibbs_free_energy_result(
        calculation_paths=dft_paths,
        save_path=dft_save_path,
        file_name="OUTCAR"
    )
    

    uma_paths = Path("/home/ytk/HAE/DFTs/CoNiCuRu/prediction_224")
    uma_save_path = Path("/home/ytk/HAE/DFTs/CoNiCuRu/results/uma_result.json")
    assemble_gibbs_free_energy_result(
        calculation_paths=uma_paths,
        save_path=uma_save_path,
        file_name="uma_prediction.traj"
    )

def run_draw():
    from ..draw.figures import (draw_parity_plot, 
                                draw_gibbs_free_energy_distribution,
                                draw_volcano_plot
    )

    dft_result = Path("/home/ytk/HAE/DFTs/CoNiCuRu/results/dft_result.json")
    uma_result = Path("/home/ytk/HAE/DFTs/CoNiCuRu/results/uma_result.json")
    draw_gibbs_free_energy_distribution(
        dft_result_path=dft_result,
        uma_result_path=uma_result
    )
    #draw_volcano_plot(uma_result_path=uma_result)

def run_prediction():
    from ..tests.tests import relaxation_prediction
    from ..utils.globals import UMA_CKPT_PATH

    calculation_paths = Path("/home/ytk/HAE/DFTs/CoNiCuRu/prediction_224")
    
    for calculation_path in calculation_paths.glob("*"):
        slab_path = calculation_path / "slab"
        adslab_path = calculation_path / "adslab"

        relaxation_prediction(ckpt_path=UMA_CKPT_PATH,
                              calculation_path=slab_path,
        )
        relaxation_prediction(ckpt_path=UMA_CKPT_PATH,
                              calculation_path=adslab_path,
        )

def run_evaluation():
    import json

    from ..tests.tests import singlepoint_prediction, relaxation_prediction
    from ..utils.globals import UMA_CKPT_PATH
    
    calculation_path = Path("/home/ytk/HAE/DFTs/CoNiCuRu/slab_224/Cu_38/adslab/")


    '''
    Example of singlepoint_prediction
    ---------------------------------
    singlepoint_prediction(
           ckpt_path=UMA_CKPT_PATH,
           dataset_path=Path("/home/ytk/HAE/DFTs/CoNiCuRu/dataset_224/test/test.xyz"),
           result_path=Path("/home/ytk/HAE/DFTs/CoNiCuRu/dataset_224/singlepoint_prediction_result.pt")
    )
    
    Example of relaxation_prediction
    --------------------------------
    split_log = json.load(open("/home/ytk/HAE/DFTs/CoNiCuRu/dataset_224/split_log.json", "r"))
    
    for calculation_path in split_log["test"].keys():
        relaxation_prediction(ckpt_path=UMA_CKPT_PATH,
                              calculation_path=Path(calculation_path).parent,
        )

    '''

def run_slab_generation():
    from ..utils.logging import setup_logger
    from ..core.generate_atoms import generate_random_slab
    from ..VASP.write_vasp_inputs import write_vasp_inputs
    
    adsorbate_path = Path(__file__).parents[1] / "VASP" / "NO_CONTCAR"
    calculation_path = Path("/home/ytk/HAE/DFTs/CoNiCuRu/prediction_224") 
    generate_random_slab(constituent_atoms=["Al", "Fe", "Co", "Ni", "Cu"],
                         calculation_path=calculation_path,
                         size=(2, 2, 4),
                         constraint=True,
                         adsorbate=adsorbate_path,
                         adsorbate_positions=(0, 0),
                         count=2500,
                         )

    """
    Example of generate_random_slab
    -------------------------------

    adsorbate_path = Path(__file__).parents[1] / "VASP" / "NO_CONTCAR"
    calculation_path = Path("/home/ytk/HAE/DFTs/CoNiCuRu/slab_336") 
    generate_random_slab(constituent_atoms=["Co", "Ni", "Cu", "Ru"],
                         calculation_path=calculation_path,
                         size=(3, 3, 6),
                         constraint=True,
                         adsorbate=None,
                         adsorbate_positions=None,
                         count=500,
                         )
    """    

def run_bulk_generation():
    from ..utils.logging import setup_logger
    from ..core.generate_atoms import generate_random_bulk
    from ..VASP.write_vasp_inputs import write_vasp_inputs

    logger = setup_logger("HEA")
    calculation_path = Path("/home/ytk/hea_test2")

    for i in range(2):
        vasp_path = calculation_path / f"calculation_{i}"
        vasp_path.mkdir(exist_ok=True, parents=True)
        
        bulk = generate_random_bulk(constituent_atoms=["Co", "Ni", "Cu", "Ru"],
                                    size=(2, 2, 2),
                                    )
        
        vasp_path.mkdir(exist_ok=True, parents=True)
        
        write_vasp_inputs(bulk, vasp_path, calculation_type="bulk")
        logger.info(f"calculation_{i} saved in {calculation_path}")


def run_adslab_zpe_generation():
    from ..utils.logging import setup_logger
    from ..core.generate_atoms import generate_adslab_zpe
    
    logger = setup_logger("HEA")
    
    """
    Example of generate_adslab_zpe
    ------------------------------
    Co_0 = Path("/home/ytk/HAE/DFTs/CoNiCuRu/slab_224/Al_0/adslab/CONTCAR")
    Co_save_path = Path("/home/ytk/HAE/DFTs/CoNiCuRu/ZPE/")
    generate_adslab_zpe(Co_0, Co_save_path)

    """
def run_zpe_calculation():
    from ..core.energy_calculation import calculate_zpe, calculate_entropy
    
    """
    Example of calculate_zpe and calculate_entropy
    Co_save_path = Path("/home/ytk/HAE/DFTs/CoNiCuRu/ZPE/Co")
    calculate_zpe(Co_save_path)
    calculate_entropy(Co_save_path)
    """

def run_job_submit():
    from ..VASP.submit_job import submit_a_job
    
    '''
    calculation_path = Path("/home/ytk/HAE/DFTs/CoNiCuRu/slab_336")
    for vasp_path in calculation_path.glob("*"):
        submit_a_job(vasp_path / "slab", core_64=True)
    '''

def run_preprocessing():
    from pathlib import Path
    from ..preprocessing import collect_outcar
    
    """
    Example of collect_outcar
    -------------------------
    collect_outcar(calculation_path=Path("/home/ytk/HAE/DFTs/CoNiCuRu/slab_224"),
                   save_path=Path("/home/ytk/HAE/DFTs/CoNiCuRu/dataset_224"),
                   )
    """

if __name__ == "__main__":
    #run_slab_generation()
    #run_preprocessing()
    #run_adslab_zpe_generation()
    #run_job_submit()
    #run_zpe_calculation()
    #run_evaluation()
    run_draw()
    #run_prediction()
    #run_assemble_result()

