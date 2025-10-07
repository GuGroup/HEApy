import subprocess
from pathlib import Path

def submit_a_job(job_path: Path,
                 core_64: bool = False,
    ) -> None:
    '''
    Submit a job in job_path
    
    Args:
        job_path: Path

    Returns:
        None
    '''
    qsub_file_path = Path(__file__) / "submit_vasp.sh"
    
    if core_64:
        qsub_file_path = Path(__file__) / "submit_vasp_64.sh"

    process = subprocess.Popen(
            ["qsub"] + str(qsub_file_path),
            cwd=job_path,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
    )

    for line in process.stdout:
        print(line, end="")
    
    process.wait()
    print(f"Job submitted in {job_path}")
