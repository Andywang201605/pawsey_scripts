#!/bin/usr/env python
# functions to make dspsr files for each epoch

import glob
import os

import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)

logger = logging.getLogger(__name__)

def copy_par_file(date, ):
    pardir = "/scratch/pawsey0272/zwan4817/J1646/SCRIPTS/dspsr/J1646.ell.par"
    dspsrdir = f"/scratch/pawsey0272/zwan4817/J1646/preprocess/{date}/dspsr"

    cmd = f"cp {pardir} {dspsrdir}"
    logger.info(f"Copying par file to {dspsrdir}")
    os.system(cmd)
    
def make_sbatch_file(date):
    dspsrdir = f"/scratch/pawsey0272/zwan4817/J1646/preprocess/{date}/dspsr"

    ### find sf file...
    files = glob.glob(f"/scratch/pawsey0272/zwan4817/J1646/preprocess/{date}/add/J1646_p*.sf")
    if len(files) == 0:
        logger.error(f"No sf files found for date {date}")
        return
    sffile = files[0]

    ### make sbatch file
    cmd = f"""#!/bin/bash
#SBATCH --account=pawsey0272
#SBATCH --output=./logs/%x_%A_%a.out
#SBATCH --error=./logs/%x_%A_%a.err
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --job-name=hdspsr{date[2:]}
#SBATCH --cpus-per-task=12



module use /software/projects/pawsey0272/setonix/2025.08/modules/zen3/gcc/14.2.0
module load singularity/4.1.0-psr

export psrsrc="singularity exec /software/projects/pawsey0272/zwan4817/container/psr-search.sif"
export psraly="singularity exec /software/projects/pawsey0272/zwan4817/container/psr-analysis.sif" 

$psraly dspsr -E J1646.ell.par -scloffs -U 4096 -t 12 -K -r -b 32 -A -L 30 -nsub 10 -O J1646.{date}.30s.ddm {sffile}
$psraly psradd -v -o J1646.{date}.30s.ddm.ar J1646.{date}.30s.ddm_*.ar
$psraly pdmp -pr 0.005 -ps 0.00001 -dr 5 -ds 0.1 -mc 256 -g pdmp.1s.png/png J1646.{date}.30s.ddm.ar
"""

    sbatch_file = f"{dspsrdir}/dspsr_{date}.30s.sbatch"
    with open(sbatch_file, "w") as f:
        f.write(cmd)
    logger.info(f"Created sbatch file: {sbatch_file}")

if __name__ == "__main__":
    submit_scripts = "#!/bin/bash\n"

    epochfolders = glob.glob("/scratch/pawsey0272/zwan4817/J1646/preprocess/20*")
    for folder in epochfolders:
        date = os.path.basename(folder)
        if not date.startswith("2024"):  # only process 2024 dates
            continue
        # if date == "20251230": continue
        dspsrdir = f"/scratch/pawsey0272/zwan4817/J1646/preprocess/{date}/dspsr"
        copy_par_file(date)
        make_sbatch_file(date)
        submit_scripts += f"#--------------{date}------------------\n"
        submit_scripts += f"cd {dspsrdir}\n"
        submit_scripts += f"sbatch dspsr_{date}.30s.sbatch\n\n"

    with open("submit_dspsr.sh", "w") as f:
        f.write(submit_scripts)
