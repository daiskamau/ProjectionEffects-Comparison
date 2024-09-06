#!/bin/bash
#SBATCH -p bsudfq         # queue (partition) -- bsudfq, short
#SBATCH -t 03:00:00
#SBATCH -N 1    # number of nodes you want to run on
#SBATCH -n 48		  # number of cores on the node you want to use
#SBATCH -J Cardinal60         # job name
#SBATCH --mail-type=ALL
#SBATCH -o test60.sbatch.out


## source /global/common/software/lsst/common/miniconda/setup_current_python.sh
. ~/.bashrc
conda activate bsu

python /bsuhome/gladyskamau/BSU-Research/Cardinal/depth-richness.py "/bsuhome/gladyskamau/BSU-Research/Cardinal/Data/DepthData/FixedDepth-60/rlambda=2_035_05/" 60



##python /global/homes/k/kamau/SE-CLMM-LSSTDESC/cosmoDC2/DS_JOBS/py_files/parallelize-halos-calc-3.py  2>&1 | tee test.sbatch.out
## file_path1 = r"/bsuhome/gladyskamau/BSU-Research/Cardinal/Data/"
##SBATCH --qos debug
##SBATCH --qos regular
##SBATCH --time=00:30:00
##SBATCH --mem=4GB
## 124,125,126,127,129
##SBATCH --exclusive

## Activate the environment
## Replace environmentName with your environment name
##. ~/.bashrc
##conda activate environmentName
