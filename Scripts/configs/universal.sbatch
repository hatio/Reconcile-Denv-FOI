#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=15gb
#SBATCH --time=96:00:00
pwd; hostname; date
echo "slurm jobname:" ${SLURM_JOB_NAME}
echo "Prepping job:" ${SLURM_JOB_ID}



num_cpus=${SLURM_CPUS_ON_NODE:-1}
export SLURM_CPUS_ON_NODE=${num_cpus}

modules_to_load="${modules_to_load:-'R'}"
echo "loading modules:" "${modules_to_load}"
module load "${modules_to_load}"



echo -e ${cmd}

${cmd}

date
