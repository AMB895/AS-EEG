#!/usr/bin/env bash

#SBATCH --partition=RM
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128

#SBATCH --time 17:00:00
echo "[$(date)] STARTING JOB"
module load matlab
export NUMBER_OF_PROCESSORS=128
cd /ocean/projects/soc230004p/shared/antisaccade_eeg/
matlab -nodisplay -r "try,run('permute_tfce'),catch e, e, end; quit"
echo "[$(date)] STOP JOB"
