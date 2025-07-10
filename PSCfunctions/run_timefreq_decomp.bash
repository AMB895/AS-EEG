#!/usr/bin/env bash

#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32

#SBATCH --time 50:00:00

echo  "[$(date)] STARTING JOB"
module load matlab
cd /ocean/projects/soc230004p/shared/antisaccade_eeg/code/
matlab -nodisplay -r "try,run('gen_ersp_timefreq_psc'),catch e, e, end; quit"
echo "[$(date)] STOP JOB"
