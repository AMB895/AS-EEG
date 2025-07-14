#!/usr/bin/env bash
#
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#
#SBATCH --time 00:07:00
#
# 20250714WF - init
#
echo "[$(date)] STARTING JOB"
module load matlab
cd /ocean/projects/soc230004p/shared/antisaccade_eeg/code/
matlab -nodisplay -r "try,run('tfce_clusters_ersp_psc'),catch e, e, end; quit"
echo "[$(date)] STOP JOB"
