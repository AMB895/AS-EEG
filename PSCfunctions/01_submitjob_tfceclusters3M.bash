mkdir log/ -p
sbatch -J tfce_ersp_clusters_3M -o log/tfce_ersp_clusters_3M_o-%j.log -e log/tfce_ersp_clusters_3M_e-%j.log PSCfunctions/run_tfceERSPclusters3M.bash

