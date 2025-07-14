mkdir log/ -p
sbatch -J tfce_ersp_clusters -o log/tfce_ersp_clusters_o-%j.log -e log/tfce_ersp_clusters_e-%j.log PSCfunctions/run_tfceERSPclusters.bash
