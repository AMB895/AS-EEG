mkdir log/ -p
currentdate = $(date + %Y%m%d)
sbatch -J tfce_ersp_clusters -o log/tfce_ersp_clusters_o_${currentdate}-%j.log -e log/tfce_ersp_clusters_e_${currentdate}-%j.log PSCfunctions/run_tfceERSPclusters.bash
