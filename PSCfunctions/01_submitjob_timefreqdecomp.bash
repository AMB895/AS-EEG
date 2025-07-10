mkdir log/ -p
sbatch -J ersp_timefreq_decomp -o log/ersp_timefreq_o-%j.log -e log/ersp_timefreq_e-%j.log run_timefreq_decomp.bash
