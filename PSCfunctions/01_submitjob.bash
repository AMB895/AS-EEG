mkdir log/ -p
sbatch -J run_permutation -o log/permutation-%j.log -e log/permutation-%j.log run_permutations.bash
