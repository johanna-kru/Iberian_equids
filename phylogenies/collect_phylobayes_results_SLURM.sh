#!/bin/bash

#SBATCH --time=00-02:00
#SBATCH --job-name=bp_trace
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --partition compute
#SBATCH --output=log/array_%A-%a.out
#SBATCH --error=log/array_%A-%a.err

# severeal chains started for method optimization, chain 4 and 5 were the final ones used for analysis.
# np 16 spread on two nodes iteration time of 11s, one one node 40 s, np 24 spread on two nodes iteration time of 40s, one one node 50 s 


module load mpi/openmpi/4.1.7-gnu-14.2

cd phylobayes/pliocene

mpirun -np 4 bin/pbmpi-1.9/data/tracecomp -x 1000 tree_pliocene_4 tree_pliocene_5 > tracecomp_tree_pliocene_4_5.txt

mpirun -np 4 bin/pbmpi-1.9/data/bpcomp -x 1000 10 tree_pliocene_4 tree_pliocene_5 > bpcomp_tree_pliocene_4_5.txt