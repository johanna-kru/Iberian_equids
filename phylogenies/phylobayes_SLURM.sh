#!/bin/bash

#SBATCH --time=10-00:00
#SBATCH --job-name=pb_plio5
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=90G
#SBATCH --partition compute
#SBATCH --output=pliocene/log/array_%A-%a.out
#SBATCH --error=pliocene/log/array_%A-%a.err

# severeal chains started for method optimization, chain 4 and 5 were the final ones used for analysis.
# np 16 spread on two nodes iteration time of 11s, one one node 40 s, np 24 spread on two nodes iteration time of 40s, one one node 50 s 


module load mpi/openmpi/4.1.7-gnu-14.2

cd pliocene

mpirun -np 16 bin/pbmpi-1.9/data/pb_mpi -cat -gtr -d selected_pliocene.phylip -x 1 50000 tree_pliocene_5
