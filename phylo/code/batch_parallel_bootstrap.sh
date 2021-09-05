#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=2gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=rax-bs



##### load dependencies #####
module load raxmlng

##### Assign vars #####
SEED=$1$1
THREADS=$2
MODEL=$3
ALN=$4
ODIR=$5


##### run raxml-ng bootstrap #####

for i in $(seq 1 3)
do
raxml-ng --bootstrap --msa $ALN --model $MODEL --prefix "$ODIR/bs${1}n${i}" \
--threads $THREADS --seed $i$SEED --tree pars{2} --bs-trees 1  &
done

for i in $(seq 4 6)
do
raxml-ng --bootstrap --msa $ALN --model $MODEL --prefix "$ODIR/bs${1}n${i}" \
--threads $THREADS --seed $i$SEED --tree rand{2} --bs-trees 1  &
done

wait

