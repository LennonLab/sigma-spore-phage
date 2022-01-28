#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1:30:00
#SBATCH --mem=1gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=rax-bs

##### load dependencies #####
module load raxmlng

##### Assign vars #####
SEED=$1$1$1
THREADS=$2
TREES=$3 # type{number}
MODEL=$4
ALN=$5
ODIR=$6
BOOTS=$7

##### run raxml-ng bootstrap #####

raxml-ng --bootstrap --msa $ALN --model $MODEL  --prefix "$ODIR/bs$1" \
--seed $SEED --threads $THREADS \
--tree $TREES --bs-trees $BOOTS