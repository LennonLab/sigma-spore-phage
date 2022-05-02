#------------------------------
# bootstrapping with  RAxML-NG
#------------------------------


#Interactive job on Carbonate
srun -p interactive -N 1 --ntasks-per-node=1 --cpus-per-task=8 --time=07:59:00 --pty bash

#### load dependencies ####

module load raxmlng
  # Flex 2.6.4 loaded.
  # CMake version 3.8.0 loaded.
  # raxmlng version 0.9.0-pthreads loaded.

##### Define paths #####
PARENT=~/GitHub/sigma-spore-phage/phylo-clust

ODIR=${PARENT}/data/align-trim-tree/tree/bootstraps
mkdir -p $ODIR
cd $ODIR

ALN=$PARENT/data/align-trim-tree/check_msa/parse-msa.raxml.rba

#modeltest-ng selected model : 
MODEL="LG+G4"

# From check_msa:
# * Recommended number of threads / MPI processes: 4
threads=4

# parallelized bootstraps

### Vars in batch_bootstrap.sh
# SEED=$1$1$1 (set by loop)
# THREADS=$2 (set above)
# TREES=$3 # type{number}
# trees="pars{5},rand{5}"
# MODEL=$4 (set above)
# ALN=$5 (set above)
# ODIR=$6 (set above)
# BOOTS=$7
boots=5
 

# 500 bootstrap trees with parsimony starts (100 runs x 5 per run)
seeds=($(seq 1 100))
trees="pars{5}"
for i in ${seeds[@]}; do
sbatch --job-name="bs$i" --time=4:59:00 --cpus-per-task=4 \
$PARENT/code/batch_bootstrap.sh $i $threads $trees $MODEL $ALN $ODIR $boots
done

# 500 bootstrap trees with random starts (100 runs x 5 per run)
seeds=($(seq 101 200))
trees="rand{5}"
for i in ${seeds[@]}; do
sbatch --job-name="bs$i" --time=4:59:00 --cpus-per-task=4 \
$PARENT/code/batch_bootstrap.sh $i $threads $trees $MODEL $ALN $ODIR $boots
done

# # while waiting for those trees to compute I am trying 
# # to run a more parellized script suggested by Kozolov&Stamatakis
# # this should compute 6 trees simultaneously: 3 random start and 3 parsimony start
# ## Vars to batch script
# # SEED=$1$1$1
# # THREADS=$2
# # MODEL=$3
# # ALN=$4
# # ODIR=$5
# seeds=($(seq 105 105))
# for i in ${seeds[@]}; do
# sbatch --job-name="bsP$i" --time=0:09:00 --cpus-per-task=4 \
# $PARENT/code/batch_parallel_bootstrap.sh $i $threads $MODEL $ALN $ODIR 
# done


#check for convergence
cat *.bootstraps > allbootstraps
 wc -l allbootstraps
raxml-ng --bsconverge --bs-trees allbootstraps \
--prefix check_bs_conv --seed 123 --threads 4 --bs-cutoff 0.03

# Bootstopping test converged after 450 trees!!!

# Consensus tree building
# raxml-ng --consense MRE --tree allbootstraps --prefix consMRE
# # this yielded a tree with multifurications :(

#### compute support for best scoring ML tree ####

cd ..

raxml-ng --support \
--tree bestScoreML.tree \
--bs-trees bootstraps/allbootstraps \
--prefix ML_TBE_tree --threads 4 --bs-metric tbe,fbp 



# # #reruns to complete time out
# seeds=( $( ls $ODIR |grep "ckp" | sed 's/n.*//' | egrep -o [0-9]+ | sort -u ) )
# echo ${seeds[@]}
# # # seeds=($(seq 55 57))
# # # seeds=( 59 60 62 65 68 69 70 )
# seeds=( 101 102  )
# trees="rand{5}"
# 
# seeds=( 100 11 48 49 50 51 79 98 99 )
# trees="pars{5}"
# 
# 
# for i in ${seeds[@]}; do
# sbatch --job-name="bs$i" --time=1:59:00 --cpus-per-task=4 \
# $PARENT/code/batch_bootstrap.sh $i $threads $trees $MODEL $ALN $ODIR $boots
# done
# squeue -u danschw 
