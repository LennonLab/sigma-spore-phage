#--------------------------
# phylogeny with  RAxML-NG
#--------------------------


#Interactive job on Carbonate
srun -p interactive -N 1 --ntasks-per-node=1 --cpus-per-task=8 --time=07:59:00 --pty bash

#### load dependencies ####

module load raxmlng
  # Flex 2.6.4 loaded.
  # CMake version 3.8.0 loaded.
  # raxmlng version 0.9.0-pthreads loaded.

##### Define paths #####
PARENT=~/GitHub/sigma-spore-phage/phylo
ODIR=${PARENT}/data/align-trim-tree/check_msa
mkdir -p $ODIR
ALN=$PARENT/data/align-trim-tree/sigmas_MafftEinsi.trim

#modeltest-ng selected model : 
MODEL="LG+G4"
cd $ODIR

raxml-ng --check --msa $ALN  \
--model $MODEL --data-type AA \
--prefix check-msa 

# WARNING: Duplicate sequences found: 41
# 
# NOTE: Reduced alignment (with duplicates and gap-only sites/taxa removed) 
# NOTE: was saved to: /geode2/home/u020/danschw/Carbonate/GitHub/sigma-spore-phage/phylo/data/align-trim-tree/check_msa/check-msa.raxml.reduced.phy

# checked if sequences removed as duplicates do not include sequences directly discussed in MS
# duplicated-seqs.R
# swapped MSA headers as needed

ALN=$ODIR/check-msa.raxml.reduced.phy

# For large alignments, we recommend using the --parse command after, or, instead of
# --check:
raxml-ng --parse --msa $ALN --data-type AA --model $MODEL --prefix parse-msa 

# Binary MSA file created: parse-msa.raxml.rba
ALN=$ODIR/parse-msa.raxml.rba

# * Estimated memory requirements                : 153 MB
# 
# * Recommended number of threads / MPI processes: 4

# get tree
ODIR=$PARENT/data/align-trim-tree/tree
mkdir -p $ODIR
cd $ODIR
# raxml-ng --msa $ALN --model $MODEL --threads 4 --seed 123 --prefix $ODIR/MafftEinsi.trim

# parallelized ML search
mkdir -p $ODIR/ml_search
cd $ODIR/ml_search

#-- Vars in batch_raxML-ng.sh
# SEED=$1
# THREADS=$2
# TREES=$3 # type{number}
# MODEL=$4
# ALN=$5
# ODIR=$6 

# 100  ML searches starting with parsimony trees (10 runs x 10 per run)
seeds=($(seq 11 10 101))
for i in ${seeds[@]}; do
sbatch --job-name=MLpars$i --time=5:30:00 --cpus-per-task=4 \
$PARENT/code/batch_raxML-ng.sh $i 4 "pars{10}" $MODEL $ALN "$ODIR/ml_search"
done

# 100  ML searches starting with random trees (10 runs x 10 per run)
seeds=($(seq 13 10 103))
for i in ${seeds[@]}; do
sbatch --job-name=MLrand$i --time=5:30:00 --cpus-per-task=4 \
$PARENT/code/batch_raxML-ng.sh $i 4 "rand{10}" $MODEL $ALN "$ODIR/ml_search"
done


grep "Final LogLikelihood:" *.raxml.log | sort -k 3 
cat *.raxml.mlTrees > mltrees
raxml-ng --rfdist --tree mltrees --redo --prefix RF
# Loaded 200 trees with 661 taxa.
# Average absolute RF distance in this tree set: 289.487739
# Average relative RF distance in this tree set: 0.219975
# Number of unique topologies in this tree set: 200
#does not seem to converge



# Adding more trees to reach 500
# 150  ML searches starting with parsimony trees (10 runs x 10 per run)
seeds=($(seq 311 10 451))
for i in ${seeds[@]}; do
sbatch --job-name=MLpars$i --time=5:59:00 --cpus-per-task=4 \
$PARENT/code/batch_raxML-ng.sh $i 4 "pars{10}" $MODEL $ALN "$ODIR/ml_search"
done

# 150  ML searches starting with random trees (10 runs x 10 per run)
seeds=($(seq 313 10 453))
for i in ${seeds[@]}; do
sbatch --job-name=MLrand$i --time=5:59:00 --cpus-per-task=4 \
$PARENT/code/batch_raxML-ng.sh $i 4 "rand{10}" $MODEL $ALN "$ODIR/ml_search"
done

grep "Final LogLikelihood:" *.raxml.log | sort -k 3 
cat *.raxml.mlTrees > mltrees
raxml-ng --rfdist --tree mltrees --redo --prefix RF
# Loaded 500 trees with 661 taxa.
# Average absolute RF distance in this tree set: 290.622012
# Average relative RF distance in this tree set: 0.220837
# Number of unique topologies in this tree set: 500
#does not seem to converge

best_trees=($( grep "Final LogLikelihood:" *.raxml.log | sort -k 3 | head -n 5 | cut -c-8| tr -d . ))

best_trees=( "${best_trees[@]/%/.raxml.bestTree}" )
cat ${best_trees[@]} > mltrees_top
raxml-ng --rfdist --tree mltrees_top --redo --prefix RF_top
# Loaded 5 trees with 661 taxa.
# Average absolute RF distance in this tree set: 230.600000
# Average relative RF distance in this tree set: 0.175228
# Number of unique topologies in this tree set: 5
    
# using Best scoring ML tree
BEST=( $(grep "Final LogLikelihood:" *.raxml.log | sort -k 3 | head -n 1 | grep -o ".*.raxml")) 
cat "$BEST.bestTree" > ../bestScoreML.tree

# #reruns to complete time out
# i=391
# sbatch --job-name=MLpars$i --time=2:59:00 --cpus-per-task=4 \
# $PARENT/code/batch_raxML-ng.sh $i 4 "pars{10}" $MODEL $ALN "$ODIR/ml_search"
# 
# i=343
# sbatch --job-name=MLrand$i --time=0:59:00 --cpus-per-task=4 \
# $PARENT/code/batch_raxML-ng.sh $i 4 "rand{10}" $MODEL $ALN "$ODIR/ml_search"
# 
# cat *.raxml.mlTrees | wc -l
