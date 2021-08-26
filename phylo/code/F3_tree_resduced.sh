#--------------------------
# phylogeny with  RAxML-NG
#--------------------------
# # Following preprint "Using RAxML-NG in Practice"
# # https://doi.org/10.20944/preprints201905.0056.v1

#Interactive job on Carbonate
srun -p interactive -N 1 --ntasks-per-node=1 --cpus-per-task=8 --time=07:59:00 --pty bash

#### load dependencies ####

module load raxmlng
  # Flex 2.6.4 loaded.
  # CMake version 3.8.0 loaded.
  # raxmlng version 0.9.0-pthreads loaded.

##### Define paths #####
PARENT=~/GitHub/spore-phage-sigma/phylo
ODIR=${PARENT}/data/reduced_set_to_align/check_msa
mkdir -p $ODIR
ALN=$PARENT/data/reduced_set_to_align/sigmas_MafftEinsi.trim


#modeltest-ng selected model : LG+G4
cd $ODIR

raxml-ng --check --msa $ALN  \
--model LG+G4 --data-type AA \
--prefix check-msa 

# WARNING: Duplicate sequences found: 19
# 
# NOTE: Reduced alignment (with duplicates and gap-only sites/taxa removed) 
# NOTE: was saved to: /geode2/home/u020/danschw/Carbonate/GitHub/spore-phage-sigma/phylo/data/reduced_set_to_align/check_msa/check-msa.raxml.reduced.phy

# checked that sequenced removed as duplicates do not include sequences directly discussed in MS
# reduced_duplicated-seqs.R

ALN=${ODIR}/check-msa.raxml.reduced.phy

# For large alignments, we recommend using the --parse command after, or, instead of
# --check:
raxml-ng --parse --msa $ALN --data-type AA --model LG+G4 --prefix parse

ALN=${ODIR}/parse.raxml.rba

# * Estimated memory requirements                : 27 MB
# 
# * Recommended number of threads / MPI processes: 2

ODIR=${PARENT}/data/reduced_set_to_align/tree
mkdir -p $ODIR
cd $ODIR


#### ML tree ####
SEED=123
# 25 ML searches starting with parsimony trees
sbatch --job-name=MLpars$SEED $PARENT/code/batch_ML_parsimony.sh $SEED 25
# 25 ML searches starting with random trees
sbatch --job-name=MLrand$SEED $PARENT/code/batch_ML_random.sh $SEED 25

# Check if an optimum has been reached
grep "Final LogLikelihood:" *$SEED.raxml.log
# pars-123.raxml.log:Final LogLikelihood: -28077.266208
# rand-123.raxml.log:Final LogLikelihood: -28075.232741
# random is doing a little better
cat *.raxml.mlTrees > mltrees
raxml-ng --rfdist --tree mltrees --prefix RF

# Average absolute RF distance in this tree set: 76.006531
# Average relative RF distance in this tree set: 0.212309
# Number of unique topologies in this tree set: 50

# increase ML search
mkdir -p $ODIR/ml_search
cd $ODIR/ml_search

# 500  ML searches starting with parsimony trees (50 runs x 10 per run)
seeds=($(seq 11 10 501))
for i in ${seeds[@]}; do
sbatch --job-name=MLpars$i --time=0:45:00 \
$PARENT/code/batch_ML_parsimony.sh $i 10
done

# 500  ML searches starting with random trees (50 runs x 10 per run)
seeds=($(seq 13 10 503))
for i in ${seeds[@]}; do
sbatch --job-name=MLrand$i --time=0:45:00 \
$PARENT/code/batch_ML_random.sh $i 10
done

cd ..
grep "Final LogLikelihood:" *.raxml.log | sort -k 3 
cat *.raxml.mlTrees > mltrees2
raxml-ng --rfdist --tree mltrees2 --redo --prefix RF2
#does not seem to converge
# using Best scoring ML tree
grep "Final LogLikelihood:" *.raxml.log | sort -k 3 | head -n 3

BestTree=$ODIR/rand-313.raxml.bestTree

#### Branch support ####
# 1000 bootstrap trees (40 runs x 25 per run)
seeds=($(seq 1 40))
for i in ${seeds[@]}; do
sbatch --job-name=bs$i $PARENT/code/batch_bootstrap.sh $i
done

#check for convergence
cat bootstraps/*.bootstraps > allbootstraps
raxml-ng --bsconverge --bs-trees allbootstraps --prefix 1000BS --seed 2 \
--prefix check_bs_conv --seed 123 --threads 2 --bs-cutoff 0.03

# Bootstopping test converged after 950 trees !!!

# Consensus tree building
# raxml-ng --consense MRE --tree allbootstraps --prefix consMRE
# # this yielded a tree with multifurications :(

#### compute support for best scoring ML tree ####
raxml-ng --support \
--tree $BestTree \
--bs-trees allbootstraps \
--prefix ML_TBE_tree --threads 2 --bs-metric tbe,fbp



# old code

# raxml-ng --all --msa $ALN --model LG+G4  --prefix tree-all --seed 123 --threads 2 \
# --tree pars{10} --bs-trees autoMRE{200} --bs-metric fbp,tbe
# 
# # check for covergrence
# raxml-ng --bsconverge --bs-trees tree-all.raxml.bootstraps \
# --prefix check_bs_conv --seed 123 --threads 2 --bs-cutoff 0.01
# 
# #no convergence run 200 more
# raxml-ng --bootstrap --msa $ALN --model LG+G4  --prefix bs-to400 --seed 321 --threads 2 \
# --tree pars{10} --bs-trees autoMRE{200} --bs-metric fbp,tbe
# 
# #check for convergence
# cat *.bootstraps > allbootstraps
# raxml-ng --bsconverge --bs-trees allbootstraps --prefix T12 --seed 2 \
# --prefix check_bs_conv --seed 123 --threads 2 --bs-cutoff 0.01