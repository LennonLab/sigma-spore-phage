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
PARENT=~/GitHub/spore-phage-sigma/phylo
ODIR=${PARENT}/data/reduced_set_to_align
ALN=${ODIR}/sigmas_MafftEinsi.trim

#modeltest-ng selected model : LG+G4
cd $ODIR

raxml-ng --check --msa $ALN  \
--model LG+G4 --data-type AA \
--prefix check-msa 

# WARNING: Sequences YP_009818180.1-phage and YP_009798107.1-phage are exactly identical!
# WARNING: Duplicate sequences found: 1
# 
# NOTE: Reduced alignment (with duplicates and gap-only sites/taxa removed) 
# NOTE: was saved to: /geode2/home/u020/danschw/Carbonate/GitHub/spore-phage-sigma/phylo/data/reduced_set_to_align/check-msa.raxml.reduced.phy


ALN=${ODIR}/check-msa.raxml.reduced.phy

# For large alignments, we recommend using the --parse command after, or, instead of
# --check:
raxml-ng --parse --msa $ALN --data-type AA --model LG+G4 --prefix parse

ALN=${ODIR}/parse.raxml.rba

# * Estimated memory requirements                : 21 MB
# 
# * Recommended number of threads / MPI processes: 2


raxml-ng --all --msa $ALN --model LG+G4  --prefix tree-all --seed 123 --threads 2 \
--tree pars{10} --bs-trees autoMRE{200} --bs-metric fbp,tbe





# # get tree. Following preprint "Using RAxML-NG in Practice"
# # https://doi.org/10.20944/preprints201905.0056.v1
# 
# # coarse-grain parallelization can easily be emulated by executing multiple RAxML-NG instances, but with distinct random seeds
# 
# #make directory for slurm scripts
# sh_dir="$PARENT/code/tree_slurm"
#     if [ -f $sh_dir ]; then
#       rm -f $sh_dir
#     fi
# mkdir -p $sh_dir
# 
# #random seeds for 50 trees to find best ML
# seeds=($(seq 123 173))
# for seed in ${seeds[@]}; do
#   bash_out="$sh_dir/MLsearch_seed_$seed.sh"
# 
#   echo '#!/bin/bash' >> $bash_out
#   echo '#SBATCH --mail-user=danschw@iu.edu' >> $bash_out
#   echo '#SBATCH --nodes=1' >> $bash_out
#   echo '#SBATCH --ntasks-per-node=1' >> $bash_out
#   echo '#SBATCH --cpus-per-task=2' >> $bash_out
#   echo '#SBATCH --time=1:59:00' >> $bash_out
#   echo '#SBATCH --mem=1gb' >> $bash_out
#   echo '#SBATCH --mail-type=FAIL,BEGIN,END' >> $bash_out
#   echo "#SBATCH --job-name=${seed}" >> $bash_out
#   echo '' >> $bash_out
#   echo '##### load dependencies #####' >> $bash_out
#   echo 'module load raxmlng' >> $bash_out
#   echo '' >> $bash_out
#   echo '##### tree search #####' >> $bash_out
#   echo "raxml-ng --search --msa rbcl.raxml.rba"> $bash_out
# 
#    sbatch $bash_out
# 
# raxml-ng --msa $ALN --model LG+G4 --threads 1 --seed 123

