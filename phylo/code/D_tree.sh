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
ODIR=${PARENT}/data/align-trim-tree
ALN=$PARENT/data/align-trim-tree/sigmas_MafftEinsi.trim

#modeltest-ng selected model : LG+G4
cd $ODIR

raxml-ng --check --msa $ALN  \
--model LG+G4 --data-type AA 
#--prefix check-msa 

# WARNING: Duplicate sequences found: 71
# 
# NOTE: Reduced alignment (with duplicates and gap-only sites/taxa removed) 
# NOTE: was saved to: /geode2/home/u020/danschw/Carbonate/GitHub/spore-phage-sigma/phylo/data/align-trim-tree/sigmas_MafftEinsi.trim.raxml.reduced.phy
ALN=$PARENT/data/align-trim-tree/sigmas_MafftEinsi.trim.raxml.reduced.phy

# For large alignments, we recommend using the --parse command after, or, instead of
# --check:
raxml-ng --parse --msa $ALN --data-type AA --model LG+G4

ALN=$PARENT/data/align-trim-tree/sigmas_MafftEinsi.trim.raxml.reduced.phy.raxml.rba

# * Estimated memory requirements                : 24 MB
# 
# * Recommended number of threads / MPI processes: 1

# get tree
raxml-ng --msa $ALN --model LG+G4 --threads 1 --seed 123

