#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1:59:00
#SBATCH --mem=3gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=cd-85-95
#SBATCH --output=CD-HIT-all.%J.out
#SBATCH --error=CD-HIT-all.%J.err


module load cd-hit

PARENT="/geode2/home/u020/danschw/Carbonate/GitHub/sigma-spore-phage"
cd $PARENT

mkdir -p "vogdb/data/cd-clusts"

# MIUViG guide lines
# 95% average nucleotide identity 
# over 85% alignment fraction (relative to the shorter sequence)
cd-hit-est -aS .85 -c .95 -T 16 -M 0 -g 1 -sc 1 -d 0 \
-i vogdb/data/vog_all.fna \
-o vogdb/data/cd-clusts/vog_cd-clust_MIUViG 

 # -i   input filename in fasta format, required
 # -o   output filename, required
 # -aS  alignment coverage for the shorter sequence, default 0.0
 #      if set to 0.9, the alignment must covers 90% of the sequence
 # -c   sequence identity threshold, default 0.9
 #      this is the default cd-hit's "global sequence identity" calculated as:
 #      number of identical amino acids in alignment
 #      divided by the full length of the shorter sequence
 # -T   number of threads, default 1; with 0, all CPUs will be used
 # -M   memory limit (in MB) for the program, default 800; 0 for unlimitted;
 # -g	1 or 0, default 0
 #      by cd-hit's default algorithm, a sequence is clustered to the first 
 #      cluster that meet the threshold (fast cluster). If set to 1, the program
 #      will cluster it into the most similar cluster that meet the threshold
 #      (accurate but slow mode)
 #      but either 1 or 0 won't change the representatives of final clusters     
 # -sc	sort clusters by size (number of sequences), default 0, output clusters by decreasing length
 # 	if set to 1, output clusters by decreasing size
 # -d	length of description in .clstr file, default 20
 # 	if set to 0, it takes the fasta defline and stops at first space
