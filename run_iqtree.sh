#!/bin/bash

#BSUB -W 40:00             	# How much time does your job need (HH:MM)
#BSUB -q long            	# Which queue to use {short, long, parallel, GPU, interactive}
#BSUB -R span[hosts=1]		# Keep on one CPU cluster
#BSUB -R rusage[mem=250]	# How much memory
#BSUB -n 9            # Where X is in the set {1..X}
#BSUB -J IQtree     	# Job Name
#BSUB -o out.%J           	# Append to output log file
#BSUB -e err.%J           	# Append to error log file

#################
#Note optimum perfromace found with 8 threads ( -nt 8 ) and 16 CPU nodes #BSUB -n 16 
#may need to tone that down for smaller trees - IQtree does't like lots of threads on small trees
##############

module load noisy/1.5.12 
module load iqtree/1.6.3
module load MAFFT/7.313

run_name_file=run_name*
run_name=$(echo $run_name_file | sed 's/.*=//') #only return that portion after the '=', save that as run name variable

	cat *geneIDs.txt | awk '!a[$0]++' > $run_name.gene_names.input.txt #consolidate gene list

#	cat $(ls ~/supertree/databases/primary_transcript_databases/*.oneline.fa) > combined_onelines.fa #build database
	database=combined_onelines.fa #assign name
	
#	fgrep -w -A 1 -f  $run_name.gene_names.input.txt --no-group-separator $database > $run_name.input.seqs.fa #collect all sequences
	
	mafft $run_name.input.seqs.fa > $run_name.input.align.fasta #align all sequences




noisy -s --noconstant $run_name.input.align.fasta #clean with noisy

#remove noisy log files
rm *typ.eps
rm *sta.gr
rm *idx.txt

echo "Running IQtree"
# fast option, initial tree based on parsimony and BIONJ, don't check for threads, evo model, stop after 100 iterations, no bootstrap
#	iqtree -s *out.fas -fast -m JTT -nt 8 

#medium version
	iqtree -s *out.fas -bb 1000 -nt AUTO -m JTT -bcor 0.97


#Long and careful version: If want to find model of evolution, no shortcuts
#	iqtree -s *out.fas -bb 3000 
	
#careful and use partition
#	iqtree -s *out.fas -bb 1000 -nt 8 -spp *.nex
	
	
#remove iqtree log files
rm *splits.nex
rm *ckp.gz
rm *bionj
rm *treefile
rm *mldist
rm *.log
rm *.iqtree
rm *model.gz
rm *uniqueseq.phy
