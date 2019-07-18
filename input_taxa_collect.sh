#!/bin/bash

#BSUB -W  4:00             	# How much time does your job need (HH:MM)
#BSUB -R rusage[mem=1000]	# How much memory
#BSUB -n 1                	# Where X is in the set {1..X}
#BSUB -J inp_coll     	# Job Name
#BSUB -o out.%J           	# Append to output log file
#BSUB -e err.%J           	# Append to error log file
#BSUB -q short            	# Which queue to use {short, long, parallel, GPU, interactive}

##USAGE#######################
#move this script to the directory that your inputs file is contained in
#modify the input_file line to reflect your inputs file name
#it should look like X.round3.lrr.input
#This  will convert the input list to an exact match of database ID > X.round3.lrr.input.taxa.txt
#Then it will pull their sequences from the database > X.round3.lrr.input.seqs.fa
# and align > X.round3.lrr.input.align.fasta


module load MAFFT/7.313 

run_name_file=run_name*
run_name=$(echo $run_name_file | sed 's/.*=//') #only return that portion after the '='


#construct database
cat $(ls ~/supertree/databases/primary_transcript_databases/*.oneline.fa) > combined_onelines.fa
database=combined_onelines.fa

#inputs files Should be <group>.<round>.input
for input_file in *.input
	do
	#Modify the input taxa to reflect exactly the taxa ID found in the database, and then remove duplicates
	grep -f $input_file $database | awk '{print$1}' | cut -c2- | awk '!a[$0]++' > $input_file.taxa.txt
	grep -A 1 -f $input_file.taxa.txt --no-group-separator $database > $input_file.seqs.fa

	echo "Aligning..."
	mafft $run_name.seqs.fa >  $input_file.align.fasta
	echo "Finished with  $input_file"

done

rm $database


