#!/bin/bash

#BSUB -W 24:00             	# How much time does your job need (HH:MM)
#BSUB -R rusage[mem=15000]	# How much memory
#BSUB -n 1                	# Where X is in the set {1..X}
#BSUB -J search     	# Job Name
#BSUB -o out.%J           	# Append to output log file
#BSUB -e err.%J           	# Append to error log file
#BSUB -q long            	# Which queue to use {short, long, parallel, GPU, interactive}


##Usage
#Dependancies:
#ensure .fasta subalignments are in working directory
#This needs eval_helper file, also in working directory, with .fasta files listed and their evalues specified
#Modify the run name file, should be something like run_name=X.round3
#modify the outgroup to use below, should be something like alignment_outgroup=AT1G75820.1 #CLV1

############################################


blast_identity_percent_cutoff=60
#Set blasting cutoff, this is %identity to target seqs, use 60% for subdomains, #50% finds all known clade XI core taxa,  60% misses about 4% of those but only returns an enriched list 1/3 the size of 50%


module load MAFFT/7.313
module load hmmer/3.1b2
module load blast/2.2.22

known_taxa=*input.taxa.txt #read the gene names of known members
echo "$(wc -l < $known_taxa) known genes are used as input"

#If CLV1 is in the search terms, use RPK2, otherwise use CLV1 as outgroup
if grep -q AT1G75820 $known_taxa
	then
	alignment_outgroup=AT3G02130.1 #RKP2
	else
	alignment_outgroup=AT1G75820.1 #CLV1
	fi
#alignment_outgroup=AT2G25880 #kinase-like outgroup


run_name_file=run_name*
run_name=$(echo $run_name_file | sed 's/.*=//') #only return that portion after the '='

#construct sequence database - inputs needs to be oneline
echo "combining oneline databses and formatting for blasting..."
cat $(ls ~/supertree/databases/primary_transcript_databases/*.oneline.fa) > combined_onelines.fa
database=combined_onelines.fa
formatdb -i $database  

#Report databse entries quantity
let lines_in_combined_database=$(wc -l < $database)
let lines_in_combined_database=lines_in_combined_database/2
echo "$lines_in_combined_database genes being searched in $(ls ~/supertree/databases/primary_transcript_databases/*.oneline.fa | wc -l) genomes"

#set up BLast directory
blast_dir="Blast_results"
mkdir $blast_dir

 #set up HMM directory
HMMdir="HMM_results"
mkdir $HMMdir

echo "HMM search results will be thresholded at $alignment_outgroup, then thresholded to $blast_identity_percent_cutoff% identity matches to one of the input sequences"

######################### Main search algorithm.#############################

echo "Begin search loops..."
for input_alignment in *.fasta
do
	
	#save the run name as the file name without extension
	subalignment_name=${input_alignment%.fasta}	
#------HMM SECTION-------
	#make run name directory in the outgroup directory
	mkdir $HMMdir/$subalignment_name


	##build and run the HMM search.
	
	hmmbuild -o hmmout.txt $HMMdir/$subalignment_name/model.hmm $input_alignment
	hmmsearch -o hmmout.txt --noali --tblout $HMMdir/$subalignment_name/search.output.txt $HMMdir/$subalignment_name/model.hmm $database
	rm hmmout.txt #report file generated, couldn't figure out how to prevent it


	#using HMMsearch output file search.output.txt, make clean search output for processing called search.output.clean
	sed -n -e "4,$(($(wc -l < $HMMdir/$subalignment_name/search.output.txt) - 10))p" $HMMdir/$subalignment_name/search.output.txt > $HMMdir/$subalignment_name/search.output.clean

	#find line with alignment_outgroup | grep to isolate just the line #
		#if not found, use entire list
	if grep -q $alignment_outgroup $HMMdir/$subalignment_name/search.output.clean; then
			alignment_outgroup_position=$(grep -n $alignment_outgroup $HMMdir/$subalignment_name/search.output.clean | grep -Eo '^[^:]+')
			#echo "$alignment_outgroup found at line $alignment_outgroup_position"
		else 
			alignment_outgroup_position=$(cat $HMMdir/$subalignment_name/search.output.clean | wc -l)
			#echo "$alignment_outgroup not found, setting cutoff to end of results at line $alignment_outgroup_position"
	fi
			
		
		#abidge results to outgroup position and add a space after
	head -$alignment_outgroup_position $HMMdir/$subalignment_name/search.output.clean | awk '{print $1}' > $HMMdir/$subalignment_name/abridged_hits.txt


#	get full sequences of HMM hits 

			grep -A 1 -f $HMMdir/$subalignment_name/abridged_hits.txt --no-group-separator $database | awk '{print $1}' >> $HMMdir/$subalignment_name/HMM_abridged_seqs.fa
	
#	once the hmm enriched seqs are collected:	
####	blast the  seqs to threshold again - blast all HMM hits against input
#	format HMM results as a database
	formatdb -i $HMMdir/$subalignment_name/HMM_abridged_seqs.fa
#	Find the percent ID for each seach term to the HMM results databse
	blastall -m 8 -p blastp -d $HMMdir/$subalignment_name/HMM_abridged_seqs.fa -i $input_alignment -e 1e-10 | awk '{print $2,$3}' > $HMMdir/$subalignment_name/HMM_hits_blast_results.txt
	
#	collect only hits that pass the blast %identity threshold cutoff
	awk " \$2 > $blast_identity_percent_cutoff " $HMMdir/$subalignment_name/HMM_hits_blast_results.txt | awk '{print $1}' | awk '!a[$0]++' > $HMMdir/$subalignment_name/HMM_hits_blast_ID_percent_thresholded.txt
	
	#add these thresholded hits to cumulative list
	cat $HMMdir/$subalignment_name/HMM_hits_blast_ID_percent_thresholded.txt >> $subalignment_name.cumulative_enriched_taxa.txt
	
	###what is this?
	hmm_count=$(wc -l < $HMMdir/$subalignment_name/HMM_hits_blast_ID_percent_thresholded.txt )
	echo "HMM $subalignment_name results: $(wc -l < $HMMdir/$subalignment_name/abridged_hits.txt) hits based on outgroup, re-thresholded by blast ID% to $hmm_count hits"
	echo "Out of $(wc -l < $known_taxa) known taxa, $subalignment_name HMM found $( grep -f $known_taxa $HMMdir/$subalignment_name/HMM_hits_blast_ID_percent_thresholded.txt |  awk '{print $1}' | awk '!a[$0]++' | wc -l)"

#	clean up formatdb files
	rm $HMMdir/$subalignment_name/HMM_abridged_seqs.fa.phr
	rm $HMMdir/$subalignment_name/HMM_abridged_seqs.fa.pin
	rm $HMMdir/$subalignment_name/HMM_abridged_seqs.fa.psq
	rm $HMMdir/$subalignment_name/search.output.txt
	





	#-------BLAST SECTION-----------
	
	#find the input alignment file name in eval_helper
		blast_evalue_cutoff=$(grep $input_alignment *eval_helper.txt | awk '{print $2}')

		
	#Blast each sequence in fasta, print relevant results -m 8 says only print table, awk separates just the columns to see, full column list can be found at http://www.pangloss.com/wiki/Blast
	# method: "grep e- | awk '{print $1, substr($2,4); }'" processes 7e-19 to 19
	blastall -m 8 -p blastp -d $database -i $input_alignment | awk '{print $2,$11}' | grep e- | awk '{print $1, substr($2,4); }' > $blast_dir/$subalignment_name.raw_blast_results
			#awk '{print $2,$3} will return search hit ID and % identity
			#awk '{print $2,$11} will return search hit ID and e-value
	
	#clean up blast results: 
		#remove hits with ID% less than blast cutoff 
		#print only the gene name
		#remove duplicates
		#add results to  to cumuluative list
		awk " \$2 > $blast_evalue_cutoff " $blast_dir/$subalignment_name.raw_blast_results | awk '{print $1}' | awk '!a[$0]++' >> $subalignment_name.cumulative_enriched_taxa.txt
	
	echo "Blast $subalignment_name results after e-value thresholding at 1e-$blast_evalue_cutoff: $(awk " \$2 > $blast_evalue_cutoff " $blast_dir/$subalignment_name.raw_blast_results | awk '{print $1}' | awk '!a[$0]++'| wc -l) hits"
	echo "Out of $(wc -l < $known_taxa) known taxa, $subalignment_name BLAST found $(awk " \$2 > $blast_evalue_cutoff " $blast_dir/$subalignment_name.raw_blast_results | awk '{print $1}' | awk '!a[$0]++' | grep -f $known_taxa |  awk '{print $1}' | awk '!a[$0]++' | wc -l)"
	
#-------------END BLAST SECTION-------------


	# consolidate the results from each search and pull their full sequences from the database

	#sort and remove duplicates using magic awk command | use sed to add space at end of each line so Amborella doesn't pull a million  matches
	cat $subalignment_name.cumulative_enriched_taxa.txt | awk '!a[$0]++' > $subalignment_name.enriched_taxa.txt
	#get rid of cumulative file
	rm $subalignment_name.cumulative_enriched_taxa.txt

	echo "$subalignment_name unique enriched taxa found: $(wc -l < $subalignment_name.enriched_taxa.txt)"
	echo "Out of $(wc -l < $known_taxa) known taxa, $subalignment_name searches found $(grep -f $known_taxa $subalignment_name.enriched_taxa.txt |  awk '{print $1}' | awk '!a[$0]++' | wc -l)"
	echo "-----------end $subalignment_name subalignment searches-----------"

done
###################END blast and HMM search loop here##########################

#clean up files for blast formatted database
rm $database.p*
rm formatdb.log

#report on inputs found in all searches
#consolidate hits, remove duplicates
cat *enriched_taxa.txt | awk '!a[$0]++' > $run_name.all_runs_enriched_taxa.txt
echo "Unique hits from all searches: $(wc -l < $run_name.all_runs_enriched_taxa.txt)"
echo "out of $(wc -l < $known_taxa) known taxa, combined searches found $(cat $run_name.all_runs_enriched_taxa.txt | grep -f $known_taxa |  awk '{print $1}' | awk '!a[$0]++' | wc -l)"


#this section checks for genes to keep an eye on
	
	if grep -q AT1G75820 $run_name.all_runs_enriched_taxa.txt; then
	echo "CLV1 is found"
	else
	echo "CLV1 is not found"
	fi

if grep -q AT1G65380 $run_name.all_runs_enriched_taxa.txt; then
	echo "CLV2 is found"
	else
	echo "CLV2 is not found"
	fi

if grep -q Zm00001d040130 $run_name.all_runs_enriched_taxa.txt; then
	echo "FEA3 is found"
	else
	echo "FEA3 is not found"
	fi

if grep -q AT5G13290 $run_name.all_runs_enriched_taxa.txt; then
	echo "CRN is found"
	else
	echo "CRN is not found"
	fi

if grep -q AT5G48740 $run_name.all_runs_enriched_taxa.txt; then
	echo "clade I RLK in a rlp branch is found"
	else
	echo "clade I RLK in a rlp branch is not found"
	fi

if grep -q AT2G31880 $run_name.all_runs_enriched_taxa.txt; then
	echo "evershed is found"
	else
	echo "evershed is not found"
	fi

#######Grep and align section_##############
### creates a taxa list, sequence list, and alignment for LRR finds, Kinase finds, and total finds###

#pull all the unique taxa found from searches from the database


#consolidate kinase results, pull seqs, align
echo "preparing alignment from kinase search results..."
cat $(ls kin*enriched_taxa.txt) | awk '!a[$0]++' > $run_name.kinase_enriched_taxa.txt
rm $(ls kin*enriched_taxa.txt)
grep -w -A 1 -f $run_name.kinase_enriched_taxa.txt --no-group-separator $database | awk '{print $1}' > $run_name.kinase.enriched_seqs.fa #awk command will clean up extra junk after taxa name
echo "manually adding $alignment_outgroup to the kinase alignment"
grep -A 1 $alignment_outgroup --no-group-separator $database >> $run_name.kinase.enriched_seqs.fa
	#make sure all sequences are found. if not return an error and exit.
	ID_count=$(cat $run_name.kinase_enriched_taxa.txt | wc -l)
	seq_count=$(cat $run_name.kinase.enriched_seqs.fa | wc -l)
	seq_count=$(($seq_count / 2)) #sequences take up two lines so divide the number in half
	seq_count=$(($seq_count -1)) #remove the outgroup that was added after IDcount taken
	if [ $ID_count != $seq_count ]
	then
		echo "error finding all sequences, exiting."
		exit
	else
		echo "$ID_count sequences collected"
	fi		
mafft $run_name.kinase.enriched_seqs.fa > $run_name.kinase.enriched_align.fasta
echo "finished with kinase alignment"

#consolidate lrr results, pull seqs, align
echo "preparing alignment from lrr search results..."
cat $(ls lrr*enriched_taxa.txt) | awk '!a[$0]++' > $run_name.lrr_enriched_taxa.txt
rm $(ls lrr*enriched_taxa.txt)
grep -w -A 1 -f $run_name.lrr_enriched_taxa.txt --no-group-separator $database | awk '{print $1}' > $run_name.lrr.enriched_seqs.fa #awk command will clean up extra junk after taxa name
echo "manually adding $alignment_outgroup to the lrr alignment"
grep -A 1 $alignment_outgroup --no-group-separator $database >> $run_name.lrr.enriched_seqs.fa
	#make sure all sequences are found. if not return an error and exit.
	ID_count=$(cat $run_name.lrr_enriched_taxa.txt | wc -l)
	seq_count=$(cat $run_name.lrr.enriched_seqs.fa | wc -l)
	seq_count=$(($seq_count / 2)) #sequences take up two lines so divide the number in half
	seq_count=$(($seq_count -1)) #remove the outgroup that was added after IDcount taken
	if [ $ID_count != $seq_count ]
	then
		echo "error finding all sequences, exiting."
		exit
	else
		echo "$ID_count sequences collected"
	fi		

mafft --bl 45 $run_name.lrr.enriched_seqs.fa > $run_name.lrr.enriched_align.fasta
echo "finished with lrr alignment"

# consolidate all results, pull seqs, align
# echo "preparing alignment from all search results..."
# cat  all_lrr_enriched_taxa.txt all_kinase_enriched_taxa.txt  | awk '!a[$0]++' > all_enriched_taxa.txt
# grep -A 1 -f all_enriched_taxa.txt --no-group-separator $database | awk '{print $1}' > all_enriched_seqs.fa #awk command will clean up extra junk after taxa name
# add CLV1
# grep -A 1 AT1G75820 --no-group-separator $database >> all_enriched_seqs.fa
# mafft --bl 45 all_enriched_seqs.fa> all_enriched_align.fasta
# echo "finished with combined search hits alignments"

echo "finished with all alignments"
rm $database
rm *.log






exit





