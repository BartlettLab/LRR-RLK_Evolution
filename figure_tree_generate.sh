	#!/bin/bash

#BSUB -W 48:00             	# How much time does your job need (HH:MM)
#BSUB -q long            	# Which queue to use {short, long, parallel, GPU, interactive}
#BSUB -R span[hosts=1]		# Keep on one CPU cluster
#BSUB -R rusage[mem=300]	# How much memory
#BSUB -n 17		            # Where X is in the set {1..X}
#BSUB -J figure_tre     	# Job Name
#BSUB -o out.%J           	# Append to output log file
#BSUB -e err.%J           	# Append to error log file

module load noisy/1.5.12 
module load iqtree/1.6.3
module load MAFFT/7.313
module load hmmer/3.1b2


#####this script will take the leaves of the extracted branch from a clade's tree and combine them and make a new tree with all######
#####Also performs a Pfam domain scan and outputs results as a table for reading in R########

#define outgroups
#	alignment_outgroup=AT1G75820.1 #CLV1, clade XI
#	alignment_outgroup=AT3G02130.1 #RKP2 clade XV
	alignment_outgroup=AT2G20850.1 #SRF1, clade V
	pelle_outgroup=DmPelle #Pelle kinase from Drosophila

#define run name
run_name_file=run_name*
run_name=$(echo $run_name_file | sed 's/.*=//') #only return that portion after the '='


echo "combine lists, add outgroup ($alignment_outgroup), consolidate"
	cat *taxa.txt > $run_name.input.temp #collect all genes from previous searches to be included in this tree
	echo $alignment_outgroup >> $run_name.input.temp #add the alignment outgroup to the gene list
	echo $pelle_outgroup >> $run_name.input.temp #add pelle outgroup to the gene list
	cat $run_name.input.temp | awk '!a[$0]++' > $run_name.input.taxa.txt #remove duplicates from gerne list
	rm $run_name.input.temp #remove the temp file
echo "done"

echo "constructing database"
	#construct primary transcript database
	cat $(ls ~/supertree/databases/primary_transcript_databases/*.oneline.fa) > combined_onelines.fa
	database=combined_onelines.fa
echo "done"


echo "grep all sequence variants for pfam domain table"
	#Strip transcript variants from IDs
	grep Zm*. $run_name.input.taxa.txt | awk -v FS="_" '{ print $1 }' >> $run_name.no_transcript_ID.txt
	grep Solyc*. $run_name.input.taxa.txt | awk -v FS="." '{ print $1 }' >> $run_name.no_transcript_ID.txt
	grep Pp3*. $run_name.input.taxa.txt | awk -v FS="." '{ print $1 }' >> $run_name.no_transcript_ID.txt
	grep Potri.*. $run_name.input.taxa.txt | awk -v FS="." -v OFS="." '{ $3 = "" ; print }' | sed 's/.$//' >> $run_name.no_transcript_ID.txt
	grep LOC_Os*. $run_name.input.taxa.txt | awk -v FS="." '{ print $1 }' >> $run_name.no_transcript_ID.txt
	grep Brad*. $run_name.input.taxa.txt | awk -v FS="." '{ print $1 }' >> $run_name.no_transcript_ID.txt
	grep AT.G*. $run_name.input.taxa.txt | awk -v FS="." '{ print $1 }' >> $run_name.no_transcript_ID.txt

	#these do not have additional transcript variants, so just push the IDs right over
	grep Smo_ $run_name.input.taxa.txt >> $run_name.no_transcript_ID.txt
	grep evm_27 $run_name.input.taxa.txt >> $run_name.no_transcript_ID.txt	
	
	#add outgroups with transcript variant stripped
	echo $(echo $alignment_outgroup |  sed 's/.\..*//')"0" >> $run_name.no_transcript_ID.txt	
	echo $pelle_outgroup >> $run_name.no_transcript_ID.txt	
	
	#Identify the Pre-constructed database with transcript variants delimed by @
	pfam_database=~/supertree/databases/pfam_all_transcripts/pfam_alltranscript_database.fa
	
	#Grab all the transcript variants of each taxon from the special pfam database
	fgrep -A 1 -f $run_name.no_transcript_ID.txt --no-group-separator $pfam_database | awk '{ print $1 }' > $run_name.pfam_transcript_variants.seqs.fa

echo "done"

echo "Building Pfam domain table" 	##in case need to regenerate pfam searchable database: $ hmmpress Pfam-A.hmm 
	hmmscan --noali -o pfamout.temp --cut_tc --tblout $run_name.pfamout.tsv ~/pfam/hmmfiles/Pfam-A.hmm $run_name.pfam_transcript_variants.seqs.fa
	rm pfamout.temp
echo "done"



echo "Collect primary transcript sequences for tree"
	fgrep -A 1 -f $run_name.no_transcript_ID.txt --no-group-separator $database > $run_name.input.seqs.fa
#	ID_count=$(cat $run_name.input.taxa.txt | wc -l)
#	seq_count=$(cat $run_name.input.seqs.fa | wc -l)
#	seq_count=$(($seq_count / 2)) #sequences take up two lines so divide the number in half
	#make sure all sequences are found. if not return an error and exit.
echo "done"



echo "aligning"
	mafft $run_name.input.seqs.fa > $run_name.input.align.fasta
echo "done"





echo "running noisy to clean up alignment"
	noisy -s --noconstant $run_name.input.align.fasta
	#remove noisy log files
	rm *typ.eps
	rm *sta.gr
	rm *idx.txt
echo "done"


echo "Running IQtree with auto find model of evo"
	iqtree -s *out.fas -bb 1000 -nt 8 -m JTT+F+R9
echo "done"	

echo "cleaning up"
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
	rm $database #remove the full oneline concatenated database
	rm $run_name.pfam_transcript_variants.seqs.fa
echo "exiting script"
exit