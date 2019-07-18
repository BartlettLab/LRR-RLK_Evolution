#!/bin/bash

#BSUB -W  4:00             	# How much time does your job need (HH:MM)
#BSUB -R rusage[mem=1000]	# How much memory
#BSUB -n 4                	# Where X is in the set {1..X}
#BSUB -J eval_find   	# Job Name
#BSUB -o out.%J           	# Append to output log file
#BSUB -e err.%J           	# Append to error log file
#BSUB -q short            	# Which queue to use {short, long, parallel, GPU, interactive}
#BSUB -R span[hosts=1]

###Usage
#move this script to your working directory
#ensure that the only .fasta files in the directory are sub-alignments of your targets
#ensure that a list of your targets is also present in the directory, and ends in .input.taxa.txt

module load MAFFT/7.313
module load hmmer/3.1b2
module load blast/2.2.22 

run_name_file=run_name*
run_name=$(echo $run_name_file | sed 's/.*=//') #only return that portion after the '='

#set up eval_helper document
echo "Usage of this document:" >> $run_name.eval_helper.txt
echo "Run blast_eval_finder.sh to generate eval_finder_results.csv showing eval curve for data" >> $run_name.eval_helper.txt
echo "Pick an eval on the curve that is near the point between steep drop-off and leveling off. Hint - import to excela dn select the entire column and click line graph" >> $run_name.eval_helper.txt
echo "For each input sub-alignment, write the full filename here on one line, then a space, then the e-value chosen above, as a two digit number, ie, an e-value of 1e-17 would be reported here as 17" >> $run_name.eval_helper.txt
echo "This file is grepped by the search script, so make sure each sub-alignment has its own line and the file name is not present twice" >> $run_name.eval_helper.txt
echo " " >> $run_name.eval_helper.txt
echo $run_name >> $run_name.eval_helper.txt

#find input file
input_taxa=$(ls *input.taxa.txt)

csv_output_file=$run_name.eval_finder_results.csv


#create a single searchable database for all species genomes
cat $(ls ~/supertree/databases/primary_transcript_databases/*.oneline.fa) > combined_onelines.fa
database=combined_onelines.fa
formatdb -i $database #required for blast function

#count the number of taxa used as input
input_taxa_total=$(wc -l < $input_taxa)
echo "input taxa found: $input_taxa_total"

#set up results file
 echo "Input sequence,E-value cutoff,Number of hits retained,Inputs retained,Non-inputs retained" > $csv_output_file


#Blast each subalignment against the consolidated database
for search_terms_subalignment in *.fasta #use each sub alignment in the current directory
do
	#add subalignment name to eval helper list
	echo "$search_terms_subalignment " >> $run_name.eval_helper.txt
	
	echo "Blasting entire database with $search_terms_subalignment"
	#Blast all sequences in subalignment to database, return all hits as table|print the taxa and eval,| only with real eval, | remove e-	> save		
	blastall -m 8 -p blastp -d $database -i $search_terms_subalignment -a 8 | awk '{print $2,$11}' | grep e- | awk '{print $1, substr($2,4); }'  > $search_terms_subalignment.unfiltered_hits.txt

	 #Loop to catalog the blast hits collection by e-val threshold. Good idea to set near # residues in longest subalignment
	 eval_test_number=1
	 while [ $eval_test_number -lt 150 ]
	 do
		 #reduce file to only hits that clear threshold
		 hits_number=$(awk " \$2 > $eval_test_number " $search_terms_subalignment.unfiltered_hits.txt | awk '{print $1}' | awk '!a[$0]++' | wc -l) #report total number of blast hits
		 input_taxa_found=$(awk " \$2 > $eval_test_number " $search_terms_subalignment.unfiltered_hits.txt | awk '{print $1}' | awk '!a[$0]++' | grep -f $input_taxa | wc -l ) #number of input taxa found in the blast hits
		 delta=$((hits_number-input_taxa_found)) #How many blast hits are not part of the input group
			
		#report finds in standard out	
		 echo "thresholded at 1e-$eval_test_number, total hits = $hits_number, and $input_taxa_found hits from the known $input_taxa_total, difference of $delta"
		 #report finds by adding to file eval_finder_results.csv
		 echo "$search_terms_subalignment,$eval_test_number,$hits_number,$input_taxa_found,$delta" >> $csv_output_file
		 eval_test_number=$[$eval_test_number+1] #move to next eval number
	 done
done

#clean up blast database files
rm *.phr
rm *.pin
rm *.psq
rm formatdb.log
rm error.log
#rm *unfiltered_hits.txt
rm $database

exit

