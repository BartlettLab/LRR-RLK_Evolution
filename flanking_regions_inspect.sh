#!/bin/bash

#BSUB -W 4:00             	# How much time does your job need (HH:MM)
#BSUB -q short            	# Which queue to use {short, long, parallel, GPU, interactive}
#BSUB -R span[hosts=1]		# Keep on one CPU cluster
#BSUB -R rusage[mem=1000]	# How much memory
#BSUB -n 16                  # Where X is in the set {1..X}
#BSUB -J flanks          	# Job Name
#BSUB -o out.%J           	# Append to output log file
#BSUB -e err.%J           	# Append to error log file


module load noisy/1.5.12 
module load iqtree/1.6.3
module load MAFFT/7.313
module load bedtools/2.26.0	
module load hmmer/3.1b2 
module load EMBOSS/6.6.0
module load R/3.6.0
module load gcc/8.1.0



#generate the fasta nucleotide sequence of intergenic spaces from the bed file    
    #bedtools getfasta -fi Athaliana_167_TAIR9.fa -bed target_intergenic_space.bed -fo target_intergenic_space.fa -name #Amborella version
    #bedtools getfasta -fi Athaliana_167_TAIR9.fa -bed target_intergenic_space.bed -fo target_intergenic_space.fa -name #arabidopsis version
    #bedtools getfasta -fi Osativa_204_v7.0.fa -bed target_intergenic_space.bed -fo target_intergenic_space.fa -name #rice version
    bedtools getfasta -fi ../*.fa -bed target_intergenic_space.bed -fo target_intergenic_space.fa -name #generalized version, get fasta from parent folder

#remove coordinates in name, otherwise translation will only keep these and forget geneID
    sed -i 's/::.*//g' target_intergenic_space.fa


#translate in all 6 frames
    transeq target_intergenic_space.fa target_intergenic_space_translations.fa -auto -stdout -frame 6 -table 0 #all 6 frames
    #transeq target_intergenic_space.fa target_intergenic_space_translations.fa -auto -stdout -frame=F -table 0 #if you want only Forward or reverse frames

#scan all regions with pfam with VERY low cutoff threshold (bit score -100 or greater)
	hmmscan -T -100 --noali -o pfamout.temp --tblout pfamout.tsv ~/pfam/hmmfiles/Pfam-A.hmm target_intergenic_space_translations.fa
    rm pfamout.temp #log of every sequence hits, any found will be in table

#generate a report of the findings
    hits_number=$(cat pfamout.tsv | sed '1,3d' | tac | sed '1,10d' | wc -l)
    echo "Total number of hits found: $hits_number" >> log.txt
    echo "LRR and RLK hits shown below" >> log.txt
    head -3 pfamout.tsv >> log.txt
    grep "LRR" pfamout.tsv >> log.txt
    grep "Pkinase" pfamout.tsv >> log.txt
        lrr_score=$(grep "LRR" pfamout.tsv |  awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }')
        rlk_score=$(grep "Pkinase" pfamout.tsv |  awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }')
    echo "Average bit score of LRR domains found: $lrr_score (trusted cutoff < 27.0)" >> log.txt
    echo "Average bit score of RLK domains found: $rlk_score (trusted cutoff < 23.1)">> log.txt
