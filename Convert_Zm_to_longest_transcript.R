setwd("/Users/jarrettman/Documents/Bartlett_Lab/Phylogenetics/LRR_RLKs/supertree_scratch")
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
library(Biostrings)

## read your fasta in as Biostrings object
fasta.s <- readDNAStringSet("Zea_mays.AGPv4.pep.all.fa.oneline.fa")

## get the read names (in your case it has the isoform info)
names.fasta <- names(fasta.s)

## extract only the relevant gene and isoform id (split name by the period symbol)
gene.iso <- sapply(names.fasta,function(j) cbind(unlist(strsplit(j,'\\_'))[1:2]))

## convert to good data.frame = transpose result from previous step and add relevant column names
gene.iso.df <- data.frame(t(gene.iso))
colnames(gene.iso.df) <- c('gene','isoform')

## and length of isoforms
gene.iso.df$width <- width(fasta.s)

## split data.frame into list with entry for each gene
gene.iso.df.split <- split(gene.iso.df,gene.iso.df$gene)

## pull out the longest isoform ID for each gene (in case of a tie just take the first one)
#this will need to be cleaned up in text editor, but is working correctly
best.id <- sapply(gene.iso.df.split,function(x) row.names(x)[order(x$width, decreasing = TRUE)[1]])
write.csv(best.id, file = "bestID.csv")

#AFTER CLEANING, you can use this to pull the correct transcript variants from the all transcript variants file
grep -A 1 -f Zm_longest_transcript_list --no-group-separator Zea_mays.AGPv4.pep.all.fa.oneline.fa | awk '{print$1}' > Zea_mays.AGPv4.pep.longest.oneline.fa




#below are functions from internets to actually print the gene sequences, for a full database this is unadvisable in R and on a personal machine

## subset your original reads with the subset
fasta.s.best <- fasta.s[best.id]

## export new fastafile containing longest isoform per gene
writeXStringSet(fasta.s, filepath='sample_best_isoform.fasta')
