#### this script uses a genome .gff3 file and takes as argument a list of gene IDs, and produces a set of coordinates to scan for missing LRR and RLK domains. 
#install.packages("tidyverse")
#library(tidyverse)

#update this section for each run
setwd("~/supertree/Annotation_checks/") #go to the directory for this run.
domain_found_version <- "LRR_only"
domain_found_version <- "RLK_only"





####this section for creating the coordinates table for the CDS of all coding genes in the genome.####
  raw_gff <- read.delim(paste("../", list.files(path = "../", pattern = "*.gff3"), sep = ""), header = F) #pull the gff file from the parent directory
  cds_table <- raw_gff[raw_gff$V3 == "CDS",] #extract only the CDS rows
  cds_table <- cds_table[, c("V1", "V4", "V5", "V7", "V9")]#remove extra columns
  

#choose which genome, necessary for extracting the correct gene ID from the gene column
  #cds_table$gene_ID <- substr(cds_table[,"V9"], 4, 12) #Add new column with the processed Arabidopsis gene ID
  #cds_table$gene_ID <- substr(cds_table[,"V9"], 4, 15) #Add new column with the processed Brachy gene ID
  #cds_table$gene_ID <- substr(cds_table[,"V9"], 4, 17) #Add new column with the processed Rice gene ID
  #cds_table$gene_ID <- substr(cds_table[,"V9"], 4, 19) #Add new column with the processed Poplar gene ID
  #cds_table$gene_ID <- substr(cds_table[,"V9"], 4, 17) #Add new column with the processed maize gene ID
  #cds_table$gene_ID <- substr(cds_table[,"V9"], 4, 17) #Add new column with the processed tomato gene ID
 
#some gffs need further modification before the script can work on them

  #### #For Physco only, modify CDS entries to reflect the gene they derive from
  #   cds_table$V9 <- str_extract(cds_table$V9, "pacid=.*") # replace each with just the pacid number
  #   cds_table$V9 <- substr(cds_table$V9, 7,14) #just number
  #   mRNA_table <- raw_gff[raw_gff$V3 == "mRNA",] #extract only the mRNA rows
  #   mRNA_table %>% mutate_if(is.factor, as.character) -> mRNA_table #convert to character only
  #   cds_table$gene_ID <- cds_table$V9
  #   #do this after the all_gene_coordinates_table is set up, reduces computational load
  #   #find line with pacID number
  #   for (i in 1:nrow(all_gene_coordinates_table)) {
  #     mRNA_number <- as.numeric(grep(pattern = all_gene_coordinates_table[i,"gene_ID"], mRNA_table$V9, fixed = TRUE)) #find line with matching pacID
  #     all_gene_coordinates_table[i,"gene_ID"] <- str_extract(mRNA_table[mRNA_number, "V9"], pattern = "Parent=.*") %>% str_replace(".*=", "") #pull just the gene ID out from mRNA table and store it as the geneID in the CDS table
  #     print(i)
  #   }
  #   all_gene_coordinates_table <-  all_gene_coordinates_table %>% distinct() %>% arrange(gene_ID)
  #   #save coordinate table for physco because it took a long time to compile
  #   write.table(all_gene_coordinates_table, file = "all_gene_corrdinates_for_physco.tsv")
  
  #For amborella only, convert shorter gene IDs to underscores  
  # cds_table$gene_ID <- substr(cds_table[,"V9"], 4, 43) #Add new column with the processed Amborella gene ID
  #    for (i in 1:nrow(cds_table)) {
  #      last_digit <- substring(cds_table[i,"gene_ID"], 40,40)
  #      second_to_last_digit <- substring(cds_table[i,"gene_ID"], 39,39)
  #      if (last_digit == ".") {
  #        substr(cds_table[i,"gene_ID"], 40,40) <- "_"
  #      } else if (second_to_last_digit == ".") {
  #          substr(cds_table[i,"gene_ID"], 39,40) <- "__"
  #      }
  #     } 

  #for Salagenela only, convert shorter gene IDs to underscores
  # cds_table$gene_ID <- substr(cds_table[,"V9"], 4, 9) #Add new column with the processed Amborella gene ID
  #    for (i in 1:nrow(cds_table)) {
  #      last_digit <- substring(cds_table[i,"gene_ID"], 6,6)
  #      second_to_last_digit <- substring(cds_table[i,"gene_ID"], 5,5)
  #      if (last_digit == ".") {
  #        substr(cds_table[i,"gene_ID"], 6,6) <- "_"
  #      } else if (second_to_last_digit == ".") {
  #          substr(cds_table[i,"gene_ID"], 5,6) <- "__"
  #      }
  #     }
  # cds_table$gene_ID <- paste("Smo_", cds_table$gene_ID, sep = "")#add Smo_ prefix to each gene (my addition for clarity on trees)
    
#polish tables
  cds_table <- cds_table[, c("V1", "V4", "V5", "V7", "gene_ID")] #remove original column that contained the ID
  cds_table$max_coord <- apply(cds_table[,c("V4", "V5")], 1, max) #add column for whichever is max value
  cds_table$min_coord <- apply(cds_table[,c("V4", "V5")], 1, min) #add column for whichever is min value

#create grouped table for consolidating CDS annotations to their master gene
  by_geneID <- cds_table %>% group_by(gene_ID)
  max_table <- distinct(by_geneID %>% filter(max_coord == max(max_coord)))[,c(1,4,5,6)] #filter by gene ID, remove duplicates, and retain only the max value for each gene
  min_table <- distinct(by_geneID %>% filter(min_coord == min(min_coord)))[,c(5,7)]#filter by gene ID, remove duplicates, and retain only the min value for each gene

#make all_gene_coordinates_table with 1 line per gene, and the min and max coordinates for all CDSs found for that gene
  all_gene_coordinates_table <- merge(max_table, min_table, by = "gene_ID", all.x = TRUE)
  colnames(all_gene_coordinates_table) <- c("gene_ID", "Chr", "strand", "max_coord", "min_coord")
  all_gene_coordinates_table <- distinct(all_gene_coordinates_table) #remove duplicate rows
  all_gene_coordinates_table <- all_gene_coordinates_table[,c("gene_ID", "Chr", "strand", "min_coord", "max_coord")] #change order so min is first
  all_gene_coordinates_table <- all_gene_coordinates_table %>% arrange(Chr, min_coord)#sort to make sure gene order is correct for analysis
  ####/-------------------------------------------------/#####
  #end section to create gene coordinates table

  
#  write.table(all_gene_coordinates_table$gene_ID, file = "At_all_genes_IDs.txt", quote=F, row.names = F) #used to get a list of all IDs in genome
  
  
#this section takes in a list of gene IDs and outputs a bed file with coordinates relative to those genes

#grab the list of genes to explore and format it as a dataframe
  gene_list <- data.frame(read.table(file = list.files(pattern = "*.geneIDs.txt"))) %>% distinct() 
  #gene_list <- data.frame(all_gene_coordinates_table[sample(1:nrow(all_gene_coordinates_table), 10000), "gene_ID"]) #control, generate 1000 intergenic regions at random

    colnames(gene_list) <- c("gene_ID")
#  gene_list$"strand" <- NA

#set up blank dataframe to work with
  bed.df <- data.frame("chr" = character(0),"start" = integer(0),"end" = integer(0), "gene" = character(0)) #set up the df to print to the bed file
  bed.df %>% mutate_if(is.factor, as.character) -> bed.df #make sure everything is a character

#generate bed coordinates based on genes from your gene list
  for (i in 1:nrow(gene_list)){
    current_geneID <- as.character(gene_list[i,"gene_ID"]) #look at each gene row by row
    gene_number <- grep(current_geneID, all_gene_coordinates_table$gene_ID, fixed = TRUE)[1] #what line the gene is found in the all gene coordinates table
    
    current_chr <- as.character(all_gene_coordinates_table[gene_number, "Chr"]) #retrieve the Chr number for the gene
    current_strand <- as.character(all_gene_coordinates_table[gene_number, "strand"]) #retrieve the strand for the gene
    gene_list[i,"strand"] <- current_strand #fill in the strandedness on the gene list dataframe
  
  # bed dataframe construction for each of the 4 variants: lrr only, rlk only, and - and + strands
    if (domain_found_version == "LRR_only") {
         if (current_strand == "+") {
          intergenic_start <- all_gene_coordinates_table[gene_number , "max_coord"] 
          intergenic_end <-  all_gene_coordinates_table[gene_number + 1, "max_coord"] 
        } else if (current_strand == "-") {
          intergenic_start <- all_gene_coordinates_table[gene_number - 1, "min_coord"] 
          intergenic_end <- all_gene_coordinates_table[gene_number , "min_coord"] 
        } 
    } else if (domain_found_version == "RLK_only") {
        if (current_strand == "+") {
          intergenic_start <- all_gene_coordinates_table[gene_number - 1, "min_coord"] 
          intergenic_end <-  all_gene_coordinates_table[gene_number, "min_coord"] 
        } else if (current_strand == "-") {
          intergenic_start <- all_gene_coordinates_table[gene_number, "max_coord"] 
          intergenic_end <- all_gene_coordinates_table[gene_number + 1, "max_coord"] 
          }
    }
    
      #remove any with overlapping CDS with another gene. No intergenic space exists
    if (intergenic_end < intergenic_start){
      #Don't write anything to line if irrational intergenic space
    } else if (intergenic_end > intergenic_start + 100000){
      #Don't write anything to line if space is greater than 100,000 nt
    } else {
      # write min and max to dataframe
      bed.df[i,] <- c(current_chr, intergenic_start, intergenic_end, current_geneID)
    }
  }

  
#print all lines to bed file
  bed.df <- na.omit(bed.df) #remove NA rows 
  write.table(bed.df, file = "target_intergenic_space.bed", quote = F, row.names = F, sep = "\t", col.names = F)

#record stats in a log file
  nt_scan_total <- sum(as.numeric(bed.df[,"end"]) - as.numeric(bed.df[,"start"]))
  write(paste("Number of genes in input list: ", nrow(gene_list)), append = TRUE, file = "log.txt")
  write(paste("Number of genes in scan: ", nrow(bed.df)), append = TRUE, file = "log.txt")
  write(paste("Total nucleotides in scan: ", nt_scan_total), append = TRUE, file = "log.txt")
  write(paste("Average length scanned per gene: ", round(nt_scan_total / nrow(bed.df)), " nucleotides"), append = TRUE, file = "log.txt")
  



