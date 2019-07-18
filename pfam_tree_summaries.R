#Get libraries and setwd----
setwd("/Users/jarrettman/Google Drive/Umass/Bartlett Lab/Shared from Madelaine/CLV_Network/LRR_RLK_Paper/Jarretts_working_stuff/supertree_scratch_JM")

#install.packages("ape")
#install.packages("ggplot2",  dependencies = TRUE)
#install.packages("tidyverse")
#install.packages('rhmmer')
#install.packages("stringr")

library(reshape2)
library(RColorBrewer)
library(ape)
library(ggtree)
library(tidyverse)

library(rhmmer) #required for read_tblout()
library(dplyr)

clades <- c("I", "II", "III_VIIa", "IV","V","VI","VIIb","VIII1","VIII2", "IX", "X","XI_XIIb","XIIa", "XIIIa","XIIIb","XIV", "XV") #need to put all here and put their respective files in working directory.

#Set up the summary table for using in paper
summary_table <- data.frame(matrix(0, ncol = 7, nrow = length(clades)))
colnames(summary_table) <- c("clade", "lrr_domain_only", "kinase_domain_only", "both_domains", "neither_domain", "other_domain", "total")
summary_table$clade <- clades


#set up the supplimental table that has all genes
supp_table1 <- data.frame(matrix(0, ncol = 6, nrow = 0))
colnames(supp_table1) <- c("locus","domains found","LRR domain","Kinase domain","Other domain type","clade")


for(run in clades) {
  #if want to run just a single clade
  #    run <- "XI_XIIb"
  #      input_tree_nwk <- "XI_XIIb.input.align_out.fas.contree"
  
  
  #get the appropriate tree file
    input_tree_nwk <- paste("nwk_trees_and_pfam_calls/",run, ".input.align_out.fas.contree", sep = "")


#get the appropriate pfam table
    raw_pfam_df <- data.frame(read_tblout(paste("nwk_trees_and_pfam_calls/",run, ".pfamout.tsv", sep = "")))

  
  #stuff for buiding a tree easier to read----
    #Read tree but replace branch length on pelle to 1
    temptree <- gsub("DmPelle:\\d.\\d+","DmPelle:1.0", read_file(input_tree_nwk))
    write(temptree, file = "temptree.tre")
    tree <- read.tree("temptree.tre")
    file.remove("temptree.tre")
  
  
  
  #root the tree----
  tree <- root(tree, outgroup = "DmPelle", edgelabel = T, resolve.root = F) #root on Pelle kinase
  
  #comprehensive list of all full gene IDs shown  
  transcripts_in_tree <- data.frame(tree$tip.label) #dataframe for back-converting the pfam gene names after it strips variants.
  colnames(transcripts_in_tree) <- "locus"
  
  
  ##colored genes currently broken, not sure why
  #define tree groups----
  clade_I.group <- as.character(read.table("clades_input_IDs/clade_I.input")[,1])
  clade_II.group <- as.character(read.table("clades_input_IDs/clade_II.input")[,1])
  clade_III.group <- as.character(read.table("clades_input_IDs/clade_III.input")[,1])
  clade_IV.group <- as.character(read.table("clades_input_IDs/clade_IV.input")[,1])
  clade_IX.group <- as.character(read.table("clades_input_IDs/clade_IX.input")[,1])
  clade_V.group <- as.character(read.table("clades_input_IDs/clade_V.input")[,1])
  clade_VI.group <- as.character(read.table("clades_input_IDs/clade_VI.input")[,1])
  clade_VIIa.group <- as.character(read.table("clades_input_IDs/clade_VIIa.input")[,1])
  clade_VIIb.group <- as.character(read.table("clades_input_IDs/clade_VIIb.input")[,1])
  clade_VIII1.group <- as.character(read.table("clades_input_IDs/clade_VIII1.input")[,1])
  clade_VIII2.group <- as.character(read.table("clades_input_IDs/clade_VIII2.input")[,1])
  clade_Xa.group <- as.character(read.table("clades_input_IDs/clade_Xa.input")[,1])
  clade_Xb.group <- as.character(read.table("clades_input_IDs/clade_Xb.input")[,1])
  clade_XI.group <- as.character(read.table("clades_input_IDs/clade_XI.input")[,1])
  clade_XIIa.group <- as.character(read.table("clades_input_IDs/clade_XIIa.input")[,1])
  clade_XIIb.group <- as.character(read.table("clades_input_IDs/clade_XIIb.input")[,1])
  clade_XIIIa.group <- as.character(read.table("clades_input_IDs/clade_XIIIa.input")[,1])
  clade_XIIIb.group <- as.character(read.table("clades_input_IDs/clade_XIIIb.input")[,1])
  clade_XIV.group <- as.character(read.table("clades_input_IDs/clade_XIV.input")[,1])
  clade_XV.group <- as.character(read.table("clades_input_IDs/clade_XV.input")[,1])
  
  #set color table
  #make a new dataframe to organize colors
  color.guide = data.frame(OTU = tree$tip.label, Color = "black", stringsAsFactors =F)
  color.guide[color.guide[,1] %in% clade_I.group,2] = "blue4"
  color.guide[color.guide[,1] %in% clade_II.group,2] = "blue2"
  color.guide[color.guide[,1] %in% clade_III.group,2] = "deeppink"
  color.guide[color.guide[,1] %in% clade_IV.group,2] = "gray"
  color.guide[color.guide[,1] %in% clade_IX.group,2] = "hotpink"
  color.guide[color.guide[,1] %in% clade_V.group,2] = "green2"
  color.guide[color.guide[,1] %in% clade_VI.group,2] = "green4"
  color.guide[color.guide[,1] %in% clade_VIIa.group,2] = "blue"
  color.guide[color.guide[,1] %in% clade_VIIb.group,2] = "green"
  color.guide[color.guide[,1] %in% clade_VIII1.group,2] = "orange"
  color.guide[color.guide[,1] %in% clade_VIII2.group,2] = "darkmagenta"
  color.guide[color.guide[,1] %in% clade_Xa.group,2] = "red"
  color.guide[color.guide[,1] %in% clade_Xb.group,2] = "cyan"
  color.guide[color.guide[,1] %in% clade_XI.group,2] = "purple"
  color.guide[color.guide[,1] %in% clade_XIIa.group,2] = "pink"
  color.guide[color.guide[,1] %in% clade_XIIb.group,2] = "forestgreen"
  color.guide[color.guide[,1] %in% clade_XIIIa.group,2] = "yellow"
  color.guide[color.guide[,1] %in% clade_XIIIb.group,2] = "yellow3"
  color.guide[color.guide[,1] %in% clade_XIV.group,2] = "turquoise1"
  color.guide[color.guide[,1] %in% clade_XV.group,2] = "turquoise4"
  
  
  
  #separate gene name and transcript variant into separate columns
  pfam_df_full <- separate(raw_pfam_df, into = c("locus", "transcript_variant"), col = query_name, sep = "@", fill = "right")
  pfam_df_full[is.na(pfam_df_full)] <- "" #get rid of NAs caused by Smo and AmTr genes that don't have transcript variants
  
  
  #clean unused colums from pfam df and consolidate by gene family,
  pfam_df <- pfam_df_full %>% 
    select(domain_name, locus) %>%
    group_by(locus)  %>% 
    dplyr::summarise("domains found" = paste(domain_name, collapse=", "))  #summarize() is shared by several packages. dplyer::summarize() is necessary to make this work
  

  
  #The pfam domains strip transcript variant from all genes. This adds it back, according to whats on the tree
  pfam_df$locus <- transcripts_in_tree$locus[pmatch(pfam_df$locus, transcripts_in_tree$locus)] #pmatch returns the position in the tree transcripts list that matches each gene from the Pfam table. Then find that postion and return its contents, <- to the pfam table

  
  
  
  
  #create a dataframe of genes that have other types of domains
  other_domains_df <- raw_pfam_df %>%
    dplyr::select(query_name, domain_name) %>% #remove unneeded columns
    dplyr::filter( !grepl('LRR|kinase', `domain_name`))  #remove lines that are either LRR or Kinase calls

    #make a list (held as a dataframe) of all genes that have another domain found  
  other_domains_df <- separate(other_domains_df, into = c("locus", "transcript_variant"), col = query_name, sep = "@") %>% # split by transcript variant
    dplyr::distinct(locus)

    #update gene name to match tree gene name
  other_domains_df$locus <- transcripts_in_tree$locus[pmatch(other_domains_df$locus, transcripts_in_tree$locus)]
  colnames(other_domains_df) <- c("query_name") #make sure column names match up later
  
  
  #add rows to pfam database for genes with no domains called
  domaintable <- merge(transcripts_in_tree, pfam_df, by = "locus", all = T)
  
  #Convert NAs in domain table to blank, otherwise grep errors
  domaintable[is.na(domaintable)] <- ""
  
  #add new columns for describing if specifc domains are found
  domaintable$"LRR domain" <- NA
  domaintable$"Kinase domain" <- NA
  domaintable$"Other domain type" <- NA
  
  
  #create a searchable concatenation of the gene names with other domains for using in search
  genes_with_other_domains <- paste(other_domains_df$query_name, collapse = " ")
  
  
  #change factors to character strings
  domaintable %>% mutate_if(is.factor, as.character) -> domaintable
  
  #Add boolean values to the found variables -
  for (row in 1:nrow(domaintable)) {
    domaintable[row,"LRR domain"] <-  str_detect(domaintable[row, "domains found"], "LRR*")
    domaintable[row,"Kinase domain"] <-  str_detect(domaintable[row, "domains found"], "Pkinase")
    domaintable[row,"Other domain type"] <-  grepl(domaintable[row,"locus"], genes_with_other_domains)
  }
  
  
  #replace all the TRUE and FALSE with Found and Not Found
  domaintable[domaintable=="TRUE"] <- "Found"
  domaintable[domaintable=="FALSE"] <- "Not Found"
  
 
  
  
   #Before removing, make a copy of the domain table for counts at the end
  copy_domain_table <- domaintable
  
    
    
    
  #convert tibble to standard dataframe
  domaintable <- data.frame(domaintable)
  
  #set up the table so that the row names are the gene names rather than numbers
  rownames(domaintable) <- domaintable$locus
  #remove gene name column
  domaintable <- domaintable[,3:5]
  #rename columns for clarity on figure
  colnames(domaintable) <- c("LRR domain","Kinase domain","Other domain type") 
  
  

  
  #decribe the tree
  pfam_tree <- ggtree(tree) + 
    geom_nodelab(aes(label = label),size = 2.5, vjust = .4, hjust = -.15) +
  #  geom_nodelab(aes(label = node), size = 2, vjust = 1.8, hjust = -.15) + 
    ggtitle(paste(run, " clade")) +
    theme(plot.title = element_text(hjust = .5)) +
    geom_treescale(y=-3) +
    geom_tiplab(size=1.8, align=TRUE, hjust = 0, offset = .08, linetype = "dotdash", linesize = .2) +
    scale_x_continuous(expand = expand_scale(mult = c(0, .1))) + 
    scale_y_continuous(expand = c(.08,0))
  
  #render the tree with domains annotated
 gheatmap(pfam_tree, domaintable, offset = .3, width=0.15,  colnames_position = "top", colnames_angle = 35, font.size=2, colnames_offset_x = 0, colnames_offset_y = 0, hjust = 0) + 
    scale_fill_manual(breaks=c("Found", "Not Found"), values=c("blue", "gray")) +
    theme(legend.position = c(.1,.9)) 
  

    #save pdf scaled to 10 genes per inch
 # ggtree::ggsave(plot = last_plot(), filename = paste("supp_table1_and_figures/",run,".tree.pdf",sep = ""), height = length(tree$tip.label)/10+2, limitsize = FALSE)
  ggtree::ggsave(plot = last_plot(), filename = paste("supp_table1_and_figures/",run,".tree.pdf",sep = ""), height = length(tree$tip.label)/10+2, limitsize = FALSE)
  
#### counting section
  #First add all genes and their calls to the supplimental table
  copy_domain_table$clade <- run #modify the domain table to include clade identifier
  supp_table1 <- rbind(supp_table1, copy_domain_table)#add to summary table
  
  #Second, summarize the number of genes with each domain combo
  #these counts need to be reset for each run (clade)
  both_domains <- 0
  lrr_domain_only <- 0
  kinase_domain_only <- 0
  neither_domain <- 0
  other_domain <- sum(copy_domain_table$`Other domain type` == "Found") #this is just a count of other domains
  for (row in 1:nrow(copy_domain_table)) {
    if (copy_domain_table[row,"LRR domain"] == "Found" & copy_domain_table[row,"Kinase domain"] == "Found") {
      both_domains <- both_domains + 1
    } else if (copy_domain_table[row,"LRR domain"] == "Found" & copy_domain_table[row,"Kinase domain"] != "Found") {
      lrr_domain_only <- lrr_domain_only + 1
    } else if (copy_domain_table[row,"LRR domain"] != "Found" & copy_domain_table[row,"Kinase domain"] == "Found") {
      kinase_domain_only <- kinase_domain_only + 1
    } else if (copy_domain_table[row,"LRR domain"] != "Found" & copy_domain_table[row,"Kinase domain"] != "Found") {
      neither_domain <- neither_domain + 1
    }
  }
  #print results of current clade (run) to the summary table
  summary_table[summary_table$clade == run, "lrr_domain_only"] <- lrr_domain_only
  summary_table[summary_table$clade == run, "kinase_domain_only"] <- kinase_domain_only
  summary_table[summary_table$clade == run, "both_domains"] <- both_domains
  summary_table[summary_table$clade == run, "neither_domain"] <- neither_domain
  summary_table[summary_table$clade == run, "other_domain"] <- other_domain
  summary_table[summary_table$clade == run, "total"] <- nrow(domaintable)
  

}#end of run loop
  


write.csv(supp_table1, file = paste("supp_table1_and_figures/","supp_table1.csv", sep = ""))




##below this are for additional figures



#Plot the members from each clade and their domains
ggplot() + 
  geom_bar(data = melt(summary_table[,1:5]), aes(x = clade, y = value, fill = variable),stat = 'identity') + 
  scale_fill_brewer("gene types",palette = "Set2") +
  labs(y = "# genes found" ) +
  geom_col(data = summary_table[,c(1,6)], aes(clade, other_domain), width = .25, alpha = .5 ) + 
  theme(legend.position = c(.4,.9)) 

  



#also plot how many genes each species has
Arabidopsis <- sum(grepl(pattern = "AT.*", supp_table1[,"locus" ]))
Brachy <- sum(grepl(pattern = "Bradi.*", supp_table1[,"locus" ]))
Amborella <- sum(grepl(pattern = "evm.*", supp_table1[,"locus" ]))
Rice <- sum(grepl(pattern = "LOC_Os.*", supp_table1[,"locus" ]))
Poplar <- sum(grepl(pattern = "Potri.*", supp_table1[,"locus" ]))
Physco <- sum(grepl(pattern = "Pp.*", supp_table1[,"locus" ]))
Salaginella <- sum(grepl(pattern = "Smo.*", supp_table1[,"locus" ]))
Tomato <- sum(grepl(pattern = "Solyc.*", supp_table1[,"locus" ]))
Maize <- sum(grepl(pattern = "Zm.*", supp_table1[,"locus" ]))

species_list <- c("Arabidopsis", "Brachy", "Amborella", "Rice", "Poplar", "Physco", "Salaginella", "Tomato", "Maize" )

genes_per_species_df <- data.frame(matrix(0, ncol = 2, nrow = length(species_list)))
#rownames(genes_per_species_df) <- species_list
colnames(genes_per_species_df) <- c("species", "genes")
genes_per_species_df$species <- species_list
genes_per_species_df[1,2] <- Arabidopsis
genes_per_species_df[2,2] <- Brachy
genes_per_species_df[3,2] <- Amborella
genes_per_species_df[4,2] <- Rice
genes_per_species_df[5,2] <- Poplar
genes_per_species_df[6,2] <- Physco
genes_per_species_df[7,2] <- Salaginella
genes_per_species_df[8,2] <- Tomato
genes_per_species_df[9,2] <- Maize

genes_per_species_df <- genes_per_species_df[order(genes_per_species_df$genes),]

ggplot(genes_per_species_df, aes(x=1, y=genes, fill=species)) +  geom_bar(stat="identity") 
