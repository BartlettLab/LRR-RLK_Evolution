library(dplyr)


tbl.raw <- read.csv("/Users/jarrettman/Google Drive/Umass/Bartlett Lab/Shared from Madelaine/CLV_Network/LRR_RLK_Paper/Jarretts_working_stuff/supertree_scratch_JM/supp_table1_and_figures/supp_table1.csv")

tbl.full <- tbl.raw %>%
  select(locus,LRR.domain,Kinase.domain,Other.domain.type,clade)

#remove Pelle rows
tbl <- tbl.full[which(tbl.full$locus != "DmPelle"),]

#total_genes
nrow(tbl)

#regular cannonical LRR-RLks
nrow(tbl[which(tbl$LRR.domain == "Found" & tbl$Kinase.domain == "Found" & tbl$Other.domain.type == "Not Found"),])
nrow(tbl[which(tbl$LRR.domain == "Found" & tbl$Kinase.domain == "Found" & tbl$Other.domain.type == "Not Found"),]) / nrow(tbl) #percent of total

#canonical LRR-RLKs without testing for other domains in clades I and VIII-2
nrow(tbl[which(tbl$LRR.domain == "Found" & tbl$Kinase.domain == "Found" & tbl$clade == "I"),]) + nrow(tbl[which(tbl$LRR.domain == "Found" & tbl$Kinase.domain == "Found" & tbl$clade == "VIII2"),]) +  #canonical genes in clades I and VIII2
 nrow(tbl[which(tbl$LRR.domain == "Found" & tbl$Kinase.domain == "Found" & tbl$Other.domain.type == "Not Found" & tbl$clade != "I" & tbl$clade != "VIII2"),]) #canonical genes found outside clades I and VIII2  
(nrow(tbl[which(tbl$LRR.domain == "Found" & tbl$Kinase.domain == "Found" & tbl$clade == "I"),]) + nrow(tbl[which(tbl$LRR.domain == "Found" & tbl$Kinase.domain == "Found" & tbl$clade == "VIII2"),]) +  #canonical genes in clades I and VIII2
  nrow(tbl[which(tbl$LRR.domain == "Found" & tbl$Kinase.domain == "Found" & tbl$Other.domain.type == "Not Found" & tbl$clade != "I" & tbl$clade != "VIII2"),])) / nrow(tbl) #percent of total


nrow(tbl[which(tbl$LRR.domain == "Found" & tbl$Kinase.domain == "Found" & tbl$Other.domain.type == "Not Found"),]) / nrow(tbl) #percent of total


#LRR only genes
nrow(tbl[which(tbl$LRR.domain == "Found" & tbl$Kinase.domain == "Not Found" & tbl$Other.domain.type == "Not Found"),])
nrow(tbl[which(tbl$LRR.domain == "Found" & tbl$Kinase.domain == "Not Found" & tbl$Other.domain.type == "Not Found"),]) / nrow(tbl) #percent of total


#kinase loss genes
nrow(tbl[which(tbl$LRR.domain == "Found" & tbl$Kinase.domain == "Not Found"),])
nrow(tbl[which(tbl$LRR.domain == "Found" & tbl$Kinase.domain == "Not Found"),]) / nrow(tbl) #percent of total

#kinase only genes
nrow(tbl[which(tbl$LRR.domain == "Not Found" & tbl$Kinase.domain == "Found" & tbl$Other.domain.type == "Not Found"),])
nrow(tbl[which(tbl$LRR.domain == "Not Found" & tbl$Kinase.domain == "Found" & tbl$Other.domain.type == "Not Found"),]) / nrow(tbl) #percent of total

#LRR loss genes
nrow(tbl[which(tbl$LRR.domain == "Not Found" & tbl$Kinase.domain == "Found"),])
nrow(tbl[which(tbl$LRR.domain == "Not Found" & tbl$Kinase.domain == "Found"),]) / nrow(tbl) #percent of total

#genes with other domain type
nrow(tbl[which(tbl$Other.domain.type == "Found"),])
nrow(tbl[which(tbl$Other.domain.type == "Found"),]) / nrow(tbl) #percent of total

#genes with other domain types (excluding clades I and VII-2)
nrow(tbl[which(tbl$Other.domain.type == "Found" & tbl$clade != "I" & tbl$clade != "VIII2" ),])
nrow(tbl[which( tbl$clade != "I" & tbl$clade != "VIII2" ),]) #total genes outside this clade
nrow(tbl[which(tbl$Other.domain.type == "Found" & tbl$clade != "I" & tbl$clade != "VIII2" ),]) / nrow(tbl[which( tbl$clade != "I" & tbl$clade != "VIII2" ),]) #percent of total

#genes with neither LRR or Kinase domain
nrow(tbl[which(tbl$LRR.domain == "Not Found" & tbl$Kinase.domain == "Not Found") ,])
nrow(tbl[which(tbl$LRR.domain == "Not Found" & tbl$Kinase.domain == "Not Found") ,]) / nrow(tbl) #percent of total

#gene fission pairs
19  / nrow(tbl) #percent of total
