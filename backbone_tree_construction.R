#code for building LRR-RLK backbone

#Backbone was constructed in CIPRES (https://www.phylo.org/) using RAxML with the following call:
#raxmlHPC-HYBRID -T 4 -n RepSeq_190610 -s infile.txt -q part.txt -p 12345 -m PROTGAMMALG -g constraint.tre -f a -N 1000 -x 12345 --asc-corr lewis

#change to working directory with "RAxML_bipartitions.*" file

#load libraries
library(ape)
library(ggtree)

#read in file
tree <- read.tree("RAxML_bipartitions.*")#your file name here

#root tree
#In order for bootstrap to track branches correctly when rerooted, the edgelabel=TRUE option is vital
root <- "Dm_pelle"
tree <- root(tree, outgroup = root, resolve.root = T, edgelabel=TRUE)

#write out tree in Newick format
write.tree(tree,file = "backbone_rooted.tre")

#metadata for tips
#sets were assigned based on (1) clades and (2) constraint groups
tree.guide = data.frame(OTU = tree$tip.label, Color = 0, Shape = 1, stringsAsFactors =F)
colorSetFiles <- list.files(path = "ColorSets",pattern = "clade")
colorSetList <- list()

#This loop reads in lists of sets to be grouped together and assigned values for clades (Color) and constraint groups (Shape)
for(i in 1:length(colorSetFiles)){
	temp_vec <- as.character(read.table(paste0("ColorSets/",colorSetFiles[i]),header=F)[,1])#these are included in github
	cladeNum <- gsub(colorSetFiles[i],pattern = "_set.*",rep = "")
	tree.guide$Color[tree.guide$OTU %in% temp_vec] <- cladeNum
	setNum <- gsub(gsub(colorSetFiles[i],pattern = ".*set",rep = ""),pat = ".txt", rep = "")
	tree.guide$Shape[tree.guide$OTU %in% temp_vec] <- setNum
	colorSetList[[colorSetFiles[i]]] = temp_vec
}
tree.guide$Color[tree.guide$Color == 0] = "OG"

#make clades (Color) a factor
factor(tree.guide$Color)

#make a tree with boostrap supports and no tip labels
p <- ggtree(tree, size = 0.3) + geom_text2(aes(subset=!isTip, label=label), size = 1.5, hjust=1.25, vjust = -0.4)

#make tree with tip labels (i.e. supplementary fig. S18)
p2 <- p %<+% tree.guide + geom_tiplab(size = 0.75)
plot(p2)

#write out tree to PDF
pdf("backbone_with_names.pdf",width = 8.5, height = 11)
p2
dev.off()

#add tips with color and shape based on clade and constraint group
p3 <- p %<+% tree.guide + geom_tippoint(aes(color = Color,shape = Shape),size = 1)
plot(p3)

#here, we collapse constrained branches to focus on deep nodes
#first, we draw the tree with node numbers to establish which nodes to collapse
p4 <- ggtree(tree, size = 0.3) + geom_text2(aes(subset=!isTip, label=node), size = 1.5, hjust=1.25, vjust = -0.4)
plot(p4)

#These were the nodes we collapsed based on comparing trees p3 (colored nodes) and p4 (node numbers)
p5 <- p %>% collapse(node=458)
p5 <- p5 %>% collapse(node = 485)
p5 <- p5 %>% collapse(node = 453)
p5 <- p5 %>% collapse(node = 445)
p5 <- p5 %>% collapse(node = 436)
p5 <- p5 %>% collapse(node = 409)
p5 <- p5 %>% collapse(node = 427)
p5 <- p5 %>% collapse(node = 406)
p5 <- p5 %>% collapse(node = 489)
p5 <- p5 %>% collapse(node = 540)
p5 <- p5 %>% collapse(node = 532)
p5 <- p5 %>% collapse(node = 369)
p5 <- p5 %>% collapse(node = 548)
p5 <- p5 %>% collapse(node = 338)
p5 <- p5 %>% collapse(node = 353)
p5 <- p5 %>% collapse(node = 328)
p5 <- p5 %>% collapse(node = 304)
p5 <- p5 %>% collapse(node = 306)
p5 <- p5 %>% collapse(node = 316)
p5 <- p5 %>% collapse(node = 557)
p5 <- p5 %>% collapse(node = 573)
p5 <- p5 %>% collapse(node = 572)
p5 <- p5 %>% collapse(node = 570)
p5 <- p5 %>% collapse(node = 567)
p5 <- p5 %>% collapse(node = 574)

#write out tree
pdf("backbone_collapsed_branches.pdf",width = 8.5, height = 11)
p5
dev.off()
