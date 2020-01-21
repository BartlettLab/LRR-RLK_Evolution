# LRR-RLK_Evolution
code for LRR-RLK evolution paper

backbone_tree_construction_191125.R - R script for building backbone tree from RAxML run

backbone_color_sets/ - text files containing sets of genes, use with backbone_tree_construction.R

input_taxa_collect.sh - script for collecting and creating an alignment of known prior genes in a clade

blast_eval_finder.sh - script for determining the typical E-values for priors in a clade subalignment

hmm_blast_search.sh - script for ranking all genes from a database of sequences by both HMM and BLAST search strength to a clade subalignment, and thresholding the list according to the typical E-values found in blast_eval_finder.sh

run_iqtree.sh - script for constructing a phylogenetic tree based on the consolidated hit from hmm_blast_search.sh

figure_tree_generate.sh - script for both constructing a phylogenetic tree based on the consolidated hit from hmm_blast_search.sh
as well as performing a Pfam domain scan on all results

pfam_tree_summaries.R - script for rendering a tree from the output of figure_tree_generate.sh that analyszes the results from all Pfam domain scans and annotates them on the tree

Convert_Zm_to_longest_transcript.R - script to generate a database of only the longest transcript variant for each maize gene

expression_correlation.R - script to perform expression correlation tests on two genes

Find_number_of_variants.R - script to generate tables of gene structural variants by clade or species

flanking_region_inspect.sh - script to run pfam domain scans on translated nucleotide sequence flanking genes with apparent truncations

bed_generator.R - R script to find the genomic coordinates for a possible cryptic domain near genes with apparent truncations
