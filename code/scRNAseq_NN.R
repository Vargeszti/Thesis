setwd("../Desktop/MSc/Szakdoga/Thesis/code/")
library(nichenetr)
library(tidyverse)

#Read expression data
scrna_expression = gene_perturb = read.csv('../data/scrnaseq/test_mat.csv', sep=',', header=TRUE, row.names=1)
#expression = scrna_expression$expression
sample_info = scrna_expression$sample_info # contains meta-information about the cells

#which genes are expressed
tumors_remove = c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")

CAF_ids = sample_info %>% filter(`Lymph node` == 0) %>% filter((tumor %in% tumors_remove == FALSE)) %>% filter(`non-cancer cell type` == "CAF") %>% .$cell
malignant_ids = sample_info %>% filter(`Lymph node` == 0) %>% filter(`classified  as cancer cell` == 1) %>% filter((tumor %in% tumors_remove == FALSE)) %>% .$cell

expressed_genes_CAFs = expression[CAF_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_malignant = expression[malignant_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()

#Load the ligand-target model
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5]# target genes in rows, ligands in columns

#NichNet single-cell ligand activity analysis

#1)defining potentially active ligands (expressed in CAFs or bind to cancer cell receptors)
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
ligands = lr_network$from %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_CAFs)
receptors = lr_network$to %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_malignant)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% .$from %>% unique()
head(potential_ligands)

#2)scaling the data
background_expressed_genes = expressed_genes_malignant %>% .[. %in% rownames(ligand_target_matrix)]
expression_scaled = expression %>% .[malignant_ids,background_expressed_genes] %>% scale_quantile()

#3)Ligand activity analysis on 10 cells from HN5 tumor
malignant_hn5_ids = sample_info %>% filter(tumor == "HN5") %>% filter(`Lymph node` == 0) %>% filter(`classified  as cancer cell` == 1)  %>% .$cell %>% head(10)
ligand_activities = predict_single_cell_ligand_activities(cell_ids = malignant_hn5_ids,expression_scaled = expression_scaled,
                                                          ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
#Ligand prioritization by regression analysis

#1)scoring malignant cells on expression of "TGFBI"
cell_scores_tbl = tibble(cell = malignant_hn5_ids, score = expression_scaled[malignant_hn5_ids,"TGFBI"])
#2)normalizing the ligand activities with modZ to make diff. cells comparable
normalized_ligand_activities = normalize_single_cell_ligand_activities(ligand_activities)
#3)perform correlation and regression analysis
output_correlation_analysis = single_ligand_activity_score_regression(normalized_ligand_activities,cell_scores_tbl)
output_correlation_analysis %>% arrange(-pearson_regression) %>% select(pearson_regression, ligand)
#visualizing relation between lig.act+cell property score
inner_join(cell_scores_tbl,normalized_ligand_activities) %>% ggplot(aes(score,TNC)) + geom_point() + geom_smooth(method = "lm")
