setwd("C:/Users/nidik/Desktop/MSc/Szakdoga/Thesis/code")
#only for first installation
#install.packages("devtools")
#install.packages("BiocManager")
#BiocManager::install("multtest")
#BiocManager::install("limma")
#devtools::install_github("saeyslab/nichenetr")
library(nichenetr)
#install.packages("tidyverse")
library(tidyverse)
rm(list=ls()) 
#Getting the NicheNet matrix
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

write.csv(ligand_target_matrix,'../data/ligand_target_matrix.csv')
write.csv(lr_network,'../data/lr_network.csv')


