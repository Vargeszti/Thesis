setwd("~/Documents/Projects/VargaEszter/Thesis/code")
library(GEOquery)
library(limma)
library(Biobase)


drug_perturb = read.csv('../results/drug_perturb.csv', sep=',', header=TRUE, row.names=1)
for (i in rownames(drug_perturb)[1:3]){
  print(paste0(i,' experiment'))
  geo_id = drug_perturb[i, "geo_id" ]
  ctrl = drug_perturb[i, "ctrl_ids" ]
  stim = drug_perturb[i,  "pert_ids" ]
  sign = drug_perturb[i,  "sign" ]
  fname = i
  ctrl = strsplit(ctrl, '|', fixed = TRUE)[[1]]
  stim = strsplit(stim, '|', fixed = TRUE)[[1]]
  data = getGEO(geo_id, GSEMatrix=TRUE, AnnotGPL = TRUE)
  data = data[[1]]
  eset = exprs(data)
  fil = apply(is.na(eset), 1, sum) == 0
  eset = eset[fil,]
  if (max(eset) > 1000){
    eset = log2(eset+1)
  }
  eset = as.data.frame(eset)
  annot = fData(data)$`Gene symbol`
  annot = annot[fil]
  eset$MEAN = apply(eset, 1, mean)
  eset$GENE = annot
  
  fil=eset$GENE!=''
  eset=eset[fil,]
  
  eset = eset[order(eset$MEAN, decreasing = TRUE),]
  eset = eset[!duplicated(eset$GENE),]
  rownames(eset) = eset$GENE
  n = dim(eset)[2]
  eset = eset[,c(-n, -n+1)]
  
  x = rep(0, length(ctrl) + length(stim))
  names(x) = c(ctrl,stim)
  x[stim] = 1
  
  design=model.matrix(~1 + x)
  rownames(design) = names(x)
  eset = eset[ ,rownames(design)]
  
  fit <- lmFit(eset, design) 
  fit = eBayes(fit)
  results = topTable(fit, coef = 'x', adjust="BH",number = 1000000)
  results$logFC = results$logFC * as.integer(sign)
  results$t = results$t * as.integer(sign)
  write.csv(results, paste0('../results/expression_data/',fname,'.csv'))
}
### first experiment: 
#geo_id = commandArgs(TRUE)[1]
#ctrl = commandArgs(TRUE)[2]
#stim = commandArgs(TRUE)[3]
#sign = commandArgs(TRUE)[4]
#fname = commandArgs(TRUE)[5]

#geo_id = 'GSE4028'
#ctrl = 'GSM92214|GSM92215|GSM92216'
#stim = 'GSM92217|GSM92218|GSM92219|GSM92220'
#sign = '1'
#fname = '1374'


