setwd("../Desktop/MSc/Szakdoga/Thesis/code/")
library(GEOquery)
library(limma)
require(Biobase)

### second experiment: 
geo_id = 'GSE14411'
ctrl = 'GSM360098|GSM360099|GSM360100'
stim = 'GSM360101|GSM360102|GSM360103'
ctrl = strsplit(ctrl, '|', fixed = TRUE)[[1]]
stim = strsplit(stim, '|', fixed = TRUE)[[1]]
data = getGEO(geo_id, GSEMatrix=TRUE, AnnotGPL = TRUE)
data = data[[1]]
eset = exprs(data) #eset=expession set
eset = as.data.frame(eset)
annot = fData(data)$`Gene symbol`
eset$MEAN = apply(eset, 1, mean)
eset$GENE = annot
eset = eset[order(eset$MEAN, decreasing = TRUE),]
eset = eset[!duplicated(eset$GENE),]
rownames(eset) = eset$GENE
n = dim(eset)[2] 
eset = eset[,c(-n, -n+1)] #utolsó 2 sort levágja

x = rep(0, length(ctrl) + length(stim))
names(x) = c(ctrl,stim)
x[stim] = 1

design=model.matrix(~1 + x)
rownames(design) = names(x)
eset = eset[ ,rownames(design)]

#Estimate the fold changes and standard errors by fitting 
  #a linear model for each gene. The design matrix indicates which arrays are dye-swaps
fit <- lmFit(eset, design) #limma t test correction
#Apply empirical Bayes smoothing to the standard errors
fit = eBayes(fit)
#Show statistics for the top 10 genes.
results = topTable(fit, coef = 'x', adjust="BH",number = 1000000)
write.csv(results, '../results/GSE14411.csv')
