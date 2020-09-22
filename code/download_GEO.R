setwd("~/Documents/Projects/VargaEszter/Thesis/code")
library(GEOquery)
library(limma)

### first experiment: 
geo_id = 'GSE32316'
ctrl = 'GSM800590|GSM800591|GSM800592|GSM800593|GSM800594|GSM800595'
stim = 'GSM800596|GSM800597|GSM800598|GSM800599|GSM800600|GSM800601'
ctrl = strsplit(ctrl, '|', fixed = TRUE)[[1]]
stim = strsplit(stim, '|', fixed = TRUE)[[1]]
data = getGEO(geo_id, GSEMatrix=TRUE, AnnotGPL = TRUE)
data = data[[1]]
eset = exprs(data)
eset = as.data.frame(eset)
annot = fData(data)$`Gene symbol`
eset$MEAN = apply(eset, 1, mean)
eset$GENE = annot
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
write.csv(results, '../results/GSE32316.csv')
