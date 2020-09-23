setwd("~/Documents/Projects/VargaEszter/Thesis/code")

par1= commandArgs(TRUE)[1]
par2 = commandArgs(TRUE)[2]

print(par1)
print(par2)
write.csv(par1,'../results/dummy.csv')