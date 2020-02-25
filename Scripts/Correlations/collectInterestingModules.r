"Collects modules with a high correlation and low pval in the highCor DataFrame"

corMod = cbind(ConmoduleTraitCor, OrdmoduleTraitCor, BinmoduleTraitCor)
pvalMod = cbind(ConmoduleTraitPvalue, OrdmoduleTraitPvalue, BinmoduleTraitPvalue)

modNames = rownames(corMod)
colLen = dim(corMod)[[2]]
rowLen = dim(corMod)[[1]]

highCor = data.frame(module = 0, trait = 0, cor = 0, pval = 0)

for (i in 1:colLen) {
  
  
  colName = colnames(corMod)[i]
  for (j in 1:rowLen) {
    name = modNames[j]
    if (abs( corMod[j, i] ) >= 0.195 && pvalMod[j, i] <= 0.05) {
      corVal = corMod[j, i]
      pval = pvalMod[j, i]
      
      highCor = rbind(highCor, c(name, colName, corVal, pval))
    }
  }
  
}