"Merges highly correlated modules into meta-modules"

library(WGCNA)

load('./Data/Modules/Module_Genes.RDATA')

file = list.files('./Data/Clusters RDATA')[4]
file = paste('./Data/Clusters RDATA', file, sep = '/')
load(file)

n = dim(datExpr)[2]
colors = rep(x = 'grey', n)

#Meta-Genes list
metaModules = list()
metaModules[[1]] = c('purple', 'lightyellow', 'greenyellow', 'black', 'blue', 'plum2', 'turquoise')
metaModules[[2]] = c('paleturquoise', 'brown4', 'maroon', 'darkgrey', 'skyblue3', 'darkturquoise', 'darkorange2', 'darkgreen', 'cyan', 'brown', 'green', 'saddlebrown', 'magenta', 'red')
metaModules[[3]] = c('royalblue', 'tan', 'yellowgreen', 'palevioletred3', 'darkolivegreen', 'grey60')
metaModules[[4]] = c('pink', 'darkorange', 'lightgreen', 'salmon', 'yellow', 'orangered4', 'thistle2', 'violet')

for (i in 1:length(metaModules)) {
  modules = metaModules[[i]]
  name = modules[1]
  
  genes = c()
  for (module in modules) {
    mGenes = moduleGenes[[module]]
    
    genes = c(genes, mGenes)
  }
  
  pos = which(colnames(datExpr) %in% genes)
  
  colors[pos] = name
}

eigengenes = moduleEigengenes(datExpr, colors = colors, excludeGrey = TRUE)