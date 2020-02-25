"Groups modules by ZSummary and extracts the gene names"

Data = read.csv('./Data/Metrics LUSC_Tumor.csv')

#Filter Preserved
Filter = Data[Data[, 10] > 10, ]
modules = Filter$module
modules = as.character(modules)

dir.create('./Enrichment/Preserved')

for (i in modules) {
  
  file = list.files('./Data/Modules/LUSC_Tumor')
  pos = grepl(sprintf(' %s ', i), file)
  file = file[pos]
  file = paste('./Data/Modules/LUSC_Tumor', file, sep = '/')
  
  load(file)
  
  geneNames = rownames(adjacencyModule)
  geneNames = as.character(geneNames)
  geneNames = sapply(geneNames, function (x) strsplit(x, split = '\\.')[[1]][1])
  geneNames = unname(geneNames)
  
  write.csv(geneNames, sprintf('./Enrichment/Preserved/%s.csv', i))
}