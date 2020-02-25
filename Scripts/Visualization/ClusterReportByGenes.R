"Finds and counts the occurrences of a specific set of genes in every module and writes them in a CSV file"

library(ggplot2)

ClusterReportByGenes <- function(file='.', classF='.', dfile='.') {
  
  lname = load(file = file)
  
  #Set variabes
  order = unname( unlist(geneTree[3]) )
  names = colnames( adjacency )
  moduleColors
  
  #change order with names
  order = names[order]
  
  #read info and find classes
  #transformInfo returns dataframe
  #gene class
  info = read.delim(classF)
  info = transformInfo(info)
  
  #TODO test ggplot
  #TODO test 38
  #Turn to data.frame first
  
  #Set variables
  moduleD = vector() #holds number of samples by module as found in moduleColors
  classD = vector() #holds number of samples by class by module
  
  uniqueModules = unique(moduleColors)
  
  ModuleClass = data.frame(matrix(NA, ncol=800))
  
  pos = 1
  
  for (i in uniqueModules)
  {
    module = order[moduleColors == i]
    
    moduleD = c(moduleD, length(module))
    uniqueClasses = unique(info$class)
    
    ModuleClass[pos, ] = c(i, rep(NA, 799))
    pos = pos + 1
    
    classD = vector(mode='list')
    for (j in uniqueClasses)
    {
      m2 = info[info$class == j,]
      m2 = m2[m2$gene %in% module,] #returns all genes in module by class
      
      ModuleClass[pos, ] = c(j, rep(NA, 799))
      pos = pos + 1
      ModuleClass[pos, ] = c(as.character(m2$gene), rep(NA, 800 - length(m2$gene)))
      pos = pos + 1
      
    }
    
  }
  
  write.csv(ModuleClass, dfile)
}

transformInfo <- function(info)
{
  newinfo = data.frame()[1:782, 1]
  
  newinfo$gene <- paste(info$Metagene, info$EntrezGene.ID, sep='.')
  newinfo$class <- info$Cell.type
  newinfo <- data.frame(newinfo)
  
  return(newinfo)
}

file = 'C:\\Users\\Kanan\\Google Drive\\Cancer Research\\Scripts\\Networks\\NormalTissue_Combined.tsvall.RDATA\\NormalTissue_Combined.tsvall.RDATA.RDATA'
classF = 'C:\\Users\\Kanan\\Google Drive\\Cancer Research\\Data\\Trajanoski_CellRep_immune_signatures.txt'
dfile = 'C:\\Users\\Kanan\\Google Drive\\Cancer Research\\Scripts\\Networks\\NormalTissue_Combined.tsvall.RDATA\\ByModule.csv'
