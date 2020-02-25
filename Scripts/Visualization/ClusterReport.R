"Finds and counts the occurrences of a specific set of genes in every module"

library(ggplot2)

#file - network subdivided in modules
#classF - files with names of genes
#dfile - destination folder where the plots should be saved
ClusterReport <- function(file='.', classF='.', dfile='.') {
  
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
  
  #Set variables
  moduleD = vector() #holds number of samples by module as found in moduleColors
  classD = vector() #holds number of samples by class by module
  
  uniqueModules = unique(moduleColors)
  
  for (i in uniqueModules)
  {
    module = order[moduleColors == i]
    
    moduleD = c(moduleD, length(module))
    uniqueClasses = unique(info$class)
    
    classD = vector()
    for (j in uniqueClasses)
    {
      m2 = info[info$class == j,]
      m2 = m2[m2$gene %in% module,] #returns all genes in module by class
      
      classD = c(classD, length(m2[,1]))
    }
    names(classD) = uniqueClasses
    classD = data.frame(classD)
    l = sprintf('in cluster %s', i)
    classD$name = row.names(classD)
    
    plot <- ggplot(classD, aes(x=name, y=classD)) + 
      geom_bar(stat='identity', fill='lightgreen') +
        geom_text(aes(label=classD), vjust=1.6, color='white', size=3.5) +
          theme(axis.text.x = element_text(angle = -90)) +
            ylab(l)
    
    filepath = file.path(dfile,
                         sprintf('samples by classes in %s', i))
    ggsave(sprintf('%s.eps', filepath), plot=plot, device='eps')
  }
    names(moduleD) = uniqueModules
    moduleD = data.frame(moduleD)
    moduleD$name <- row.names(moduleD)
    
    plot <- ggplot(moduleD, aes(x=name, y=moduleD)) + 
      geom_bar(stat='identity', fill='lightgreen') +
        geom_text(aes(label=moduleD), vjust=1.6, color='white', size=3.5) +
          theme(axis.text.x = element_text(angle = -90)) +
            ylab('genes by modules')
    
    filepath = file.path(dfile,
                         sprintf('samples by modules'))
    ggsave(sprintf('%s.eps', filepath), plot=plot, device='eps')
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
dfile = 'C:\\Users\\Kanan\\Google Drive\\Cancer Research\\Scripts\\Networks\\NormalTissue_Combined.tsvall.RDATA'
