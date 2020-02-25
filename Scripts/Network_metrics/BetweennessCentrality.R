"Finding betweenness centrality of hub genes"

library(igraph)
library(WGCNA)

convertToIgraph <- function(net) {
  nodes = colnames(net)
  
  #Links
  links = data.frame(fromNode = 0, toNode = 0, weight = 0)
  
  n = dim(net)[1]
  m = dim(net)[2]
  
  j = 1
  rnames = rownames(net)
  for (i in 1:n) {
    node = rnames[i]
    for (k in j:m) {
      if (i == k) next()
      
      toNode = nodes[k]
      weight = net[i, k]
      
      # cat(sprintf('%s:%s %s:%s %s\n', node, i, toNode, j, weight))
      links = rbind(links,
                    data.frame(fromNode = node, toNode = toNode, weight = weight))
      }
    j = j + 1
    }
  
  links = links[-1, ]
  
  nodes = data.frame(nodeNames = nodes)
  
  return(list(links = links, nodes = nodes))
}

hubGenes = read.csv('./HubGenes/TOM^3/Hub Genes.csv')

modules = unique(hubGenes$module)
modules = as.character(modules)

newData = hubGenes[, c(1,2,7)]
newData$betweenness = NA

for (i in modules) {
  genes = hubGenes$gene[hubGenes$module == i]
  genes = as.character(genes)
  
  file = list.files('./Data/Modules/LUSC_Tumor')
  pos = grepl(sprintf(' %s ', i), file)
  file = file[pos]
  file = paste('./Data/Modules/LUSC_Tumor', file, sep = '/')
  
  load(file)
  
  colnames(TOModule) = colnames(adjacencyModule)
  rownames(TOModule) = rownames(adjacencyModule)
  
  geneNames = rownames(adjacencyModule)
  TOModule = TOModule[, geneNames]
  
  diag(TOModule) = 0
  
  #Igraph uses edge strength as cost
  # TOModule = log2(1/TOModule)
  TOModule = 1/TOModule
  
  Tom = convertToIgraph(TOModule)
  net = graph_from_data_frame(d = Tom[[1]], vertices = Tom[[2]], directed = FALSE)
  
  result = betweenness(net, genes, directed = FALSE, normalized = TRUE)
  
  for (i in 1:length(result)) {
    newData[newData$gene == names(result)[i], 4] = result[i]
  }
}

#Visualizing networks
# plot(net,layout=layout.fruchterman.reingold,edge.width=E(net)$weight/2)