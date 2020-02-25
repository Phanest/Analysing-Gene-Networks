"Reads and calculates module metrics"

library(WGCNA)
library(igraph)

Centralization <- function(density, n, kMax) {
  
  a = n/(n-2)
  
  b = kMax/(n-1)
  b = b - density
  
  return(a*b)
}

Heterogeneity <- function(degree) {
  
  m = mean(degree)
  v = var(degree)
  
  return(sqrt(v)/m)
}

density <- function(degree) {
  
  m = mean(degree)
  n = length(degree)
  
  return(m/(n-1))
}

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

name = 'LUSC_Tumor'

#Load module eigengenes

load('./Data/Module Eigengenes/LUSC_Tumor.RDATA')
MEadj = cor(MEs)

names = colnames(MEadj)

metrics = data.frame(module = 0,
                     mean = 0,
                     var = 0,
                     diameter = 0,
                     density = 0,
                     centrality = 0,
                     ClusterCoef = 0,
                     Heterogeneity = 0,
                     betweenness = 0
)

net = MEadj
diag(net) = 0

size = dim(net)[1]

#Make Igraph Network
adj = convertToIgraph(abs(1/net))
adjNet = graph_from_data_frame(d = adj[[1]], vertices = adj[[2]], directed = FALSE)

ClusterCoef = clusterCoef(abs(net))
degree = rowSums(net)

d = diameter(adjNet)
#d = diameter(TomNet) #Calculate diameter using TOM network

den = density(degree)

kMax = max(degree)
centrality = Centralization(den, size, kMax)

heterogeneity = Heterogeneity(degree)

result = betweenness(adjNet, directed = FALSE, normalized = TRUE)

