"Exports modules to Cytoscape"

library(WGCNA)

src = './Data/Modules/Lusc_Normal'
file = list.files(src)[27]

file = paste(src, file, sep='/') 

load(file)

#TODO add colors
adj = exportNetworkToCytoscape(adjacency, weighted = 0.02)
Tom = exportNetworkToCytoscape(TOM, weighted = 0.02)

#Save adj, TOM
