"Fragments large network into module sized pieces"

src = './Data/Networks'
folder = 'THCA_Normal'

#Load Network
src = paste(src, folder, sep = '/')

file = list.files(src)[2]

file = paste(src, file, sep = '/')

load(file)

colors = unique(moduleColors)

#Load datTraits
src = './Data/Clusters RDATA'
file = list.files(src)[2]

file = paste(src, file, sep = '/')

load(file)

names = rownames(datTraits)

for (i in colors) {
  iModule = moduleColors == i
  
  newColor = i
  
  adjacencyModule = adjacency[iModule, ]
  TOModule = TOM[iModule, ]
  
  save(MEs, names, newColor, adjacencyModule, TOModule,
       file = sprintf('Network Module %s for %s.RDATA', i, folder))
  
  #old save
  #save(MEs, moduleLabels, moduleColors, geneTree, adjacency, TOM,
  #     file = sprintf("./Networks/%s.RData", file))
}