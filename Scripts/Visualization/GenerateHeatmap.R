"
Create heatmaps
"

library(WGCNA)

src = './Data/Networks'
file = 'LUSC_Normal'

file = paste(src, file, sep='/')
nFile = list.files(file)[8]

file = paste(file, nFile, sep='/')

load(file)

TOM = TOM[1:5000, 1:5000]
dissTOM = 1 - TOM

geneTree = hclust(as.dist(dissTOM), method = "average");

minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

diag(dissTOM) = NA;
sizeGrWindow(15,15)

dev.new()

tiff(file="./Heatmap.tiff", compression = 'lzw')

TOMplot(dissTOM^7, geneTree, as.character(dynamicColors))

dev.off()