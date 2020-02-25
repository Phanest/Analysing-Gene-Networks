
source('./correlateTraits.R')

name = 'LUSC_Tumor'

#READ Gene expression data (datExpr)

file = list.files('./Data/Clusters RDATA/')[4]
file = paste('./Data/Clusters RDATA', file, sep = '/')
load(file)

#Read datKME

load('./Data/datakME.RDATA')

#Read Clinical Data
library('RSQLite')

con <- dbConnect( RSQLite::SQLite(), dbname="./cancer")
query <- NULL

query <- sprintf("SELECT * FROM patient")
clinicalData = dbGetQuery( con, query)

file = name

file = sprintf('./Data/Module Eigengenes/%s.RDATA', file)
load(file)

ids = clinicalData$barcode %in% names
clinicalData = clinicalData[ids, ]

rownames(clinicalData) = clinicalData$barcode

clinicalData = nominal2numeric2(clinicalData, 'LUSC_Tumor_HubGenes')

ContinousTraits = c('clinical_data_years_to_birth','clinical_data_number_of_lymph_nodes',
                    'clinical_data_days_to_death', 'clinical_data_tumor_size', 'clinical_data_karnofsky_performance_score',
                    'clinical_data_psa_value')
OrdinalTraits = c('clinical_data_pathologic_stage', 'clinical_data_pathology_T_stage', 'clinical_data_pathology_M_stage',
                  'clinical_data_pathology_N_stage', 'clinical_data_year_of_tobacco_smoking_onset', 'clinical_data_gleason_score',
                  'clinical_data_extrathyroidal_extension', 'ips_ctla4_neg_pd1_neg', 'ips_ctla4_neg_pd1_pos',
                  'ips_ctla4_pos_pd1_neg', 'ips_ctla4_pos_pd1_pos')
BinaryTraits = c('clinical_data_gender')

#Read interesting modules
intModules = read.csv('./interestingModules.csv')
intModules = intModules[, -1]
intModules$module = sapply(intModules$module, substring, 3)

hubGeneSignificanceCollection = c()

hubGenes = data.frame(gene = 0, trait = 0, GS = 0, pVal = 0,
                      moduleMembership = 0, degree = 0, module = 0, 
                      hubGeneSignificance = 0, moduleSignificance = 0) 

for (i in 1:length(intModules$module)) {
  
  moduleColor = intModules[ i, 1]
  intTrait = intModules[ i, 2]
  intTrait = trimws(intTrait)
  
  file = list.files('./Data/Modules/LUSC_Tumor')
  file = paste('./Data/Modules/LUSC_Tumor', file[which(grepl(sprintf(' %s ', moduleColor), file) == TRUE)], sep = '/' )
  
  load(file)
  
  colnames(TOModule) = colnames(adjacencyModule)
  rownames(TOModule) = rownames(adjacencyModule)
  
  geneNames = rownames(adjacencyModule)
  adjacencyModule = adjacencyModule[, geneNames]
  TOModule = TOModule[, geneNames]
  
  diag(adjacencyModule) = 0
  diag(TOModule) = 0
  
  GS = NA
  pval = NA
  if (intTrait %in% ContinousTraits) {
    GS = cor(datExpr[, geneNames], clinicalData[, intTrait], use = 'p', method = 'p')
    pval = corPvalueStudent(GS, dim(datExpr)[[1]])
  }
  if (intTrait %in% OrdinalTraits) {
    GS = cor(datExpr[, geneNames], clinicalData[, intTrait], use = 'p', method = 's')
    pval = corPvalueStudent(GS, dim(datExpr)[[1]])
  }
  if (intTrait %in% BinaryTraits) {
    GS = biserialCor(datExpr[, geneNames], clinicalData[, intTrait])
    pval = GS[[2]]
    GS = GS [[1]]
  }
  #TODO try to shoot yourself in the foot by decimating GS using a Beta power of 3
  GS = GS^3
  
  dir.create( sprintf('./HubGenes/%s - %s', moduleColor, intTrait), showWarnings = FALSE )

  sizeGrWindow(15,15)

  dev.new()

  tiff(file= sprintf("./HubGenes/%s - %s/Degree_GS for %s.tiff", moduleColor, intTrait, moduleColor),
      width = 5, height = 5, units = 'in',
      res = 500, compression = 'lzw')
  
  #TODO
  degree = rowSums(abs(TOModule))
  
  verboseScatterplot(degree,
                     abs(GS),
                     col = moduleColor,
                     xlab = 'Connectivity',
                     ylab = sprintf('Gene Significance %s', intTrait),
                     abline = TRUE
                     )

  dev.off()
  
  modkME = datakME[geneNames, moduleColor]
  names(modkME) = geneNames
  
  
  sizeGrWindow(15,15)

  dev.new()

  tiff(file= sprintf("./HubGenes/%s - %s/MM_GS for %s.tiff", moduleColor, intTrait, moduleColor),
       width = 5, height = 5, units = 'in',
       res = 500, compression = 'lzw')

  verboseScatterplot(abs(modkME),
                     abs(GS),
                     col = moduleColor,
                     xlab = 'Module Membership',
                     ylab = sprintf('Gene Significance %s', intTrait),
                     abline = TRUE
  )

  dev.off()
  
  
  sizeGrWindow(15,15)

  dev.new()

  tiff(file= sprintf("./HubGenes/%s - %s/MM_Degree for %s.tiff", moduleColor, intTrait, moduleColor),
       width = 5, height = 5, units = 'in',
       res = 500, compression = 'lzw')

  verboseScatterplot(abs(modkME),
                     degree,
                     col = moduleColor,
                     xlab = 'Module Membership',
                     ylab = sprintf('Connectivity', intTrait),
                     abline = TRUE
  )

  dev.off()

  sizeGrWindow(15,15)

  dev.new()

  tiff(file= sprintf("./HubGenes/%s - %s/MM_Degree^3 for %s.tiff", moduleColor, intTrait, moduleColor),
       width = 5, height = 5, units = 'in',
       res = 500, compression = 'lzw')

  verboseScatterplot(abs(modkME)^3,
                     degree,
                     col = moduleColor,
                     xlab = 'Module Membership',
                     ylab = sprintf('Connectivity', intTrait),
                     abline = TRUE
  )

  dev.off()
  
  #Hub Genes
  
  hGS = rownames(GS)[abs(GS) >= 0.2] #Selection criteria
  hkMM = names(modkME)[abs(modkME) >= 0.8] #Selection criteria
  
  intGenes = intersect(hGS, hkMM)
  
  #Module Significance
  MS = mean(abs(GS))
  
  hubGeneSignificance = 0
  if (moduleColor %in% names(hubGeneSignificanceCollection)) hubGeneSignificance = hubGeneSignificanceCollection[[moduleColor]]
  else {
    hubGeneSignificance = sum(abs(GS)*degree)/sum(degree^2)
    hubGeneSignificanceCollection[moduleColor] = hubGeneSignificance
  }
  
  if (length(intGenes) > 0) {
    hubGenes = rbind(hubGenes, data.frame(gene = intGenes,
                                          trait = intTrait,
                                          GS = GS[intGenes, ],
                                          pVal = pval[intGenes, ],
                                          moduleMembership = modkME[intGenes],
                                          degree = degree[intGenes],
                                          module = rep(x = moduleColor, length(intGenes)),
                                          hubGeneSignificance = rep(x = hubGeneSignificance, length(intGenes)),
                                          moduleSignificance = rep(x = MS, length(intGenes)) )
                     )
  }
 
  while (!is.null(dev.list()))  dev.off()
  
}

hubGenes = hubGenes[-1, ]
write.csv(hubGenes, './HubGenes/Hub Genes.csv', row.names = FALSE)

print('Complete!')