"Correlates traits (continous, ordinal, binary) to MEs"

library(WGCNA)
library(ltm)

source('./RemoveMissingTraits.R')

replaceValues <- function(column, findValues, changeValues)
{
  
  for(i in 1:length( findValues ) )
  {
    pos = which( column == findValues[i] )
    column[pos] = changeValues[i]
  }
  
  return(column)
}

nominal2numeric2 <- function(dataset, name)
{
  
  a <- c('m0', 'm1')
  b <- c('n0', 'n1', 'n2', 'n3')
  c <- c('t0', 'tis', 't1', 't2', 't3', 't4')
  d <- c('i', 'ii', 'iii', 'iv') #also Clinical Stage
  e <- c("female", "male")
  
  
  if ('clinical_data_pathology_M_stage' %in% colnames(dataset)) {
    
    dataset$clinical_data_pathology_M_stage[dataset$clinical_data_pathology_M_stage == 'mx'] = NA
    
    dataset$clinical_data_pathology_M_stage <- as.numeric (factor(
      replaceValues(dataset$clinical_data_pathology_M_stage, 
                    findMStage, 
                    changeMStage),
      levels= a) )
  
  }
  
  dataset$clinical_data_pathology_N_stage[dataset$clinical_data_pathology_N_stage == 'nx'] = NA
  
  dataset$clinical_data_pathology_N_stage <- as.numeric( factor(
    replaceValues(dataset$clinical_data_pathology_N_stage, 
                  findNStage, 
                  changeNStage),
    levels = b) )
  
  dataset$clinical_data_pathology_T_stage[dataset$clinical_data_pathology_T_stage == 'tx'] = NA
  
  dataset$clinical_data_pathology_T_stage <- as.numeric( factor(
    replaceValues(dataset$clinical_data_pathology_T_stage,
                  findTStage,
                  changeTStage),
    levels = c) )
  
  if ('clinical_data_pathologic_stage' %in% colnames(dataset)) {
    
    dataset$clinical_data_pathologic_stage[dataset$clinical_data_pathologic_stage == 'tx'] = NA
    
    dataset$clinical_data_pathologic_stage <- as.numeric( factor(
      replaceValues(dataset$clinical_data_pathologic_stage, 
                    findPathologicStage, 
                    changePathologicStage),
      levels = d) )
  }
  
  dataset$clinical_data_gender <- as.numeric (factor(
    dataset$clinical_data_gender,
    levels= e) )
  
  #NEW
  EE = c('none', 'minimal (t3)', 'moderate/advanced (t4a)', 'very advanced (t4b)')
  
  if ('clinical_data_extrathyroidal_extension' %in% colnames(dataset)) {
    
    dataset$clinical_data_extrathyroidal_extension <- as.numeric (factor(
      dataset$clinical_data_extrathyroidal_extension,
      levels= EE) )
  }
  
  
  if ('clinical_data_year_of_tobacco_smoking_onset' %in% colnames(dataset)) {
    
    dataset$clinical_data_year_of_tobacco_smoking_onset = sapply(dataset$clinical_data_year_of_tobacco_smoking_onset, function (x) round(x/10)*10)
  
  }
  
  #NEWER
  
  survivalLevels = c(FALSE, TRUE)
  
  if ('survival' %in% colnames(dataset)) {
    
    dataset$survival = as.numeric (factor(dataset$survival,
                                             levels = survivalLevels) )
    
  }
  
  n <- max(length(a), length(b), length(c), length(d), length(e), length(EE))
  
  length(a) <- n
  length(b) <- n
  length(c) <- n
  length(d) <- n
  length(e) <- n
  length(EE) <- n
  
  
  num <- data.frame( id=c(1:n), 
                     clinical_data_pathology_M_stage=a,
                     clinical_data_pathology_N_stage=b,
                     clinical_data_pathology_T_stage=c,
                     clinical_data_pathologic_stage=d,
                     clinical_data_gender=e,
                     EE=EE)
  
  dir.create("./Map")
  write.csv(num, sprintf("./Map/Map %s.csv",name) )
  
  assign("Map", num, envir = .GlobalEnv)
  
  return(dataset)
}

biserialCor <- function(MEs, y, name = 'Gender') {
  n = length(MEs)
  
  moduleTraitCor = data.frame()
  moduleTraitPvalue = data.frame()
  
  for (i in 1:n) {
    x = MEs[, i]
    MEname = names(MEs)[i] 
    
    correlation = cor.test(x, y)
    
    moduleTraitCor = rbind(moduleTraitCor, correlation$estimate)
    moduleTraitPvalue = rbind(moduleTraitPvalue, correlation$p.value)

  }
  
  moduleTraitCor = data.matrix(moduleTraitCor)
  moduleTraitPvalue = data.matrix(moduleTraitPvalue)
  
  rownames(moduleTraitCor) = names(MEs)
  rownames(moduleTraitPvalue) = names(MEs)
  
  colnames(moduleTraitCor) = name
  colnames(moduleTraitPvalue) = name
  
  return(list(moduleTraitCor, moduleTraitPvalue))
}

correlateTraits <- function(MEs, ClinicalData, name) {
  MEs = orderMEs(MEs)
  
  ContinousTraits = c('clinical_data_years_to_birth','clinical_data_number_of_lymph_nodes',
                      'clinical_data_days_to_death', 'clinical_data_tumor_size', 'clinical_data_karnofsky_performance_score',
                      'clinical_data_psa_value')
  OrdinalTraits = c('clinical_data_pathologic_stage', 'clinical_data_pathology_T_stage', 'clinical_data_pathology_M_stage',
                    'clinical_data_pathology_N_stage', 'clinical_data_year_of_tobacco_smoking_onset', 'clinical_data_gleason_score',
                    'clinical_data_extrathyroidal_extension', 'ips_ctla4_neg_pd1_neg', 'ips_ctla4_neg_pd1_pos',
                    'ips_ctla4_pos_pd1_neg', 'ips_ctla4_pos_pd1_pos')
  BinaryTraits = c('clinical_data_gender')
  
  #Remove NA values
  ClinicalData = ClinicalData[, c(ContinousTraits, OrdinalTraits, BinaryTraits)]
  Traits = RemoveMissingTraits(ClinicalData, c(ContinousTraits, OrdinalTraits, BinaryTraits))
  ClinicalData = nominal2numeric2(ClinicalData[, Traits], name)
  
  #Format traits
  ContinousTraits = ContinousTraits[ContinousTraits %in% Traits]
  OrdinalTraits = OrdinalTraits[OrdinalTraits %in% Traits]
  BinaryTraits = BinaryTraits[BinaryTraits %in% Traits]
  
  samples = dim(MEs)[1]
  names = rownames(ClinicalData) 
  MEnames = names(MEs)
  
  #Continous
  ContinousData = ClinicalData[, ContinousTraits]
  ConmoduleTraitCor = cor(MEs, ContinousData, use = "p", method = 'p')
  ConmoduleTraitPvalue = corPvalueStudent(ConmoduleTraitCor, samples)
  PlotHeatmap(ConmoduleTraitCor, ConmoduleTraitPvalue, names(ContinousData), MEnames, name, 'Continous Traits')
  
  #Ordinal
  OrdinalData = ClinicalData[, OrdinalTraits]
  OrdmoduleTraitCor = cor(MEs, OrdinalData, use = "p", method = 's')
  OrdmoduleTraitPvalue = corPvalueStudent(OrdmoduleTraitCor, samples)
  PlotHeatmap(OrdmoduleTraitCor, OrdmoduleTraitPvalue, names(OrdinalData), MEnames, name, 'Ordinal Traits')
  
  #Binary
  #BinmoduleTraitCor = data.frame(row.names = 1:length(MEs))
  #BinmoduleTraitPvalue = data.frame(row.names = 1:length(MEs))
  
  BinaryData = ClinicalData[, BinaryTraits]
  
  #for (i in 1:length(BinaryTraits)) {
  #  
  #  results = biserialCor(MEs, BinaryData[, BinaryTraits[i]], name = BinaryTraits[i])
  #  BinmoduleTraitCor[, BinaryTraits[i]] = results[[1]]
  #  BinmoduleTraitPvalue[, BinaryTraits[i]] = results[[2]]
  #}
  
  #rownames(BinModuleTraitCor) = rownames(results[[1]])
  #rownames(BinModuleTraitPvalue) = rownames(results[[2]])
  
  #BinmoduleTraitCor = as.matrix(BinmoduleTraitCor)
  #BinmoduleTraitPvalue = as.matrix(BinmoduleTraitPvalue)
  
  results = biserialCor(MEs, BinaryData)
  BinmoduleTraitCor = results[[1]]
  BinmoduleTraitPvalue = results[[2]]
  
  PlotHeatmap(BinmoduleTraitCor, BinmoduleTraitPvalue, BinaryTraits, MEnames, name, 'Binary Traits')
  
}

PlotHeatmap <- function(moduleTraitCor, moduleTraitPvalue, names, MEnames, tissueName, corType) {
  sizeGrWindow(50, 50)
  
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  
  dev.new()
  
  pdf(sprintf('./Correlation Heatmap/Correlation of %s - %s.pdf', tissueName, corType), height = 15, width = 15)
  
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names,
                 yLabels = MEnames,
                 ySymbols = MEnames,
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 verticalSeparator.x = 0.25,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  
  dev.off()
}

testFunction <- function() {
  library('RSQLite')
  
  con <- dbConnect( RSQLite::SQLite(), dbname="./cancer")
  query <- NULL
  
  query <- sprintf("SELECT * FROM patient")
  clinicalData = dbGetQuery( con, query)
  
  file = 'LUSC_Tumor'
  name = file
  
  file = sprintf('./Data/Module Eigengenes/%s.RDATA', file)
  load(file)
  
  ids = clinicalData$barcode %in% names
  clinicalData = clinicalData[ids, ]
  
  rownames(clinicalData) = clinicalData$barcode
  
  #NEWER
  
  #if ('clinical_data_days_to_death' %in% colnames(clinicalData)) {
    
  #  clinicalData$survival = is.na(clinicalData$clinical_data_days_to_death)
    
  #}
  
  correlateTraits(MEs, clinicalData, name)
  
  return(clinicalData)
}

#Change data frame with these values
#Turn from nominal to numeric

findPathologicStage <- c('stage iii',
                         'stage iv',
                         'stage ii',
                         'stage i',
                         'stage ia',
                         'stage iia',
                         'stage iib',
                         'stage iiia',
                         'stage ib',
                         'stage iiic',
                         'stage iiib',
                         'stage x',
                         'stage ivb',
                         'stage iva',
                         'stage iic',
                         'stage ivc',
                         'stage 0',
                         'i/ii nos',
                         'Stage III')

changePathologicStage <- c('iii',
                           'iv',
                           'ii',
                           'i',
                           'i',
                           'ii',
                           'ii',
                           'iii',
                           'i',
                           'iii',
                           'iii',
                           'x',
                           'iv',
                           'iv',
                           'ii',
                           'iv',
                           '0',
                           'i', #i/ii nos
                           'iii')

findMStage <- c('m0',
                'mx',
                'm1',
                'cm0 (i+)',
                'm1b',
                'm1a',
                'm1c')

changeMStage <- c('m0',
                  'mx',
                  'm1',
                  'm0', #cm0
                  'm1',
                  'm1',
                  'm1')

findNStage <- c('n0',
                'n1',
                'n2',
                'nx',
                'n3',
                'n0 (i-)',
                'n1b',
                'n1mi',
                'n1a',
                'n0 (mol+)',
                'n0 (i+)',
                'n2a',
                'n3a',
                'n1c',
                'n3b',
                'n3c',
                'n2b',
                'n2c')

changeNStage <- c('n0',
                  'n1',
                  'n2',
                  'nx',
                  'n3',
                  'n0', #n0 i
                  'n1',
                  'n1', #n1mi
                  'n1',
                  'n0', #n0 i
                  'n0', #n0 i
                  'n2',
                  'n3',
                  'n1',
                  'n3',
                  'n3',
                  'n2',
                  'n2')

findTStage <- c('t3',
                't4a',
                't3b',
                't2a',
                't2b',
                't3a',
                't2',
                't4',
                't4b',
                'tx',
                't0',
                't1',
                't1c',
                't1b',
                't1a',
                't4d',
                't1b2',
                't1b1',
                't2a2',
                't2a1',
                't1a1',
                'tis',
                't3c',
                't2c')

changeTStage <- c('t3',
                  't4',
                  't3',
                  't2',
                  't2',
                  't3',
                  't2',
                  't4',
                  't4',
                  'tx',
                  't0',
                  't1',
                  't1',
                  't1',
                  't1',
                  't4',
                  't1',
                  't1',
                  't2',
                  't2',
                  't1',
                  'tis',
                  't3',
                  't2')