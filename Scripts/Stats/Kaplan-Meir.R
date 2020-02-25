"Reads from hub gene file and does a survival Kaplan-Meir analysis based on correlated trait"

#call kaplan-Meir
  #find p-val
  #plot to ./Kaplan-Meir

##
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

##
replaceValues2 <- function(column, findValues, changeValues)
{
  newColumn = rep(x = NA, length(column))
  
  for(i in 1:length( findValues ) )
  {
    pos = which( column == findValues[i] )
    newColumn[pos] = changeValues[i]
  }
  
  return(newColumn)
}

addRows <- function(Data, clinicalData, geneName = '', traitName = '') {
  values = Data[ , geneName]
  
  if (traitName %in% colnames(clinicalData)) {
    name = clinicalData[, traitName]
    
    map = read.csv('./Map/Map Kaplan-Meier plot.csv')
    map = map[, -1]
    
    if (traitName %in% colnames(map)) {
      i = which(traitName == colnames(map))
      changeValues = map[, i]
      changeValues = changeValues[!is.na(changeValues)]
      changeValues = as.character(changeValues)
      
      findValues = 1:length(changeValues)
      
      name = replaceValues2(name, findValues, changeValues)
    }
  } else {
    name = rep(x = traitName, length(values))
  }
  
  remove = which(clinicalData$clinical_data_days_to_death & clinicalData$clinical_data_days_to_last_known_alive)
  clinicalData$clinical_data_days_to_last_known_alive[remove] = NA 
    
  #Add survival time and censored data
  survival = clinicalData %>% mutate(survival = coalesce(clinical_data_days_to_death,
                                              clinical_data_days_to_last_followup, 
                                              clinical_data_days_to_last_known_alive)) %>%
                          select(barcode, survival)
  
  censored = clinicalData$clinical_data_days_to_last_followup | clinicalData$clinical_data_days_to_last_known_alive
  censored[is.na(censored)] = FALSE
  censored = !censored
  
  #order survival
  names = clinicalData$barcode
  rownames(survival) = survival$barcode
  survival = survival[names, ]
  
  #Remove NAs
  removed = !is.na(survival$survival)
  
  name = as.character(name)
  
  med = median(values)
  values = cut(values,
                    breaks = c(-Inf, med, Inf),
                    labels = c('Low', 'High'))
  
  if (traitName == 'clinical_data_karnofsky_performance_score') {
    name = as.numeric(name)
    name = cut(name,
                 breaks = c(-Inf, 0, 90, Inf),
                 labels = c('Low', 'Medium', 'High'))
  }
  
  data = data.frame(barcode = names, name = name, expression = values, survival = survival$survival, censored = censored)
  data = data[removed, ]
  
  return(data)
}

rmSmallValues <- function(data) {
  
  values = data %>% 
    group_by(name) %>%
    summarise(rows = length(name))
  
  n = dim(values)[1]
  
  values$rows = values$rows < 15
  
  values = values$name[values$rows]
  values = as.character(values)
  
  if (length(values) > 0) {
    cat(sprintf('The following factors are removed %s\n', values))
  }
  
  for (i in values) {
    data = data[-which(data$name == i), ]
  }
  
  return(data)
}

kaplanMeirAndPlot <- function(data, geneName, traitName) {
  cat(sprintf('%s %s\n', geneName, traitName))
  mdata = data
  
  if (traitName != 'clinical_data_days_to_death') {
    
    data = rmSmallValues(data)
    surv_object <- Surv(time = data$survival, event = data$censored)
    
    fit <- survfit(Surv(time = data$survival, event = data$censored) ~ name, data = data)
    
    labs = names(fit$strata)
    labs = sapply(labs, function(x) strsplit(x, '=')[[1]][2])
    labs = unname(labs)
    
    title = traitName
    
    gg = ggsurvplot(data = data, fit, pval = TRUE,
                legend.title = title, legend.labs = labs)
    
    ggsave( sprintf('./Kaplan-Meier/%s.tiff', traitName), print(gg),
            width = 5, height = 5, units = 'in', dpi = 300)
  }
  
  ###########
  data = mdata
  
  surv_object <- Surv(time = data$survival, event = data$censored)
  
  fit <- survfit(Surv(time = data$survival, event = data$censored) ~ expression, data = data)
  
  labs = names(fit$strata)
  labs = sapply(labs, function(x) strsplit(x, '=')[[1]][2])
  labs = unname(labs)
  
  title = strsplit(geneName, '\\.')[[1]][1]
  
  gg = ggsurvplot(data = data, fit, pval = TRUE,
             legend.title = title, legend.labs = labs)
  
  ggsave( sprintf('./Kaplan-Meier/%s.tiff', geneName), print(gg),
          width = 5, height = 5, units = 'in', dpi = 300)
}

##
library('RSQLite')
source('./correlateTraits.R')

con <- dbConnect( RSQLite::SQLite(), dbname="./cancer")
query <- NULL

query <- sprintf("SELECT * FROM patient")
clinicalData = dbGetQuery( con, query)


name = 'LUSC_Tumor'
file = list.files('./Data/Clusters RDATA')[4]
file = paste('./Data/Clusters RDATA', file, sep = '/')

load(file)

#load module
file = list.files('./Data/Modules/LUSC_Tumor')[1]
file = paste('./Data/Modules/LUSC_Tumor', file, sep = '/')
load(file)

ids = clinicalData$barcode %in% names
clinicalData = clinicalData[ids, ]

rownames(clinicalData) = clinicalData$barcode

clinicalData = nominal2numeric2(clinicalData, 'Kaplan-Meier plot')

Data = datExpr

hubGenes = read.csv('./HubGenes/TOM/Hub Genes TOM.csv')

n = dim(hubGenes)[1]

for (i in 1:n) {
  geneName = as.character(hubGenes[i, 1])
  traitName = as.character(hubGenes[i, 2])
  
  data = addRows(datExpr, clinicalData, geneName, traitName)
  
  kaplanMeirAndPlot(data, geneName, traitName)

}
