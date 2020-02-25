"Hubgene boxplot according to trait value, tests genes for statistical differences between differing clinical traits"

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggplot2)
library(extrafont)
library(WGCNA)
library(car)

#numerical to string
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

plotBox <- function(data, name, fileName, logV = FALSE, save = TRUE) {
  
  ylab = 'expression'
  if (logV) ylab = 'log2 value'
  
  data %>%
    ggplot( aes(x = name, y = value)) +
    geom_boxplot(width = 0.5) +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter( shape = 16, position = position_jitter(0.2)) +
    stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(name) +
    xlab("") +
    ylab(ylab)
  
  
  if (save) ggsave( sprintf('./HubGenes/%s.tiff', fileName), width = 5, height = 5, units = 'in', dpi = 300)
}

plotViolin <- function(data, name, fileName, logV = FALSE, save = TRUE) {
  
  ylab = 'expression'
  if (logV) ylab = 'log2 value'
  
  data %>%
    ggplot( aes(x = name, y = value, fill = name)) +
    geom_violin() +
    ggtitle(name) +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(name) +
    xlab("") +
    ylab(ylab)
  
  if (save) ggsave( sprintf('./HubGenes/%s Violin.tiff', fileName), width = 5, height = 5, units = 'in', dpi = 300)
}

addRows <- function(Data, clinicalData, geneName = '', traitName = '') {
  values = Data[ , geneName]

  if (traitName %in% colnames(clinicalData)) {
    name = clinicalData[, traitName]
    
    map = read.csv('./Map/Map hubGene plot.csv')
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
  
  name = as.character(name)
  
  data = data.frame(name = name, value = values)

  return(data)
}


testHypothesis <- function(data) {
  testType = 0
  result = 0
  pval = 0
  comment = 0
  
  #ANOVA
  
  #Check for varience (there should be no relation)
  
  #The levane test, tests the homogeneity of the variances, if p > significance it means
  #there isn't a big difference between variances
  #TODO
  levene = leveneTest(value ~ name, data = data)
  levene = levene$`Pr(>F)`[1]
  
  if (levene >= 0.05) {
    comment = sprintf('Levane testing finds no difference in variance (pval: %s), ', levene)
    
    res.aov <- aov(value ~ name, data = data)
    
    #Shapiro test for normality, if p > significance it's normal
    #TODO
    aov_residuals <- residuals(object = res.aov )
    normal = shapiro.test(x = aov_residuals )
    
    if (normal$p.value >= 0.05) {
      comment = sprintf('%sThe data has a normal distribution (shapiro test pval: %s), ', comment, normal$p.value)
      
      Aresult = summary(res.aov)[[1]]
      pval = Aresult$`Pr(>F)`[1]
      
      testType = 'ANOVA'
      
      if (pval < 0.05) {
        comment = sprintf('%sAnova test passed (pval: %s), there is a difference between the groups. ', comment, pval)
        
        result = 'Passed'
        
        Tukey = TukeyHSD(res.aov)
        Tukey = data.frame(Tukey$name)
        Tukey$p.adj = Tukey$p.adj < 0.05
        
        Tukey = paste(rownames(Tukey)[Tukey$p.adj], collapse = ' ')
        comment = sprintf('%sThese groups have a significant difference between them (Tukey test) -> %s. ', comment, Tukey)
        
      } else {
        result = 'Failed'
        
        comment = sprintf('%sAnova test failed (pval: %s). ', comment, pval)
      }
      
    } else {
      comment = sprintf('%sThe data has a non-normal distribution (shapiro test pval: %s), using kruskal-wallis ', comment, normal$p.value)
      
      KWresult = kruskal.test(value ~ name, data = data)
      pval = KWresult$p.value
      
      testType = 'Kruskal-Wallis'
      
      if (pval < 0.05) {
        comment = sprintf('%sKruskal-Wallis test passed (pval: %s), there is a difference between the groups. ', comment, pval)
        
        result = 'Passed'
      } else {
        result = 'Failed'
        
        comment = sprintf('%sKruskal-Wallis test failed (pval: %s). ', comment, pval)
      }
    }
    
  } else {
    comment = sprintf('Levane testing finds differences in variance (pval: %s), ', levene)
    
    values = unique(data$name)
    normal = TRUE
    
    for (i in values) {
      p.value = shapiro.test(x = data$value[data$name == i])
      p.value = p.value$p.value
      
      if (p.value < 0.05) {
        normal = FALSE
        break()
      }
    }
    
    if (normal) {
      comment = sprintf('%sThe data has a normal distribution (shapiro test), using an ANOVA analysis without assumption of equal variances. ', comment)
      
      res.aov <- oneway.test(value ~ name, data = data)
      pval = res.aov$p.value 
      
      testType = 'Welch'
      
      if (pval < 0.05) {
        comment = sprintf('%sWelch test passed (pval: %s)', comment, pval)
        
        result = 'Passed'
      } else {
        comment = sprintf('%sWelch test failed (pval: %s)', comment, pval)
        
        result = 'Failed'
      }
    } else {
      comment = sprintf('%sThe data has a non-normal distribution (shapiro test), using kruskal-wallis ', comment)
      
      result = kruskal.test(value ~ name, data = data)
      pval = result$p.value
      
      testType = 'Kruskal-Wallis'
      
      if (pval < 0.05) {
        comment = sprintf('%sKruskal-Wallis test passed (pval: %s), there is a difference between the groups. ', comment, pval)
        
        result = 'Passed'
      } else {
        comment = sprintf('%sKruskal-Wallis test failed (pval: %s). ', comment, pval)
        
        result = 'Failed'
      }
    }
  }
  
  return(list(testType = testType, result = result, pval = pval, comment = comment))
}

rmSmallValues <- function(data) {
  comment = ''
  
  values = data %>% 
    group_by(name) %>%
    summarise(rows = length(name))
  
  n = dim(values)[1]
  
  values$rows = values$rows < n + 1
  
  values = values$name[values$rows]
  values = as.character(values)
  
  if (length(values) > 0) {
    comment = sprintf('The following factors are removed %s', values)
  }
  
  for (i in values) {
    data = data[-which(data$name == i), ]
  }
  
  return(list(data = data, comment = comment))
}

#returns row (gene, trait, ANOVA result, p-value, comment)
testAndPlot <- function(data, ndata, geneName, traitName) {
  comment = ''
  result = 0
  nresult = data.frame(testType = NA, result = NA, pval = NA)
  
  dir.create(sprintf('./HubGenes/Sig/%s', geneName))
  
  if ( all( is.na(ndata) ) ) {
    
    plotBox(data, geneName, sprintf('Sig/%s/%s', geneName, traitName))
    plotViolin(data, geneName, sprintf('Sig/%s/%s', geneName, traitName))
    
    data = rmSmallValues(data)
    
    comment = sprintf('%s%s', comment, data$comment)
    data = data$data
    
    result = testHypothesis(data)
    comment = sprintf('%s. %s', comment, result$comment)
    
  } else {
    
    plotBox(data, geneName, sprintf('Sig/%s/%s', geneName, traitName))
    plotViolin(data, geneName, sprintf('Sig/%s/%s', geneName, traitName))
    
    tdata = rbind(data, ndata)
    
    plotBox(tdata, geneName, sprintf('Sig/%s/%s Normal', geneName, traitName))
    plotViolin(tdata, geneName, sprintf('Sig/%s/%s Normal', geneName, traitName))
    
    data = rmSmallValues(data)
    
    comment = sprintf('%s%s', comment, data$comment)
    data = data$data
    
    result = testHypothesis(data)
    comment = sprintf('%s. %s', comment, result$comment)
    #row
    
    tdata = rbind(data, ndata)
    
    nresult = testHypothesis(tdata)
    comment = sprintf('%s. Test with normal data returned %s', comment, nresult$comment)
    
  }
  
  row = data.frame(gene = geneName, trait = traitName, typeOfTest = result$testType, test = result$result,
                   pval = result$pval, typeOfTestWithNormal = nresult$testType, testWithNormal = nresult$result,
                   pvalWithNormal = nresult$pval, comment = comment)
  return(row)
}

automateSigTests <- function(hubGenes, datExpr, clinicalData, ndatExpr, nclinicalData) {
  n = dim(hubGenes)[1]
  log = data.frame(gene = 0, trait = 0, typeOfTest = 0, test = 0, pval = 0,
                   typeOfTestWithNormal = 0, testWithNormal = 0, pvalWithNormal = 0, comment = 0)
  
  for (i in 1:n) {
    geneName = as.character(hubGenes[i, 1])
    traitName = as.character(hubGenes[i, 2])
    
    data = addRows(datExpr, clinicalData, geneName, traitName)
    ndata = addRows(ndatExpr, nclinicalData, geneName, 'Normal')
    
    row = testAndPlot(data, ndata, geneName, traitName)
    
    log = rbind(log, row)
  }
  
  log = log[-1, ]
  write.csv(log, './HubGenes/log.csv')
}

hubGenes = read.csv('./HubGenes/TOM^3/Hub Genes.csv')

#########################

library('RSQLite')
source('./correlateTraits.R')

con <- dbConnect( RSQLite::SQLite(), dbname="./cancer")
query <- NULL

query <- sprintf("SELECT * FROM patient")
clinicalData = dbGetQuery( con, query)


name = 'LUSC_Tumor'
file = list.files('./Data/Clusters RDATA')[4]
file = paste('./Data/Clusters RDATA', file, sep = '/')

# traitName = ''
# geneName = ''

load(file)

#load module
file = list.files('./Data/Modules/LUSC_Tumor')[1]
file = paste('./Data/Modules/LUSC_Tumor', file, sep = '/')
load(file)

ids = clinicalData$barcode %in% names
clinicalData = clinicalData[ids, ]

rownames(clinicalData) = clinicalData$barcode

clinicalData = nominal2numeric2(clinicalData, 'hubGene plot')

Data = datExpr

############
#LOAD NORMAL DATA
############

con <- dbConnect( RSQLite::SQLite(), dbname="./cancer")
query <- NULL

query <- sprintf("SELECT * FROM patient")
clinicalDataNormal = dbGetQuery( con, query)

name = 'LUSC_Normal'
file = list.files('./Data/Clusters RDATA')[3]
file = paste('./Data/Clusters RDATA', file, sep = '/')

load(file)

#load module
file = list.files('./Data/Modules/LUSC_Normal')[1]
file = paste('./Data/Modules/LUSC_Normal', file, sep = '/')
load(file)

ids = clinicalDataNormal$barcode %in% names
clinicalDataNormal = clinicalDataNormal[ids, ]

rownames(clinicalDataNormal) = clinicalDataNormal$barcode

clinicalDataNormal = nominal2numeric2(clinicalDataNormal, 'Normal hubGene plot')

DataNormal = datExpr

#AFTER FORMATING

group_by(data, name) %>%
  summarise(
    count = n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE)
  )

