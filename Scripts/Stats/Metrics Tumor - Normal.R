'Test significance of network metrics'

library(tidyverse)
library(hrbrthemes)
library(viridis)
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
    geom_jitter( aes(color = color, size = 0.2, alpha = 0.9)) +
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
  
  levene = leveneTest(value ~ name, data = data)
  levene = levene$`Pr(>F)`[1]
  
  if (levene >= 0.05) {
    comment = sprintf('Levane testing finds no difference in variance (pval: %s), ', levene)
    
    res.aov <- aov(value ~ name, data = data)
    
    #Shapiro test for normality, if p > significance it's normal

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
testAndPlot <- function(data, ndata, metric) {
  comment = ''
  
  data = rmSmallValues(data)
  
  comment = sprintf('%s%s', comment, data$comment)
  data = data$data
  
  outliers = boxplot(data$value, plot = FALSE)$out
  
  if (length(outliers) > 0) data = data[-which(data$value %in% outliers), ]
  
  comment = sprintf('%s. Removing %s outliers', comment, length(outliers))
  
  ndata = rmSmallValues(ndata)
  
  comment = sprintf('%s. From ndata %s', comment, ndata$comment)
  ndata = ndata$data
  
  outliers = boxplot(ndata$value, plot = FALSE)$out
  
  if (length(outliers) > 0) ndata = ndata[-which(ndata$value %in% outliers), ]
  
  comment = sprintf('%s. Removing %s outliers', comment, length(outliers))
  
  data = rbind(data, ndata)
  
  result = testHypothesis(data)
  comment = sprintf('%s. %s', comment, result$comment)
  
  row = data.frame(metric = metric, typeOfTest = result$testType, test = result$result,
                   pval = result$pval, comment = comment)
  return(row)
}

Normal = read.csv('./Data/Metrics LUSC_Normal.csv')
rownames(Normal) = Normal$module
names = colnames(Normal)

#Read Tumor metrics
Tumor = read.csv('./Data/Metrics LUSC_Tumor.csv')
rownames(Tumor) = Tumor$module

meanValues = data.frame(preservation = 0, trait = 0, valueTumor = 0, valueNormal = 0)

#Filter Preserved
PreservedNormal = Normal#[Normal[,10] > 10, ]
PreservedTumor = Tumor#[Tumor[, 10] > 10, ]

a = data.frame(metric = 0, typeOfTest = 0, test = 0, pval = 0, comment = 0)

for (i in 2:9) {
  name = names[i]
  
  data = data.frame(name = rep(x = 'Tumor', length(PreservedTumor[, i])),
                    color = PreservedTumor[, 1],
                    value = PreservedTumor[, i])
  
  ndata = data.frame(name = rep(x = 'Normal', length(PreservedNormal[, i]) ),
                     color = PreservedNormal[, 1],
                     value = PreservedNormal[, i])
  
  pdata = rbind(data, ndata)
  
  plotBox(pdata, name, sprintf('Boxplot Metrics %s', name))
  plotViolin(pdata, name, sprintf('Violin Metrics %s', name))
  
  #Change this depending on data
  data$value = log2(data$value + 1)
  ndata$value = log2(ndata$value + 1)
  
  row = testAndPlot(data, ndata, name)
  
  a = rbind(a, row)
}

a = a[-1, ]
