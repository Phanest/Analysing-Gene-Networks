"Creates plots from the metric data and compares the Normal vs Tumor tissue"

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggplot2)
library(extrafont)
library(WGCNA)

plotBox <- function(data, name, logV = FALSE, save = FALSE) {
  
  ylab = 'value'
  if (logV) ylab = 'log2 value'
  
  if (save) {}
  sizeGrWindow(15,15)
  
  dev.new()
  
  tiff(file= sprintf("./Module Metrics/Boxplot Metrics - %s.tiff", name),
       width = 5, height = 5, units = 'in',
       res = 500, compression = 'lzw')
    
  data %>%
   ggplot( aes(x = name, y = value, fill = name)) +
    geom_boxplot(width = 0.5) +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter( aes(color = color, size = 0.2, alpha = 0.9)) +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(name) +
    xlab("") +
    ylab(ylab)
  
}

plotViolin <- function(data, name, logV = FALSE, save = FALSE) {
  
  ylab = 'value'
  if (logV) ylab = 'log2 value'
  
  if (save) {}
  sizeGrWindow(15,15)
  
  dev.new()
  
  tiff(file= sprintf("./Module Metrics/Violin Metrics - %s.tiff", name),
       width = 5, height = 5, units = 'in',
       res = 500, compression = 'lzw')
  
  data %>%
    ggplot( aes(x = name, y = value, fill = name)) +
    geom_violin() +
    ggtitle(name) +
    ylab(ylab)
  
}

#Read Normal metrics
Normal = read.csv('./Data/Metrics LUSC_Normal.csv')
rownames(Normal) = Normal$module
names = colnames(Normal)

#Read Tumor metrics
Tumor = read.csv('./Data/Metrics LUSC_Tumor.csv')
rownames(Tumor) = Tumor$module

meanValues = data.frame(preservation = 0, trait = 0, valueTumor = 0, valueNormal = 0)

#Filter Preserved
PreservedNormal = Normal[Normal[,10] < 2, ]
PreservedTumor = Tumor[Tumor[, 10] < 2, ]

for (i in 2:9) {
  name = names[i]
  
  data = data.frame(name = rep(x = 'Tumor', length(PreservedTumor[, i])),
                    color = PreservedTumor[, 1], 
                    value = PreservedTumor[, i])
  
  data = rbind(data,
               data.frame(
                 name = rep(x = 'Normal', length(PreservedNormal[, i]) ),
                 color = PreservedNormal[, 1],
                 value = PreservedNormal[, i]
               ))
  
  # plotBox(data, name)
  # plotViolin(data, name)
  logV = FALSE

  # #Log2 + 1
  # data$value = log2(data$value + 1)
  # logV = TRUE

  # #divided by max
  # maxA = max(PreservedTumor[, i])
  # maxB = max(PreservedNormal[, i])
  #
  # data[data$name == 'Tumor', 3] = data[data$name == 'Tumor', 3]/maxA
  # data[data$name == 'Normal', 3] = data[data$name == 'Normal', 3]/maxB

  plotBox(data, name, logV, save = TRUE)

  plotViolin(data, name, logV, save = TRUE)
  
  meanValues = rbind(meanValues, 
                     data.frame(
                       preservation = 'Preserved',
                       trait = name,
                       valueTumor = mean(PreservedTumor[, i]),
                       valueNormal = mean(PreservedNormal[, i])
                     ))
}