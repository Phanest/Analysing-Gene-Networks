library('RSQLite')
library(ggplot2)
library(dplyr)
library(tidyr)

CreateBarplot <- function(table, column, path = '.', na.rm = FALSE, percent = TRUE)
{
  name <- gsub( '.*\\$', '', deparse(substitute( column ) ) )
  
  if (length(column) == 0)
  {
    na.rm = TRUE
  }
  
  if (!na.rm)
  {
    data <- table %>% group_by( table[[name]] ) %>%
      tally() %>% complete(table[[name]], fill = list(n = 0))
  }
  else
  {
    data <- table %>% group_by( table[[name]] ) %>%
      tally()
  }
  
  if (percent)
  {
    data <- data %>% mutate( percentage = n/sum(n) * 100)
    
    plot <- ggplot(data, 
                   aes(data[[ 'table[[name]]' ]],
                       percentage, 
                       fill= data[[ 'table[[name]]' ]],
                       label= round( percentage, digits = 2) ) ) +
      geom_bar(stat= 'identity',position = 'dodge') + theme_bw() + 
      xlab(name) + guides(fill= guide_legend(title=name ) ) +
      geom_text(size = 3, position = position_stack(vjust = 0.5))
  }
  else
  {
    plot <- ggplot(data, 
                   aes(data[[ 'table[[name]]' ]],
                       n, 
                       fill= data[[ 'table[[name]]' ]],
                       label= round( n, digits = 2) ) ) +
      geom_bar(stat= 'identity',position = 'dodge') + theme_bw() + 
      xlab(name) + guides(fill= guide_legend(title=name ) ) +
      geom_text(size = 3, position = position_stack(vjust = 0.5))
  }
  
  ggsave(sprintf('%s.eps', name), plot = plot, device = 'eps' )
  ggsave(sprintf('%sIMG.jpeg', name), plot = plot, device = 'jpeg' )
  
}

replaceColumns <- function(table, column, find, change)
{
  
  name <- gsub( '.*\\$', '', deparse(substitute( column ) ) )
  
  for (i in 1:length(find))
  {
    table[[name]][table[[name]] %in% find[i]] <- change[i]
  }
  
  return (table)
}

histogram <- function(column)
{
  n <- 1
  
  if (length(column) == 0)
  {
    n <- 0
  }
  
  name <- gsub( '.*\\$', '', deparse(substitute( column ) ) )
  num <- ModeOccurances(column)
  
  jpeg(sprintf('%sIMG.jpeg', name))
  
  if (n == 0)
  {
    plot.new()
  }
  else
  {
    plot <- hist(column, 
                 col='grey', 
                 labels = TRUE, main = '', 
                 xlab = name)
  }
  
  graphics.off()
  
  setEPS()
  postscript(sprintf('%s.eps', name))
  
  if (n == 0)
  {
    plot.new()
  }
  else
  {
    plot <- hist(column, 
                 col='grey', 
                 labels = TRUE, main = '', 
                 xlab = name)
  }
  
  graphics.off()
}

ModeOccurances <- function(column)
{
  return((max(table(column )  ) ) )
}

numberOfTissue <- function(name, normal, cancer)
{
  table = data.frame(c('normal', 'cancer'), c(normal, cancer) )
  
  colnames(table) = c('tissue', 'sample')
  
  plot = ggplot(table,
                aes(tissue,
                    sample,
                    fill = tissue,
                    label = sample
                )
  ) + geom_bar(stat= 'identity', position = 'dodge') + theme_bw() + geom_text(size = 3, position = position_stack( vjust = 0.5 ))
  
  ggsave(sprintf('TissueSamples%s.eps', name), plot = plot, device = 'eps' )
  ggsave(sprintf('TissueSamples%sIMG.jpeg', name), plot = plot, device = 'jpeg' )
  
}

automate <- function(filter, path = '.')
{
  
  for (i in 1:length(filter) )
  {
    setwd('./Images')
    
    dir.create(file.path(getwd(), filter[i] ) )
    
    setwd(file.path(getwd(), filter[i] ) )
    
    data(filter[i])
    
    graphics.off()
  }
  
  setwd('./Images')
}

data <- function(name)
{
  filter = 'clinical_data_pathology_T_stage="tis"'

  if (name == 'All')
  {
    statement = sprintf("SELECT * FROM patient WHERE %s", filter)
  }
  else
  {
    statement = sprintf("SELECT * FROM patient WHERE disease='%s' and %s", name, filter)
  }
  
  table = dbGetQuery( con, statement)
  
  if (name == 'All')
  {
    #Tissue samples
    normal = sprintf("SELECT count(*) FROM (normal JOIN patient) WHERE %s", filter)
    cancer = sprintf("SELECT count(*) FROM (cancer JOIN patient) WHERE %s", filter)
  }
  else
  {
    #Tissue samples
    normal = sprintf("SELECT count(*) FROM (normal JOIN patient ON patient.barcode = SUBSTR(normal.barcode, 0, 13)) WHERE patient.disease='%s' and %s", name, filter)
    cancer = sprintf("SELECT count(*) FROM (cancer JOIN patient ON patient.barcode = SUBSTR(cancer.barcode, 0, 13)) WHERE patient.disease='%s' and %s", name, filter)
  }
  
  normal = dbGetQuery( con, normal)
  cancer = dbGetQuery( con, cancer)
  
  numberOfTissue(name, normal[[1]], cancer[[1]])
  
  #Histograms for: clinical_data_date_of_initial_pathologic_diagnosis,
  #ips_ctla4_neg_pd1_neg, ips_ctla4_neg_pd1_pos, ips_ctla4_pos_pd1_pos
  #ips_ctla4_pos_pd1_neg
  createHistograms(table)
  #Barplots for: clinical_data_gender, clinical_data_pathologic_stage
  #clinical_data_pathology_N_stage, clinical_data_pathology_M_stage,
  #clinical_data_pathology_T_stage
  Barplots(table)
}

createHistograms <- function(table)
{
  histogram(table$clinical_data_date_of_initial_pathologic_diagnosis)
  histogram(table$ips_ctla4_neg_pd1_neg)
  histogram(table$ips_ctla4_neg_pd1_pos)
  histogram(table$ips_ctla4_pos_pd1_pos)
  histogram(table$ips_ctla4_pos_pd1_neg)
}

Barplots <- function(table)
{
  table <- replaceColumns(table, 
                          table$clinical_data_pathologic_stage, 
                          findPathologicStage,
                          changePathologicStage)
  
  table <- replaceColumns(table,
                          table$clinical_data_pathology_N_stage,
                          findNStage,
                          changeNStage
  )
  
  table <- replaceColumns(table,
                          table$clinical_data_pathology_M_stage,
                          findMStage,
                          changeMStage)
  
  table <- replaceColumns(table,
                          table$clinical_data_pathology_T_stage,
                          findTStage,
                          changeTStage)
  
  CreateBarplot(table, table$clinical_data_gender)
  CreateBarplot(table, table$clinical_data_pathologic_stage)
  CreateBarplot(table, table$clinical_data_pathology_N_stage)
  CreateBarplot(table, table$clinical_data_pathology_M_stage)
  CreateBarplot(table, table$clinical_data_pathology_T_stage)
}

con <- dbConnect( RSQLite::SQLite(), dbname="./cancer")

filter <- c('All',
            'BLCA',
            'BRCA',
            'CESC',
            'COAD',
            'GBM',
            'HNSC',
            'KICH',
            'KIRC',
            'KIRP',
            'LIHC',
            'LUAD',
            'LUSC',
            'OV',
            'PAAD',
            'PRAD',
            'READ',
            'SKCM',
            'STAD',
            'THCA',
            'UCEC')

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
                           'i/ii nos',
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
                  'cm0 (i+)',
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
                  'n0 (i-)',
                  'n1',
                  'n1mi',
                  'n1',
                  'n0 (mol+)',
                  'n0 (i+)',
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