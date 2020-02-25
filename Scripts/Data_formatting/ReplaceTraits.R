"For some applications it may be useful to replace specific traits with simpler values.
This script replaces some trait values with more appropriate ones (t1c -> t1)"

#Replaces:
#T Stage
#N Stage
#M Stage
#Pathologic Stage

replaceTraits <- function(dataset)
{
  
  a <- unique(changeMStage)
  b <- unique(changeNStage)
  c <- unique(changeTStage)
  d <- unique(changePathologicStage)
  
  dataset$clinical_data_pathology_M_stage <- factor(
    replaceValues(dataset$clinical_data_pathology_M_stage, 
                  findMStage, 
                  changeMStage),
    levels= a)
  
  dataset$clinical_data_pathology_N_stage <- factor(
    replaceValues(dataset$clinical_data_pathology_N_stage, 
                  findNStage, 
                  changeNStage),
    levels = b)
  
  dataset$clinical_data_pathology_T_stage <- factor(
    replaceValues(dataset$clinical_data_pathology_T_stage,
                  findTStage,
                  changeTStage),
    levels = c)
  
  dataset$clinical_data_pathologic_stage <- factor(
    replaceValues(dataset$clinical_data_pathologic_stage, 
                  findPathologicStage, 
                  changePathologicStage),
    levels = d)

  #NEW
  
  PAM50 = unique(dataset$clinical_data_PAM50MRNA)[-1]
  NRAS = unique(dataset$clinical_data_NRAS)[-1]
  NF1 = unique(dataset$clinical_data_NF1)[c(-1,-9)]
  MSI = unique(dataset$clinical_data_msi_status)[-1]
  EE = unique(dataset$clinical_data_extrathyroidal_extension)[-1]
  
  dataset$clinical_data_clinical_stage <- factor(
    dataset$clinical_data_clinical_stage,
    levels= d)
  
  dataset$clinical_data_PAM50MRNA <- factor(
    dataset$clinical_data_PAM50MRNA,
    levels= PAM50)
  
  dataset$clinical_data_NRAS <- factor(
    dataset$clinical_data_NRAS,
    levels= NRAS)
  
  dataset$clinical_data_NF1 <- factor(
    dataset$clinical_data_NF1,
    levels= NF1)
  
  dataset$clinical_data_msi_status <- factor(
    dataset$clinical_data_msi_status,
    levels= MSI)
  
  #TODO
  dataset$clinical_data_extrathyroidal_extension <- factor(
    dataset$clinical_data_extrathyroidal_extension,
    levels= EE)
  
  return(dataset)
}

replaceValues <- function(column, findValues, changeValues)
{
  
  for(i in 1:length( findValues ) )
  {
    pos = which( column == findValues[i] )
    column[pos] = changeValues[i]
  }
  
  return(column)
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
                  'n0 i',
                  'n1',
                  'n1mi',
                  'n1',
                  'n0 i',
                  'n0 i',
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