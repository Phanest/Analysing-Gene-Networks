"Does a log transformation of the data and saves the new data"

library('RSQLite')
library(arules)

src = 'Data/Transposed Gene Expression'

file = list.files(src)[4]

src = paste(src, file, sep='/')

Raw_Data = read.delim(src);


Data = as.data.frame(Raw_Data[, -c(1)]); #take all columns except the first with all the names
names(Data) = names(Raw_Data)[-c(1)];
rownames(Data) = strtrim(Raw_Data$gene_id , 12);

#Data = as.data.frame(Data)

LogData = log(Data + 1)
hist(as.numeric(LogData[1,]))

#Row names are lost
LogData = do.call(data.frame,lapply(LogData, function(x) replace(x, is.infinite(x), -1000000000)))
rownames(LogData) = rownames(Data)

########################Format data for SampleNetworks

con <- dbConnect( RSQLite::SQLite(), dbname="cancer")

query <- sprintf("SELECT * FROM patient")
allTraits = dbGetQuery( con, query)

#sapply(data, is.factor)
selection= "datasource,clinical_data_gender,clinical_data_clinical_stage,clinical_data_pathologic_stage,clinical_data_pathology_M_stage,
clinical_data_pathology_N_stage,clinical_data_pathology_T_stage,clinical_data_years_to_birth
,ips_ctla4_neg_pd1_neg,ips_ctla4_neg_pd1_pos,ips_ctla4_pos_pd1_neg,ips_ctla4_pos_pd1_pos,clinical_data_days_to_death,clinical_data_number_of_lymph_nodes,clinical_data_year_of_tobacco_smoking_onset,clinical_data_tumor_size,
clinical_data_progression_free,clinical_data_PAM50MRNA,clinical_data_NRAS,clinical_data_NF1,clinical_data_msi_status,
clinical_data_karnofsky_performance_score,clinical_data_HER2,clinical_data_extrathyroidal_extension,clinical_data_psa_value,
clinical_data_gleason_score"

factors = c(2,3,4,5,6,7,8,19,20,21,22,24,25)
closerLook = c(3:8,24)

answer=readline(prompt="Do you want to include 'clinical_data_years_to_birth' values (y/n): ")
while(length(intersect(c("y","n"),answer))==0){
  answer=readline(prompt="Response (y/n): ")
}

if (answer == 'y') {
  selection = c(selection, 'clinical_data_years_to_birth')
  factors = c(factors, 26)
}

barcodes = allTraits$barcode;

attributes <- unlist(strsplit(selection, ","))
attributes <- lapply(attributes, trimws)
attributes <- c(attributes)
attributes <- unlist(attributes)

allTraits = allTraits[, c(1, which(names(allTraits) %in% attributes))];

#patient names
Samples = strtrim( rownames(LogData), 12);

traitRows = barcodes %in% Samples

datTraits = allTraits[traitRows, ]
rownames(datTraits) = allTraits[traitRows, 1]

#####################################################
#Process clinical age values

if (answer == 'y') {
  missingValues = sum(is.na(datTraits$clinical_data_years_to_birth))
  sprintf('There are %s missing "clinical_data_years_to_birth" values.
          Do you want to exclude them or use the mean? (exclude/mean)', missingValues)
  
  answer=readline(prompt="Response (exclude/mean): ")
  while(length(intersect(c("exclude","mean"),answer))==0){
    answer=readline(prompt="Response (exclude/mean): ")
  }
  
  if (answer == 'exclude') {
    datTraits = datTraits[complete.cases(allTraits$clinical_data_years_to_birth), ]
  } else {
    mean = mean(datTraits$clinical_data_years_to_birth, na.rm = T)
    datTraits$clinical_data_years_to_birth[is.na(datTraits$clinical_data_years_to_birth)] = mean 
  }
  
  #Add a new clinical age column and discretize it
  datTraits$color_clinical_data_years_to_birth = discretize(datTraits$clinical_data_years_to_birth)
  attributes = c(attributes, 'color_clinical_data_years_to_birth')
  
  barcodes = datTraits$barcode;
}

#####################################################

colLen = length(names(datTraits))

#Order columns
datTraits = datTraits[c('barcode', attributes)]

#Replace Traits
source('ReplaceTraits.R')
datTraits = replaceTraits(datTraits)
#Replace Traits

#Remove any column with a NA value
ikeep = c()
inames = names(datTraits)
for (i in 1:colLen) {
  if (any(is.na(datTraits[, i]))) {
    ikeep = c(ikeep, i)
    
    #Remove column from factors
    j = which(factors==i)
    if(sum(j) != 0) {
      factors = factors[-j]
    }
    
    #Remove column from closerLook
    j = which(closerLook==i)
    if(sum(j) != 0) {
      closerLook = closerLook[-j]
    }
  }
}
datTraits = subset(datTraits, select=-ikeep)
factors = inames[factors]
closerLook = inames[closerLook]

print('The following columns are factors:')
print(factors)
print('Multivariate analysis on:')
print(closerLook)

#Find new indexes for factors and closerLook
inames = names(datTraits)
factors = match(factors, inames)
closerLook = match(closerLook, inames)

#We need this for the SampleNetwork "grouplabel" attribute
datTraits$group = 'all'

exprRows = Samples %in% barcodes
LogData = LogData[exprRows, ]

#rows must be in the same order
LogData = LogData[match( datTraits$barcode, rownames(LogData)) , ]

#For SampleNetwork, rows should be genes and columns sample names
LogData = t(LogData)
LogData = as.data.frame(LogData)

#For testing purposes
LogData$Gene = rownames(LogData)
LogData = LogData[, c('Gene', names(LogData)[-length(names(LogData))])]

######################Sample Network

source('SampleNetwork.R')

#LogData has negative values

colLen = length(names(LogData))
traitNum = length(names(datTraits))

#btrait1: attributes used for significance testing
#trait1: attributes whose individual values should be used for significance testing
SampleNetwork(datExprT=LogData,
              method1="correlation",
              impute1=FALSE,
              subset1=NULL,
              skip1=1,
              indices1=list(c(2:colLen)),
              sampleinfo1=datTraits,
              samplelabels1=1,
              grouplabels1=traitNum,
              fitmodels1=TRUE,
              whichmodel1="multivariate",
              whichfit1="pc1",
              btrait1=c(2:(traitNum-1)),
              trait1=closerLook,
              asfactors1=factors,
              projectname1= file,
              cexlabels1=0.7,
              normalize1=TRUE,
              replacenegs1=FALSE,
              exportfigures1=TRUE,
              verbose=TRUE)

dst = 'Data/Log Transformed Transformed Gene Expression - No Outliers'
dst = paste(dst, file, sep='/')