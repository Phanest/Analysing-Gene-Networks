"This script checks which of the given attributes have non-NA values in the clinical data
and returns the clinical data with these attributes"

selection= "clinical_data_gender,clinical_data_clinical_stage,clinical_data_pathologic_stage,clinical_data_pathology_M_stage,
clinical_data_pathology_N_stage,clinical_data_pathology_T_stage,clinical_data_years_to_birth
,ips_ctla4_neg_pd1_neg,ips_ctla4_neg_pd1_pos,ips_ctla4_pos_pd1_neg,ips_ctla4_pos_pd1_pos,clinical_data_days_to_death,clinical_data_number_of_lymph_nodes,clinical_data_year_of_tobacco_smoking_onset,clinical_data_tumor_size,
clinical_data_progression_free,clinical_data_PAM50MRNA,clinical_data_NRAS,clinical_data_NF1,clinical_data_msi_status,
clinical_data_karnofsky_performance_score,clinical_data_HER2,clinical_data_extrathyroidal_extension,clinical_data_psa_value,
clinical_data_gleason_score"

CheckForValues <- function(clinicalData, attributes=selection) {
  
  #Take only values in attributes
  attributes <- unlist(strsplit(attributes, ","))
  attributes <- lapply(attributes, trimws)
  attributes <- c(attributes)
  attributes <- unlist(attributes)
  clinicalData = clinicalData[, c(1, which(names(clinicalData) %in% attributes))];
  
  samples = dim(clinicalData)[1]
  toBeRemoved = c()
  
  print('attributes: number of values')
  ColVals = names(clinicalData)
  for (val in ColVals) {
    valid = sum(!is.na(clinicalData[, val]))
    cat(sprintf('%s: %d\n', val, valid))
    
    if (valid/samples < 0.15) {
      cat(sprintf('Removing column %s\n', val))
      toBeRemoved = c(toBeRemoved, val)
    }
  }
  
  
  pos = which(ColVals %in% toBeRemoved)
  clinicalData = clinicalData[, -pos]
  return(clinicalData)
}

library('RSQLite')

con <- dbConnect( RSQLite::SQLite(), dbname="./cancer")
query <- NULL

query <- sprintf("SELECT * FROM patient")
clinicalData = dbGetQuery( con, query)

#Find file path
src = './Data/Clusters RDATA'
file = list.files(src)[8]

print(file)
file = paste(src, file, sep = '/')

#Load file
load(file)

#Get patient names
Ids = rownames(datExpr)
Ids = strtrim(Ids, 12)

clinicalData$barcode = strtrim(clinicalData$barcode, 12)

clinicalData = clinicalData[which(clinicalData$barcode %in% Ids),] #TODO
rownames(clinicalData) = clinicalData$barcode

#Check if all data comes from the same file
if (length(unique(clinicalData$disease)) > 1) {
  print('Data is not clean. Patients from multiple diseases.')
  stop('Program was stopped.')
}

CheckForValues(clinicalData)
