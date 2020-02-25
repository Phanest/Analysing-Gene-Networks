"This script removes traits with less than 15% non-NA values"

selection= "clinical_data_gender,clinical_data_clinical_stage,clinical_data_pathologic_stage,clinical_data_pathology_M_stage,
clinical_data_pathology_N_stage,clinical_data_pathology_T_stage,clinical_data_years_to_birth
,ips_ctla4_neg_pd1_neg,ips_ctla4_neg_pd1_pos,ips_ctla4_pos_pd1_neg,ips_ctla4_pos_pd1_pos,clinical_data_days_to_death,clinical_data_number_of_lymph_nodes,clinical_data_year_of_tobacco_smoking_onset,clinical_data_tumor_size,
clinical_data_progression_free,clinical_data_PAM50MRNA,clinical_data_NRAS,clinical_data_NF1,clinical_data_msi_status,
clinical_data_karnofsky_performance_score,clinical_data_HER2,clinical_data_extrathyroidal_extension,clinical_data_psa_value,
clinical_data_gleason_score"

RemoveMissingTraits <- function(clinicalData, attributes) {
  
  clinicalData = clinicalData[, which(names(clinicalData) %in% attributes)];
  
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
  return(names(clinicalData))
}
