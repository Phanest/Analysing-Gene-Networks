library('RSQLite')

source('./Scripts/DataGatherer.R')

con <- dbConnect( RSQLite::SQLite(), dbname="./cancer")
query <- NULL

query <- sprintf("SELECT * FROM patient")
clinicalData = dbGetQuery( con, query)

clinicalData = nominal2numeric2(clinicalData)

src = './Data/Preprocessed Gene Data Final'

geneFile = "BRCA_Tumor_scaled_estimate_tpm.txt_all_1081_Qnorm.csv"

dest = './'

clusters(src, geneFile, dest, clinicalData = clinicalData)