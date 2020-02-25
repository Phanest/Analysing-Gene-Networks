library('RSQLite')

source('./DataGatherer.R')
source('./ConstructNetwork.R')

con <- dbConnect( RSQLite::SQLite(), dbname="./cancer")
query <- NULL

query <- sprintf("SELECT * FROM patient")
clinicalData = dbGetQuery( con, query)

clinicalData = nominal2numeric(clinicalData)


diseases = c("LUSC")

geneFiles = "./Data/LUSC"#Transposed Gene Expression"

files <- list.files(path = geneFiles)

needed <- vector()

#Add to needed all the necessary files
for (f in files)
{
  if ( sum(grepl(diseases, f)) > 0)
  {
    needed <- c(needed, f)
  }
}

#First call whole clinicalData as first case then
#no filters

files <- needed

for (f in files)
{
  Name <- unlist(strsplit(f, "_"))[1]
  print(Name)
  clusters(geneFiles, f, Name, clinicalData = clinicalData)
}
'
for (i in 2:length(Map))
{
  select <- c( colnames(Map)[i] )
  
  for (j in 1:length(Map[, i]))
  {
    if (is.na(Map[j,i]))
    {
      break
    }
    
    value <- c( Map[j, 1] )
    
    for (f in files)
    {
      Name <- unlist(strsplit(f, "_"))[1]
      clusters(geneFiles, f, Name, select, value, clinicalData)
    }
  }
}

src = "./CLUSTERS/"

Construct(src)

Construct <- function(src = ".")
{
  dirs = list.dirs(src, full.names = FALSE, recursive = FALSE)
  files = list.files(src)[!( list.files(src) %in% dirs )]
  
  for (f in files)
  {
    if (tools::file_ext(f) != "RDATA") next
    ConstructNetwork(f, src)
  }
  
  for (dir in dirs)
  {
    path = file.path(src, dir)
    Construct(path)
  }
  #For files in argument
  #if file is not folder and is an .RDATA call constructNetwork
  #else call Construct with new folder as argument
}'