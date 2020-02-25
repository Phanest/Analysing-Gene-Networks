library('RSQLite')
require(WGCNA)

allowWGCNAThreads()

# src is folder where gene data is stored
# file is name of file
# filterName is a vector with names of variables we want to filter with
# a NULL value means all
# filterValue is a list of vectors where every list corresponds with some name
# ClinicalData is the clinical data in its numeric form in the form of a dataframe

#The variables we wish to see in the graph
#selection= "clinical_data_pathologic_stage,clinical_data_pathology_M_stage,
#clinical_data_pathology_N_stage,clinical_data_pathology_T_stage,clinical_data_gender,clinical_data_years_to_birth
#,ips_ctla4_neg_pd1_neg,ips_ctla4_neg_pd1_pos,ips_ctla4_pos_pd1_neg,ips_ctla4_pos_pd1_pos"

selection= "clinical_data_gender,clinical_data_clinical_stage,clinical_data_pathologic_stage,clinical_data_pathology_M_stage,
clinical_data_pathology_N_stage,clinical_data_pathology_T_stage,clinical_data_years_to_birth
,ips_ctla4_neg_pd1_neg,ips_ctla4_neg_pd1_pos,ips_ctla4_pos_pd1_neg,ips_ctla4_pos_pd1_pos,clinical_data_days_to_death,clinical_data_number_of_lymph_nodes,clinical_data_year_of_tobacco_smoking_onset,clinical_data_tumor_size,
clinical_data_progression_free,clinical_data_PAM50MRNA,clinical_data_NRAS,clinical_data_NF1,clinical_data_msi_status,
clinical_data_karnofsky_performance_score,clinical_data_HER2,clinical_data_extrathyroidal_extension,clinical_data_psa_value,
clinical_data_gleason_score"

clusters <- function(src, file, destination, filterName = NULL, filterValue = NULL, clinicalData)
{
  
  f = file
  
  #print(sprintf("%s %s %s", src, f, is.null(flt) ) );
  
  #TODO: Maybe add a destination file argument
  #Take info about columns from arguments
  #Test filter function
  #How to automatically cut clusters
  #Cut tree, show red line
  
  #=====================================================================================
  #
  
  #  Code chunk 1
  #
  #=====================================================================================
  
  
  # Display the current working directory
  oldwd = getwd();
  # If necessary, change the path below to the directory where the data files are stored. 
  # "." means current directory. On Windows use a forward slash / instead of the usual \.
  workingDir = src;
  setwd(workingDir); 
  # Load the WGCNA package
  library(WGCNA);
  # The following setting is important, do not omit.
  options(stringsAsFactors = FALSE);
  #Read in the female liver data set
  library(readr)
  femData = read_csv(f);
  # Take a quick look at what is in the data set:
  dim(femData);
  names(femData);
  
  #Works
  #=====================================================================================
  #
  #  Code chunk 2
  #
  #=====================================================================================
  
  #if not transposed
  
  #datExpr0 = as.data.frame(t(femData[, -c(1:2)]));
  #names(datExpr0) = femData$substanceBXH;
  #rownames(datExpr0) = names(femData)[-c(1:2)];
  
  #Remove the gene column, so that our dataframe can be numeric
  datExpr0 = femData[, -1]
  colnames(datExpr0) = strtrim(colnames(datExpr0) , 12);
  #Transpose our data so columns are genes
  datExpr0 = t(datExpr0)
  
  #Colnames and rownames
  colnames(datExpr0) = femData$Gene
  
  datExpr0 = as.data.frame(datExpr0);
  #names(datExpr0) = names(femData)[-c(1)];
  #rownames(datExpr0) = strtrim(rownames(datExpr0) , 12);

  #For testing purposes
  #View(datExpr0[, -c(1500)])
  #datExpr0 = datExpr0[ -c(3000), -c(1500)]
  
  #=====================================================================================
  #
  #  Code chunk 3
  #
  #=====================================================================================
  
  #Check if too many missing values and outliers via magic
  
  gsg = goodSamplesGenes(datExpr0, verbose = 3);
  gsg$allOK
  
  
  #=====================================================================================
  #
  #  Code chunk 4
  #
  #=====================================================================================
  
  
  if (!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  }
  
  #=====================================================================================
  #
  #  Code chunk 5
  #
  #=====================================================================================
  
  fName = ""
  
  if (is.null(filterName))
  {
    fName = "all"
  } else
  {
    for (i in 1:length(filterName))
    {
      values = paste(filterValue[[i]], collapse = ",")
      fName = sprintf("%s %s;%s", fName, filterName[i], values)
    }
  }
  
  fName = file
  print(fName)
  
  sampleTree = hclust(dist(datExpr0), method = "average");
  
  #TODO
  
  dev.new()
  
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.
  sizeGrWindow(12,9)
  #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
  #par(cex = 0.6);
  #par(mar = c(0,4,2,0))
  
  plot(sampleTree, main = sprintf("Sample clustering to detect outliers %s", fName), sub="", xlab="", cex.lab = 1.5, 
       cex.axis = 1.5, cex.main = 2)
  #=====================================================================================
  #
  #  Code chunk 6
  #
  #=====================================================================================
  
  height <- as.numeric(readline('Enter height cut-off: '))
  while(is.na(height) || !is.numeric(height)){
    height= as.numeric(readline(prompt="Enter height cut-off: "))
  }
  
  dev.off()
  
  dev.new()
  
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.
  sizeGrWindow(12,9)
  #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
  #par(cex = 0.6);
  #par(mar = c(0,4,2,0))
  dir.create('./CLUSTERS')
  dir.create(sprintf('./CLUSTERS/%s', file))
  setEPS()
  postscript(sprintf("./CLUSTERS/%s/Sample clustering to detect outliers %s.eps", file, fName))
  plot(sampleTree, main = sprintf("Sample clustering to detect outliers %s", fName), sub="", xlab="", cex.lab = 1.5, 
       cex.axis = 1.5, cex.main = 2)
  
  # Plot a line to show the cut
  abline(h = height, col = "red");
  dev.off()
  # Determine cluster under the line
  clust = cutreeStatic(sampleTree, cutHeight = height, minSize = 10)
  table(clust)
  # clust 1 contains the samples we want to keep.
  keepSamples = (clust==1)
  datExpr = datExpr0[keepSamples, ]
  
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  
  #print(sum(clust==1))
  #print(sum(clust==0))
  #=====================================================================================
  #
  #  Code chunk 7
  #
  #=====================================================================================
  
  attributes <- NULL
  
  #dim(clinicalData)
  #names(clinicalData)
  
  # remove columns that hold information we do not need.
  barcodes = clinicalData$barcode;
  
  allTraits <- NULL
  
  #if no filter don't filter
  #if filter, for all values filter clinicalData
  
  if (is.null(filterName))
  {
    allTraits <- clinicalData
  }
  else
  {
    allTraits <- clinicalData
    n = length(filterName)
    for (i in 1:n)
    {
      values = filterValue[[i]] #values for filterName[i]
      col = allTraits[, filterName[i] ] #take clinical data from filterName[i] column
      allTraits <- allTraits[col %in% values, ] #take only rows whose values are in var: values
    }
  }
  
  barcodes = allTraits$barcode;
  
  #Select attributes we wish to show and compare against
  'attributes <- unlist(strsplit(selection, ","))
  attributes <- lapply(attributes, trimws)
  attributes <- c(attributes)
  attributes <- unlist(attributes)
  allTraits = allTraits[, c(1, which(names(allTraits) %in% attributes))];
  '
  
  attributes <- unlist(strsplit(selection, ","))
  attributes <- lapply(attributes, trimws)
  attributes <- c(attributes)
  attributes <- unlist(attributes)
  allTraits = allTraits[, c(1, which(names(allTraits) %in% attributes))];
  
  # Form a data frame analogous to expression data that will hold the clinical traits.
  
  #patient names
  femaleSamples = strtrim( rownames(datExpr), 12);
  
  traitRows = barcodes %in% femaleSamples
  
  print(dim(allTraits))
  datTraits = data.frame( allTraits[traitRows, -1]);
  names(datTraits) = names(allTraits)[-1]
  rownames(datTraits) = allTraits[traitRows, 1]
  #rownames(datTraits) = gsub("\\.", "-", make.names( allTraits[traitRows, 1], TRUE) );
  
  exprRows = femaleSamples %in% barcodes
  datExpr = datExpr[exprRows, ]
  print(dim(datTraits))
  collectGarbage();
  
  #=====================================================================================
  #
  #  Code chunk 8
  #
  #=====================================================================================
  
  #print(fName)
  #View(clinicalData)
  #View(allTraits)
  #View(datTraits)
  #print(length(traitRows))
  assign("datExpr", datExpr, envir = .GlobalEnv)
  assign("datTraits", datTraits, envir = .GlobalEnv)
  #readline(prompt = "Continue?")
  #return(datTraits)

  tryCatch( {
    
    # Re-cluster samples
    sampleTree2 = hclust(dist(datExpr), method = "average")
    # Convert traits to a color representation: white means low, red means high, grey means missing entry
    
    colLen = length(colnames(datTraits))
    for (i in 1:colLen) {
      if (typeof(datTraits[,i]) != 'double') {
        datTraits[, i] = as.numeric(datTraits[, i])
      }
    }
    
    removeColumns = names(which(is.na(colSums(datTraits))))
    
    datTraits = datTraits[, -which(colnames(datTraits) %in% removeColumns)]
    
    traitColors = numbers2colors(datTraits, signed = FALSE)
    
    dev.new()
    
    setEPS()
    postscript(sprintf("./CLUSTERS/%s/Sample dendrogram and trait heatmap %s.eps", file, fName))
  # Plot the sample dendrogram and the colors underneath.
    plotDendroAndColors(sampleTree2, traitColors,
                        groupLabels = names(datTraits),
                        main = sprintf("Sample dendrogram and trait heatmap %s", fName))
    dev.off()
    #=====================================================================================
    #
    #  Code chunk 9
    #
    #=====================================================================================
    
    #dir.create("./CLUSTERS")
    #dir.create( sprintf("./CLUSTERS/%s", file) )
    save(datExpr, datTraits, file = sprintf("./CLUSTERS/%s/%s.RDATA", file, file) )
    
    setwd(oldwd)
  },
  error= function(cond) {
    
    dir.create("./CLUSTERS")
    dir.create( sprintf("./CLUSTERS/%s", file) )
    dir.create( sprintf("./CLUSTERS/%s/%s", file, fName))
    write( sprintf("filters: %s, length: %d, cond: %s", fName, length(datTraits[,1]), cond),
                 sprintf("./CLUSTERS/%s/%s/%s.txt",
                         file,
                         fName,
                         fName))
    write.table(allTraits,
                sprintf("./CLUSTERS/%s/%s/DataTraits.txt",
                        file,
                        fName),
                sep="\t")
    setwd(oldwd)
    })
  
}

nominal2numeric2 <- function(dataset)
{
  
  a <- unique(changeMStage)
  b <- unique(changeNStage)
  c <- unique(changeTStage)
  d <- unique(changePathologicStage) #also Clinical Stage
  e <- c("female", "male")
  
  dataset$clinical_data_pathology_M_stage <- as.numeric (factor(
    replaceValues(dataset$clinical_data_pathology_M_stage, 
                  findMStage, 
                  changeMStage),
    levels= a) )
  
  dataset$clinical_data_pathology_N_stage <- as.numeric( factor(
    replaceValues(dataset$clinical_data_pathology_N_stage, 
                  findNStage, 
                  changeNStage),
    levels = b) )
  
  dataset$clinical_data_pathology_T_stage <- as.numeric( factor(
    replaceValues(dataset$clinical_data_pathology_T_stage,
                  findTStage,
                  changeTStage),
    levels = c) )
  
  dataset$clinical_data_pathologic_stage <- as.numeric( factor(
    replaceValues(dataset$clinical_data_pathologic_stage, 
                  findPathologicStage, 
                  changePathologicStage),
    levels = d) )
  
  dataset$clinical_data_gender <- as.numeric (factor(
    dataset$clinical_data_gender,
    levels= e) )
  
  #NEW
  
  PAM50 = unique(dataset$clinical_data_PAM50MRNA)[-1]
  NRAS = unique(dataset$clinical_data_NRAS)[-1]
  NF1 = unique(dataset$clinical_data_NF1)[c(-1,-9)]
  MSI = unique(dataset$clinical_data_msi_status)[-1]
  EE = unique(dataset$clinical_data_extrathyroidal_extension)[-1]
  
  dataset$clinical_data_clinical_stage <- as.numeric (factor(
    dataset$clinical_data_clinical_stage,
    levels= d) )
  
  dataset$clinical_data_PAM50MRNA <- as.numeric (factor(
    dataset$clinical_data_PAM50MRNA,
    levels= PAM50) )
  
  dataset$clinical_data_NRAS <- as.numeric (factor(
    dataset$clinical_data_NRAS,
    levels= NRAS) )
  
  dataset$clinical_data_NF1 <- as.numeric (factor(
    dataset$clinical_data_NF1,
    levels= NF1) )
  
  dataset$clinical_data_msi_status <- as.numeric (factor(
    dataset$clinical_data_msi_status,
    levels= MSI) )
  
  dataset$clinical_data_extrathyroidal_extension <- as.numeric (factor(
    dataset$clinical_data_extrathyroidal_extension,
    levels= EE) )
  
  #Gender
  
  GENDER = unique(dataset$clinical_data_gender)
  
  dataset$clinical_data_gender <- as.numeric (factor(
    dataset$clinical_data_gender,
    levels= GENDER) )
  
  n <- max(length(a), length(b), length(c), length(d), length(e),
           length(PAM50), length(NRAS),length(NF1),length(MSI),length(EE))
  
  length(a) <- n
  length(b) <- n
  length(c) <- n
  length(d) <- n
  length(e) <- n
  length(PAM50) <- n 
  length(NRAS) <- n
  length(NF1) <- n
  length(MSI) <- n
  length(EE) <- n
  
  
  num <- data.frame( id=c(1:n), 
                     clinical_data_pathology_M_stage=a,
                     clinical_data_pathology_N_stage=b,
                     clinical_data_pathology_T_stage=c,
                     clinical_data_pathologic_stage=d,
                     clinical_data_gender=e,
                     NRAS=NRAS,
                     PAM50=PAM50,
                     NF1=NF1,
                     MSI=MSI,
                     EE=EE)
  dir.create("./Map")
  write.csv(num, "./Map/Map.csv")
  
  assign("Map", num, envir = .GlobalEnv)
  
  return(dataset)
}

nominal2numeric <- function(dataset)
{
  
  a <- unique(changeMStage)
  b <- unique(changeNStage)
  c <- unique(changeTStage)
  d <- unique(changePathologicStage)
  e <- c("female", "male")
  
  dataset$clinical_data_pathology_M_stage <- as.numeric (factor(
    replaceValues(dataset$clinical_data_pathology_M_stage, 
                  findMStage, 
                  changeMStage),
    levels= a) )
  
  dataset$clinical_data_pathology_N_stage <- as.numeric( factor(
    replaceValues(dataset$clinical_data_pathology_N_stage, 
                  findNStage, 
                  changeNStage),
    levels = b) )
  
  dataset$clinical_data_pathology_T_stage <- as.numeric( factor(
    replaceValues(dataset$clinical_data_pathology_T_stage,
                  findTStage,
                  changeTStage),
    levels = c) )
  
  dataset$clinical_data_pathologic_stage <- as.numeric( factor(
    replaceValues(dataset$clinical_data_pathologic_stage, 
                  findPathologicStage, 
                  changePathologicStage),
    levels = d) )
  
  dataset$clinical_data_gender <- as.numeric (factor(
    dataset$clinical_data_gender,
    levels= e) )
  
  n <- max(length(a), length(b), length(c), length(d), length(e))
  
  length(a) <- n
  length(b) <- n
  length(c) <- n
  length(d) <- n
  length(e) <- n
  
  num <- data.frame( id=c(1:n), 
                     clinical_data_pathology_M_stage=a,
                     clinical_data_pathology_N_stage=b,
                     clinical_data_pathology_T_stage=c,
                     clinical_data_pathologic_stage=d,
                     clinical_data_gender=e )
  dir.create("./Map")
  write.csv(num, "./Map/Map.csv")
  
  assign("Map", num, envir = .GlobalEnv)
  
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