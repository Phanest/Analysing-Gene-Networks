#' Divides file based on values in column
#' 
#' @param filename string, can be the name of the directory which applies
#' SubdivideFrame on all files in the directory or a specific file
#' @param Column either an index or the name of the column which SubdivideFrame
#' is going to use when dividing the table
#' @param Ignore vector, values whom we need ignore
#' @param Recursive if set to TRUE, Recursive makes the function look inside folders
#' and divides the files it finds there
#' @param skip on which row does the data begin

SubdivideFrame <- function(Column, Ignore= NULL, Recursive = FALSE, skip = 0, na.rm = FALSE)
{
  
  filename <- readClipboard()
  filename <- gsub('\\\\', '/', filename, fixed = TRUE)
  
  #Check if filename exists, if it doesn't return error
  if ( !file.exists( filename ) )
  {
    stop('No such file or folder!')
  }
  #Check if filename is a directory, if it is
  if ( dir.exists( filename ) )
  {
    files <- list.files(path= filename, recursive = Recursive)
    
    old.path = getwd()
    #Loop through files
    setwd( filename )
    lapply( files, Subdivide, Column = Column, Ignore = Ignore, skip = skip, na.rm = na.rm)
      #apply function to divide data Subdivide
    setwd( old.path )
  }
  else
  {
    Subdivide( filename, Column, Ignore, skip, na.rm)
  #if it's not apply function Subdivide
  }
}

Subdivide <- function( filename, Column, Ignore, skip, na.rm)
{
  #Read table
  data <- read.table(filename, header = TRUE, sep = '\t', skip = skip)
  
  #Check if Column is name or number
  if ( class(Column) == 'numeric')
  {
  #If name then tapply on data$Column
    by( data, data[[Column]], Export, filename, Ignore = Ignore, Column = colnames( data )[ Column ], na.rm = na.rm)
    #possible error
  }
  else
  {
    #else tapply on Column
    by( data, data[[Column]], Export, filename, Ignore, Column, na.rm)
    #check for error
  }
}

Export <- function(data, filename, Ignore, ColumnName, na.rm)
{
  filename <- file.path( filename, paste(ColumnName, data[[ColumnName]][1], sep = '' ) )
  
  name <- data[[ColumnName]][1]
  
  if ( (name %in% Ignore) || is.na(name) )
  {
    return()
  }
  
  filename <- gsub('.tsv', '', filename)
  
  dir.create( filename, recursive = TRUE )
  
  name <- paste(ColumnName, name, sep = '')
  
  name <- paste(name, '.tsv', sep = '')
  
  filename <- file.path( filename, name)
  
  write.table( data, filename, quote = FALSE, col.names = NA, sep = '\t')
}

#Are NA's read as NA's?