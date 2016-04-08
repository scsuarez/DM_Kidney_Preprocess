#Preprocess to Present Calls Functions and Prep

#working directory should contain
#3 Directories
## CEL Liver Files - 2218 CEL liver files downloaded from DrugMatrix, ensure all files are present when proceeding
## Supplementary_Files - 
### Containing target of phenodatatxt currently "Supplementary_Files/Corrected_Affymetrix_Microarray_Data_Labels_spreadsheet_02_23_2012.txt"
### eventually should include 
#### desired chemical conditions for testing .xlsx file
#### detailed treatment to control mappings text table file 
## and Workspaces_and_Objects 
### directory for storing objects created to save expression sets at checkpoints in the analysis

#store current date
#today<- Sys.Date()
#format(today, format="%m-%d-%Y")

#function to download packages that I was using
#function to download packages that JL was using
#already run for Linux server, will be needed for independent (HEM) running of the code
DownloadSomePackages <- function(){
  source("http://bioconductor.org/biocLite.R")
  biocLite("affy")
  biocLite("Biobase")
  biocLite("arrayQualityMetrics")
  biocLite("genefilter")
  biocLite("rat2302.db")
  biocLite("RankProd")
  install.packages("dynamicTreeCut")
  install.packages("XML")
  install.packages("XLConnect")
}

#should run Download some packages before running these R Files
#return start date
startdate <- date()

#load packages that are used, require is similar to require() but returns an error if the package is not installed
library(Biobase) #Bioconductor - Package for analyzing MicroArray Data
library(arrayQualityMetrics)
library(XML)
library(affy)
library(genefilter)
## Not used for preprocess to FC DF, but used in analysis
#library(RankProd)
#library(rat2302.db)
#library(dynamicTreeCut)
#not currently functional on Platinum Server
# library(XLConnect)

#reusable objects to identify treatments and controls
#data tables
#ensure to use "Corrected_Affymetrix_Microarray_Data_Labels_spreadsheet_02_23_2012.txt" for this analysis
#original "Affymetrix_Microarray_Data_Labels_spreadsheet_02_23_2012.txt" has a read issue caused by prime symbol followed by -, causes ignored delimiters
phenodatatxt <- "Supplementary_Files/Corrected_Affymetrix_Microarray_Data_Labels_spreadsheet_02_23_2012.txt" #list.files(path=getwd(),full.names = TRUE,pattern = "txt")
phenodata <- read.AnnotatedDataFrame(phenodatatxt, header=TRUE, row.names="ARRAY_ID", sep="\t")

#it may seem like phenodataDF is a repeat of phenodata from an earlier line, but the table imported for phenodata has been transformed for another purpose 
#and phenodata is an annotated dataframe,(see Preprocess_to_PresentCalls.R, line 49-51)
#this reimports the table 
phenodataDF <- read.table(phenodatatxt,header=TRUE,row.names="ARRAY_ID",sep="\t")

#creating unique list of treatments and controls
#getting phenodata for treatments to make a list of treatment/condition types
pdata_treatments <- phenodataDF[phenodataDF$CHEMICAL != "", ] # treatments
pdata_controls <- phenodataDF[phenodataDF$CHEMICAL == "", ] # controls
pdata_treatments_liver <- subset(pdata_treatments, ORGAN_OR_CELL_TYPE=="LIVER")
pdata_controls_liver <- subset(pdata_controls, ORGAN_OR_CELL_TYPE=="LIVER")
uTreatments <- unique(pdata_treatments_liver[, c("CHEMICAL", "DOSE", "DURATION", "VEHICLE", "ROUTE", "ORGAN_OR_CELL_TYPE")])
uControls <- unique(pdata_controls_liver[, c("CHEMICAL", "DOSE", "DURATION", "VEHICLE", "ROUTE", "ORGAN_OR_CELL_TYPE")])

#function to get absolute file path for all CEL files in a directory, filter out duplicate arrays, no copying of files
#a directory should be used that has includes all CEL files but does not have any CEL files that are not being tested
#returns one absolute file path for each CEL file in directory and subdirectories
GetAllCELFiles <- function(directory=getwd()){
  allCELPaths <- list.files(pattern="CEL",path=directory,recursive=TRUE,full.names=TRUE)  #get file paths for all CEL files in directory and subdirectories 
  CELNames <- basename(allCELPaths)	#get just the last part of the file path-the .CEL file name
  uCELNames <- unique(CELNames) 	#getting a unique list since the same files were in multiple folders
  uniqueCELPaths <- " " 		#creating object for use in loop
  for (count in 1:length(uCELNames)){ #once for each CEL file
    tempIndex <- grep(uCELNames[count],allCELPaths,fixed=TRUE)#finds indices for that CEL file in the CEL path names
    uniqueCELPaths[count] <- allCELPaths[tempIndex[1]]	#takes only the first instance to prevent duplicates (the files should be the same)
  }
  return(uniqueCELPaths)
}

#Older comments from JTL
#function to make mappings between experimental and control arrays based on Drugmatrix annotation data
#input-should be a data frame of Drugmatrix annotation data such as from the file "Affymetrix_Microarray_Data_Labels_spreadsheet_02_23_2012" which was obtained from the Drugmatrix downloads webpage 
#output- ArrayMatch-data frame that has info for an experimental condition matched with information for a variable number of control arrays including the name of the control array
#ArrayMatch-should have multiple lines with the same experimental condition with one line for each control match, the same control may be matched to multiple experimental conditions but will not be matched multiple times to the same condition... 
#...(continued from above) so list of controls for each condition should not contain duplicates

#new function for making matches between treatments and corresponding controls
#inputs - Dataframe objects in the order (treatments, controls)
#output - combined data frame that matches each treatment CEL file to any matching control CEL files
#uses inputs like pdata_treatments/controls or uTreatments/uControls for objects of functions
ArrayMatch <- function(treatments, controls){
  #extracting row names and pasting them to the 1st column so they will be actual data instead of data labels
  treatments <- cbind(row.names(treatments), treatments)
  #eliminating row names
  row.names(treatments) <- NULL
  #rename 1st column to ARRAY_ID
  names(treatments)[1] <- "ARRAY_ID"
  #append TREATMENT_ to all treatment columns
  colnames(treatments) <- paste("TREATMENT_", colnames(treatments), sep="")
  
  #making table that matches experimental arrays and control arrays
  #table/data frame will contain information from treatments and information from controls, current usage will be for _liver
  ArrayMatch <- data.frame(1:(ncol(treatments)+ncol(controls)))#just need to preallocate the right amount of columns, might not be best way
  ArrayMatch <- t(ArrayMatch) #transposing since it was just one column and what were the row names need to be the column names  
  ArrayMatch <- as.data.frame(ArrayMatch)
  colnames(controls) <- paste("CONTROL_",colnames(controls),sep="") #adding "CONTROL_" to the column names of the object with control information
  colnames(ArrayMatch) <- c(colnames(treatments),colnames(controls)) #replacing old column names with new names (exp names did not change)
  CONTROL_ARRAY_ID <- "preallocate" #just to preallocate
  for (count1 in 1:nrow(treatments)){ #do once for each treatment
    for (count2 in 1:nrow(controls)){ #check each control array to see if it matches the current treatment
      #if the treatment matches the control make a new line in ArrayMatch with that experimental condition matched with control information
      if (treatments[count1,"TREATMENT_DURATION"]==controls[count2,"CONTROL_DURATION"] & treatments[count1,"TREATMENT_VEHICLE"]==controls[count2,"CONTROL_VEHICLE"] & treatments[count1,"TREATMENT_ROUTE"]==controls[count2,"CONTROL_ROUTE"]){
        tempLine <- cbind(treatments[count1,], controls[count2,])   #the condition and the control information stuck together as one line
        ArrayMatch <- rbind(ArrayMatch, tempLine)    #that line added to the data frame for matched controls and conditions
        CONTROL_ARRAY_ID <- c(CONTROL_ARRAY_ID, row.names(controls)[count2])     #the name of the control array added to list of control array names, the list of control array names will later be added to the data frame
        
      }
    }
  }
  ArrayMatch <- cbind(ArrayMatch, CONTROL_ARRAY_ID)#treatment and control array names are added to the data frame
  ArrayMatch <- ArrayMatch[-1,] #taking off irrelevant row
  row.names(ArrayMatch) <- NULL #making sure we have no erroneous row names
  
  return(ArrayMatch)
}

#making data frame that matches experimental conditions with appropriate controls
##"this does not have to be done at this specific point in the workflow, it could be somewhere else" -from JTL
#adjusted from JTL function, don't know how to handle the fact that multiple conditions are repeated
#this is a step we need before the fold change matrix is calculated
#only needs to be performed once at the beginning of the analysis and will not change as long as we are using DrugMatrix
#do this before actual array analysis
Liver_ConditionsMatch_AllC <- ArrayMatch(pdata_treatments_liver, pdata_controls_liver)
Liver_ConditionsMatch <- cbind(Liver_ConditionsMatch_AllC[1:6], Liver_ConditionsMatch_AllC[20:22], Liver_ConditionsMatch_AllC[34])
uLiver_ConditionsMatch <- ArrayMatch(uTreatments, uControls)
write.table(Liver_ConditionsMatch, file = "Supplementary_Files/Liver_Treated_to_Control_Mapping.txt", sep="\t")