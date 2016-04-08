# Should have run Prep_Preprocess_to_PresentCalls.R before this script, otherwise certain objects will not be present
# working directory "/Working Directories/Preprocess to Present Calls/"

# load CEL files (ReadAffy) and read annotations
uniqueCELPaths <- GetAllCELFiles() # getting all CEL files from various folders, recursive (goes into subfolders too),user defined function-JL
eset_liver_raw <- ReadAffy(filenames=uniqueCELPaths) # load CEL files 
allCELNames <- basename(uniqueCELPaths) # no path listed if one folder has all cel files
phenodata_raw <- phenodata[allCELNames,] # getting just the rows for CEL files that are being used
pdata_raw <- pData(phenodata_raw)   # getting data frame of the phenodata
eset_liver_raw <- eset_liver_raw[,allCELNames] # selects columns for CEL files that are being used
phenoData(eset_liver_raw) <- phenodata_raw # attaching phenoData to eset_liver_raw
save(eset_liver_raw, file="Workspaces_And_Objects/eset_liver_raw.RData")

# do RMA normalization
eset_liver_raw_rma <- rma(eset_liver_raw) # rma converts to log2 as part of normalization
remove(eset_liver_raw) # remove object from environment to free up memory since it is no longer needed
# can do this because it was saved in "~Workspaces_and_Objects/"
save(eset_liver_raw_rma, file="Workspaces_And_Objects/eset_liver_raw_rma.RData")

# do ARRAYQUALITYCONTROL
inttest <- c("CHEMICAL","DURATION","DOSE")# changed-JL, old was inttest <- c("COMPOUND_NAME","SACRI_PERIOD","DOSE")
arrayQualityMetrics(expressionset = eset_liver_raw_rma, outdir = "./Quality_Metrics/", force = FALSE, intgroup = inttest)
remove(eset_liver_raw_rma) # remove object from environment to free up memory since it is no longer needed and was saved previously
# *****************************
# 
# comment from Dr. AH, not used now- # new annotation file: 3_Affy_minus150_clean2.txt

# obtain outlier data from arrayQualityMetrics using XML package 
tables <- readHTMLTable("Quality_Metrics/index.html")# JL,read outliers from html file
# tables is a list, next step checks how many rows each element has and grabs the one with the most rows which should be the desired one
outlierDF <- tables[[which.max(sapply(tables,nrow))]] # JL, grabbing table that has three outlier criteria (also should have the most rows)
outliers <- outlierDF[,4]=="x" | outlierDF[,5]=="x" | outlierDF[,6]=="x" # JL,checking whether array has any of criteria

# taking outliers out of list of names
outlierNames <- outlierDF[outliers,3] # grabbing the .CEL names of outliers
notOutlierFiles <- uniqueCELPaths # JL,outliers will be taken out but are currently included, GetAllCELFiles() was used earlier to get uniqueCELPaths
for (count in 1:length(outlierNames)){ # JL, once for each outlier
  notOutlierFiles <- notOutlierFiles[-grep(outlierNames[count],notOutlierFiles,fixed=TRUE)] # JL,take out outliers one by one
}
# transform outlierNames from a factor to a vector so we can save the list of CEL files removed as outliers
outlierNamesVec <- as.vector(outlierNames)
# save list of outliers to be removed
write(outlierNamesVec, file = "Supplementary_Files/Outlier_CELnames_removed.txt")

# load in files that are not outliers and reattach phenodata
affy_liver_qc <- ReadAffy(filenames = notOutlierFiles)
notOLNames <- basename(notOutlierFiles)
phenodata_notOL <- phenodata[notOLNames,] # getting just the rows for CEL files that are being used
pdata_notOL <- pData(phenodata_notOL)
affy_liver_qc <- affy_liver_qc[,notOLNames] # changing order
phenoData(affy_liver_qc) <- phenodata_notOL
save(affy_liver_qc, file="Workspaces_And_Objects/affy_liver_qc.RData")

# use affybatch object to get present calls
# do this immediately after qc to prevent error of unknown origin
callsESet <- mas5calls(affy_liver_qc)   # on data without outliers
callsDF <- exprs(callsESet)             # dataframe of calls

# getting percent present for each probe in our expression data
# input - callsDF, a matrix that is the output of mas5calls and exprs on an affybatch object
# output - a labeled list of probes with the percentage present for each probe
# transform callsDF from P/M/A into TRUE/FALSE matrix, TRUE if P/FALSE if other
prescallsDF <- callsDF=="P"
# transform TRUE/FALSE to 1/0
# could also use matrix(), but this preserves row and col names
numprescallsDF <- prescallsDF*1
# returns a large numeric list, sums TRUES and divides by number of CEL files
# result is a list of probes and corresponding percent present
percPresent <- rowSums(numprescallsDF)/dim(numprescallsDF)[2]
# use this list to filter our final fold change matrix
# old line to do this, based on MAHC.R line 104
# saving below for later
# outputfoldchangematrix <- foldchangematrix[percPresent>=.25,] # select for rows/genes that had greater than 25% of conditions with all present calls

# do rma for this set
eset_liver_rma_qc <- rma(affy_liver_qc)
# filters out nonspecific probes
list_liver_nsFilter_rma_qc <- nsFilter(eset_liver_rma_qc, remove.dupEntrez = FALSE)
# nsfilter returns a 2 element list, of the eset and the filter log
# use exprs to return just the expression set portion of this object
eset_liver_nsFilter_rma_qc <- list_liver_nsFilter_rma_qc$eset
save(eset_liver_nsFilter_rma_qc, file="Workspaces_And_Objects/eset_liver_nsFilter.RData")

################################
# calculate FC matrix
#################################

# going back to version created by MDAH
# adapted from DM_Liver_Preprocess.R
# Read the entire dataset
# load("eset_liver_nsFilter.RData.RData")
#creating a function to make the fold change data frame
# 3 inputs
## eset - eset after gene filter fitlering of non specific probes and probes with low variance
## mappings - DF of detailed treatment to control mappings
## filename_fc - file path of where to save the FC matrix that is created
# output - fold change data frame

# previous arguments used by MDAH
# Read treatment to control mappings
## mappings <- read.delim("Supplementary_Files/Liver_Treated_to_Control_Mapping.txt", header=TRUE, sep="\t", row.names = NULL)
# filename to be saved  
##filename_fc <- paste("liver_fc_gene_filter_rma_qc.txt", sep="")

makeFC_DF <- function(eset, mappings, filename_fc){
  pdata <- pData(eset)
  pdataSorted <- pdata[order(pdata$CHEMICAL, pdata$DOSE, pdata$DURATION, pdata$VEHICLE, pdata$ROUTE), ]
  sampleNamesSorted <- row.names(pdataSorted)
  esetSorted <- eset[,sampleNamesSorted]
  eset <- esetSorted
  pdata <- pdataSorted
  fc_final <- exprs(eset)[,0]
  num_replicates_treatments <- numeric(0)
  num_replicates_controls <- numeric(0)
  
  pdata_c <- pdata[pdata$CHEMICAL == "", ] # controls
  pdata_t <- pdata[pdata$CHEMICAL != "", ] # treatments
  
  treatments <- unique(pdata_t[, c("CHEMICAL", "DOSE", "DURATION", "VEHICLE", "ROUTE")])
  for (i in 1:nrow(treatments)) {
    print("-----------------------------------------------------------------------------------")
    i_print <- paste(treatments[i, "CHEMICAL"], treatments[i, "DOSE"], treatments[i, "DURATION"], treatments[i, "VEHICLE"], treatments[i, "ROUTE"], sep="|")
    pdata_t_current <- pdata_t[treatments[i, "CHEMICAL"] == pdata_t$CHEMICAL & treatments[i, "DOSE"] == pdata_t$DOSE & treatments[i, "DURATION"] == pdata_t$DURATION & treatments[i, "VEHICLE"] == pdata_t$VEHICLE & treatments[i, "ROUTE"] == pdata_t$ROUTE,] # data for treatment
    treatmentNames <- row.names(pdata_t_current) 
    c.mappings <- mappings[mappings$TREATMENT_ARRAY_ID %in% treatmentNames,] # Map controls
    controlNames <- as.character(unique(c.mappings$CONTROL_ARRAY_ID))
    controlNames <- controlNames[controlNames %in% row.names(pdata_c)] # make sure these are in eset
    print(paste("Number of treatments for ", i_print, ": ", length(treatmentNames), sep=""))
    print(paste("Number of controls for ", i_print, ": ", length(controlNames), sep=""))
    num_replicates_treatments <- c(num_replicates_treatments, length(treatmentNames))
    num_replicates_controls <- c(num_replicates_controls, length(controlNames))
    eset_c <- eset[,c(controlNames)] # eset with controls
    exprs_c <- exprs(eset_c)
    eset_t <- eset[,c(treatmentNames)] # eset with treatments
    exprs_t <- exprs(eset_t)
    fc <- data.frame(rowMeans(exprs_t) - rowMeans(exprs_c))
    colnames(fc) <- i_print
    fc_final <- cbind(fc_final, fc)
  }
  
  fc_final <- cbind(row.names(fc_final), fc_final)
  colnames(fc_final)[1] <- "Probe_Set_ID"
  
  print(paste("Saving: ", filename_fc, sep=""))
  write.table(fc_final, file=filename_fc, sep="\t", quote=FALSE, row.names = FALSE)
  
  return(fc_final)
#************************
}

filename_fc <- paste("Workspaces_And_Objects/liver_fc_gene_filter_rma_qc.txt", sep="")
fc_initial <- makeFC_DF(eset_liver_nsFilter_rma_qc, Liver_ConditionsMatch, filename_fc )

# function to filter rows of percPresent vector and match with rows of the fold change matrix
# will be used to filter the fc matrix based on present calls
# inputs - 
## percPresent vector, named list of probes and the % of the probes present for each probe 
## fc_final data frame - created df of probes and chemical treatment to control log ratios
# output - a filtered vector that matches the dimensions and names of the fc matrix with the percPresent vector
## ie, returns only the probes from percPresent that are also in the fold change matrix
presCallsFCMatch <- function(percPresent, fc_DF){
  probenames <- names(percPresent)
  probenames_filter <- row.names(fc_DF)
  # if the probe name(row name) of the fold change matrix is present in the list of all probes from the present calls, return TRUE
  # return a logical vector indicating if the probename from the fold change matrix is present in the complete list of probes
  percPresent_filter <- probenames %in% probenames_filter
  # transform from logical to numeric
  percPresent_filter <- percPresent_filter*1
  # preserve the probes with 0 present by transforming, add 1 to all values
  percPresent_num <- percPresent+1
  # multiply percPresent_num by percPresent_filter to set to 0 all probes not present in the fc matrix
  percPresent_removed <- percPresent_filter*percPresent_num
  # remove probes not present in fc matrix
  # values of zero are not present in our fc matrix  
  percPresent_removed <- percPresent_removed[percPresent_removed!=0]
  # retransform our percPresent vector to original values, now with only probes from fc_final
  percPresent_removed <- percPresent_removed-1

  return(percPresent_removed)
}

#match rows of fc_initial with percPresent
#filters out rows of percPresent that aren't
percPresent_FCMatch <- presCallsFCMatch(percPresent, fc_initial)

#filter inital fold change data frame probes
#filter only probes with greater than %25 present
fc_final <- fc_initial[percPresent_FCMatch >= .25, ]
filename_fc_final <- paste("Workspaces_And_Objects/liver_fc_final_filtered.txt", sep="")
write.table(fc_final, file=filename_fc_final, sep="\t", quote=FALSE, row.names = FALSE)

enddate <- date()