load("~/Coding Projects/Working Directories/DM_Kidney_Preprocess/Workspaces_And_Objects/eset_kidney_pc25_nsFilter.RData")
load("~/Coding Projects/Working Directories/DM Liver Preprocess/Workspaces_And_Objects/eset_liver_pc25_nsFilter.RData")

library(Biobase)
library(BiocGenerics)
library(parallel)

pdata_kidney <- pData(eset_kidney_pc.25_nsFilter_rma_qc)
pdata_liver <- pData(eset_liver_pc.25_nsFilter_rma_qc)

conditions_kidney <- pdata_kidney[ ,1:5]
conditions_liver <- pdata_liver[ ,1:5]

t_conditions_kidney <- conditions_kidney[conditions_kidney != "", ]
t_conditions_liver <- conditions_liver[conditions_liver != '', ]

# t_conditions returns treatment conditions from phenodata and various NA rows
# complete.cases() returns a logical vector with true/false as to whether a row has any NA values

# aa <- complete.cases(t_conditions_kidney)
# aa <- complete.cases(t_conditions_liver)
# summary(aa) will show what how many rows are true, should be at the beginning of the DF
# subset for rows without NA below

t_conditions_kidney <- t_conditions_kidney[1:996, ]
t_conditions_liver <- t_conditions_liver[1:1786, ]

t_conditions_kidney <- unique(t_conditions_kidney[, c("CHEMICAL", "DOSE", "VEHICLE", "ROUTE", "DURATION")])
t_conditions_liver <- unique(t_conditions_liver[, c("CHEMICAL", "DOSE", "VEHICLE", "ROUTE", "DURATION")])

t_conditions_kidney$ORGAN <- "KIDNEY"
t_conditions_liver$ORGAN <- "LIVER"

treatments <- rbind(t_conditions_liver, t_conditions_kidney)
treatments <- treatments[ ,c("ORGAN" ,"CHEMICAL", "DOSE", "VEHICLE", "ROUTE", "DURATION")]

dupsBetweenGroups <- function (df, idcol) {
  # df: the data frame
  # idcol: the column which identifies the group each row belongs to
  
  # Get the data columns to use for finding matches
  datacols <- setdiff(names(df), idcol)
  
  # Sort by idcol, then datacols. Save order so we can undo the sorting later.
  sortorder <- do.call(order, df)
  df <- df[sortorder,]
  
  # Find duplicates within each id group (first copy not marked)
  dupWithin <- duplicated(df)
  
  # With duplicates within each group filtered out, find duplicates between groups. 
  # Need to scan up and down with duplicated() because first copy is not marked.
  dupBetween = rep(NA, nrow(df))
  dupBetween[!dupWithin] <- duplicated(df[!dupWithin,datacols])
  dupBetween[!dupWithin] <- duplicated(df[!dupWithin,datacols], fromLast=TRUE) | dupBetween[!dupWithin]
  
  # ============= Replace NA's with previous non-NA value ==============
  # This is why we sorted earlier - it was necessary to do this part efficiently
  
  # Get indexes of non-NA's
  goodIdx <- !is.na(dupBetween)
  
  # These are the non-NA values from x only
  # Add a leading NA for later use when we index into this vector
  goodVals <- c(NA, dupBetween[goodIdx])
  
  # Fill the indices of the output vector with the indices pulled from
  # these offsets of goodVals. Add 1 to avoid indexing to zero.
  fillIdx <- cumsum(goodIdx)+1
  
  # The original vector, now with gaps filled
  dupBetween <- goodVals[fillIdx]
  
  # Undo the original sort
  dupBetween[sortorder] <- dupBetween
  
  # Return the vector of which entries are duplicated across groups
  return(dupBetween)
}

dupRows <- dupsBetweenGroups(treatments, "ORGAN")
treatments <- cbind(treatments, dup = dupRows)
common_treatments <- subset(treatments, dup == "TRUE")
common_treatments <- common_treatments[ ,2:6]
common_treatments <- unique(common_treatments[, c("CHEMICAL", "DOSE", "VEHICLE", "ROUTE", "DURATION")])
rownames(common_treatments) <- NULL