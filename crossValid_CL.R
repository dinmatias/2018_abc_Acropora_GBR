# can write the path of R to use here
# e.g., /bin/R

# NOTE:
# This procedure has an equivalent function in the abc package, which utilizes a proportion of the simulated
#     summary statistics as PODs (Pseudo Observed Data). However, it iterates over the PODs, and thus can be very slow in performing
#     the analysis for many PODs. Here, the estimation of posterior probability for each POD was decoupled allowing
#     it to run in parallel. Additionaly, this also enables utilization of newly simulated PODs, independent of
#     the simulated data used in training the model.
#
# This script performs posterior probabiltiy estimation using the PODs as empirical data. The
#     results of these posterior estimates will be used for the downstream crossvalidation analysis.
# The following files are needed for this script:
#     1) simulated data sets (this should be the same data sets used in estimating posterior probability
#             of the empirical data; see NOTE on estimateProb.R)
#     2) the Pseudo Observed Data (PODs) file. This file should have the same columns as the simulated
#             summary statistics file. The number of PODs should be equal to the number of rows of this
#             data table. The PODs file should be placed in a directory, whose file name can be distincly
#             identified (see modelPrefix below)
#
# The output of this script is a binary file (rds) containing the summary of the posterior probability estimation.

#### USAGE
# Rscripts crossValid_CL.R pathDir pathDirOut modelPrefix numSim podInd_start podInd_end sumstatPattern speciesPrefix 
# pathDir:        path to folder containing the pods directory and simulated data set
# pathDirOut:     path for output directory
# modelPrefix:    IMPORTANT: this should be distinct pattern for the model. This identifies which PODs to
#                 process. This is also used as prefix/suffix for the outputs.
# numSim:         number of simulated data to use in estimating pod posterior probability
# podInd_start:   start index of the pods to use
# podInd_end:     last index of the pods to use
# sumstatPattern: pattern to pull out the summary statistics to use
# speciesPrefix:  prefix for species

# Example:
# # Rscript crossValid_CL.R "./amilSample/" "./amilSample/" "DIV" "10000" 1 5 "K|R|FST" "amil"


# capture input
inputs <- commandArgs(trailingOnly = TRUE)

#############################################
#### parse input from command line       ####
#############################################
# parse the various inputs to different R objects

# dataDir <- "D:/tempFiles/amilFinal/"
dataDir <- as.character(inputs[1])

# outDir <- "D:/tempFiles/"
outDir <- as.character(inputs[2])

# mode
# modelValid <- "IBD"
modelValid <- as.character(inputs[3])

# number of simulations to use
# numSim <- 500000
numSim <- as.integer(inputs[4])

# bounds of index to use
# boundInd <- c(1, 10)
indLB <- as.integer(inputs[5])

indUB <- as.integer(inputs[6])

# summary statistcs pattern
# patternSS <- "K|R|FST"
patternSS <- as.character(inputs[7])

# species prefix/suffix
species <- as.character(inputs[8])

#### end inputs

#############################################
#### process starts here                 ####
#############################################

# load abc
library(abc)

# define vector of POD indices to process
useInd <- indLB:indUB

# define output file for the cross validation
# summary of abc object (posterior probabilities)
outFile <- paste(outDir, "/", 
                 modelValid, "_", 
                 "crossVal_", 
                 indLB, "_", 
                 indUB, ".", species, sep = "")

# define progress file
outProg <- paste(outDir, "/", 
                 modelValid, "_", 
                 "crossVal_", 
                 indLB, "_", 
                 indUB, ".", species, 
                 ".progress", sep = "")

# define file for the summary statistics name
outSum <- paste(outDir, "/",
                modelValid, "_",
                "crossVal_",
                indLB, "_",
                indUB, ".", species,
                ".sumNames", sep = "")

# define the pod files directory
dirPods <- list.files(
  dataDir, pattern = paste(
    modelValid, "pods", sep = ""),
  full.names = TRUE)

# find the pods summary statistics table
podsFiles <- list.files(dirPods,pattern = "ss",
                        full.names = TRUE)

# load pod sumamry statistics
podsData <- read.table(podsFiles, 
                       sep = "\t", 
                       header = TRUE)

# find files of the summary statistics
# tables of the two models
sumStatFiles <- list.files(dataDir, 
                           pattern = "Sum", 
                           full.names = TRUE)

# load the summary statistincs
sumStatData <- lapply(
  sumStatFiles, function(x)
    read.table(x, sep = "\t",
               header = TRUE))

# get the filenames of summary statistics table,  
# which will be used as model identifier
repNames <- gsub("^.*/", "", 
                 sumStatFiles)

# name the list
names(sumStatData) <- repNames

# remove name
rm(repNames)

# subset the pod table for the indices to use
empUse <- podsData[useInd, ]

# subset the summary statistics table by the number of 
# simulations to use
sumstatUse <- lapply(sumStatData, function(x) x[1:numSim, ])

# make vector of model indices
modelLabel <- lapply(1:length(sumstatUse), 
                     function(x) rep(names(sumstatUse)[x], 
                                     nrow(sumstatUse[[x]])))


modelLabel <- unlist(modelLabel)

# merge sum stats
sumStat <- do.call("rbind", sumstatUse)

# get column names of target sum stats
statUse <- grep(patternSS, names(empUse))

# generate holder for output
postRes <- vector("list", 
                  length = nrow(empUse))

# write a file of names of summary statistcs
write(names(sumStat[statUse]), 
      outSum, append = TRUE)

rm(sumstatUse)

rm(sumStatData)

gc()

# loop over the pods
for(i in 1:nrow(empUse)){
  
  # perform the same abc analysis using the pod as empirical data
  temp <- abc::postpr(target = empUse[i , statUse],
                      index = modelLabel,
                      sumstat = sumStat[statUse],
                      tol = 0.01,
                      method = "neuralnet",
                      maxit = 1500,
                      sizenet = 15,
                      numnet = 10,
                      MaxNWts = 2000)
  
  # summarize the results
  # will give the posterior probabilities
  postRes[[i]] <- summary(temp)
  
  # write a progress
  write(paste("Finised", i, "validation", Sys.time()), outProg, append = TRUE)
  
  # save the result as binary
  # this should be outside the loop
  # but placed here to save every output
  # (in case of unexpected system downs)
  saveRDS(postRes, file = outFile)

}





