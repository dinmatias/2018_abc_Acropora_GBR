# NOTE:
# This script identifies a prosterior probability theshold by which a support for a certain model
#       can be considered robust (i.e., not due to random chance). Robust is defined as correctly
#       supporting a model (i.e., properly selecting the right model) at a certain probability. 
#       In this case, a 95% probability was used, which means that PODs should be correctly supported
#       95% of the time.
#
# To identify the threshold, this script utilizes the distribution of the posterior probabilities 
#       of the PODs (obtained previously using the crossValid_CL.R) and simply implement the method 
#       described by Roux et al 2016 (see page XXX). 
# In brief, it is the probability of observing a model M POD being identified as model M with probability 
#       greater than the threshold X divided by the probability of randomly supporting a POD as model M with
#       probability greater than the threshold X. P(Pm > thresX | Pm = M) / sum(Pm > thresX | Pm = i)
#
# This script uses the binary files generated from crossValid_CL.R, with the following filename structure:
#       [modelID]_crossVal_[indexStart]_[indexEnd].[species]
#       IBD_crossVal_1_500.amil : cross validation result for IBD model, PODs 1 to 500 for Acropora millepora

# All these files should be placed into a single directory, which will define dirIn below.

#############################################
#### script starts here                  ####
#############################################

# define directory
# directory containing all the binary outputs of the crossValid_CL.R (the posterior probabiltiy estimates 
#       for each POD)
dirIn <- "./amilSample/crossValid/"

# speciesID to process
species = "amil"

# list cross validation rds files
crossFiles <- list.files(
  dirIn, full.names = TRUE,
  pattern = paste(".", species, "$", sep = ""))

# get the files names without the full path
crossID <- gsub("^.*/", "", crossFiles)

# load the cross validation results
crossList <- lapply(crossFiles, function(x){
  
  # NOTE: this will load the list containing the summary of abc
  #       for each pod
  readRDS(x)
  
  })

# parse the file name to generate identifier for the cross-validation results
# this will generate the following vector:
#       [modelName]::crossValid::[podIndexStart]::[podIndexEnd]
tempID <- lapply(strsplit(
  gsub(".amil", "", crossID), "_"), 
  function(x) unlist(x))

# merge the IDs for each cross-validation file
tempID <- do.call(rbind, tempID)

# pull out these columns and convert to data.frame
# [modelName]::[podIndexStart]::[podIndexEnd]
markerID <- data.frame(model = tempID[ , 1],
                       start = as.integer(tempID[ , 3]),
                       end = as.integer(tempID[ , 4]),
                       stringsAsFactors = FALSE)

# name each element of the list (i.e., add a pod identifier for each result)
for(i in 1:length(crossList)){
  
  # for each file, pull out the start
  start = markerID$start[i]
  
  # and pull out the last index
  end = markerID$end[i]
  
  # and pull out the model id
  model = markerID$model[i]
  
  # name each pod result by [model]_POD_[podIndex]
  names(crossList[[i]]) <- paste(model, 
                                 "_POD_", 
                                 seq(start, end), 
                                 sep = "")
}

# merge the results from the different files
crossClean <- unlist(crossList, recursive = FALSE)

# list to hold the posterior probability for
# each pod
resFinal <- vector("list",
                   length = length(crossClean))

# iterate over the merged lists to pull out the
# posterior probability of each file
for(i in 1:length(crossClean)){

  # process by pod
  # this list contains results for raw estimate of posterior probability
  # and estimate using neural net
  # if all of the posterior samples belong to one model, however, then
  # no neural network can be trained, thus default to raw estimate (which is 1)
  temp <- crossClean[[i]]

  # used length to identify neural results vs non neural results
  if(length(unlist(temp)) == 12){
    
    # pull out posterior probability of neural network 
    resFinal[[i]] <- temp$neuralnet$Prob
    
    # name by the podID
    names(resFinal)[i] <- names(crossClean)[i]
    
    
  }else{
    
    # this is the raw posterior probability
    resFinal[[i]] <- temp$Prob
    
    # name with podID
    names(resFinal)[i] <- names(crossClean)[i]
    
    
  }
  
}

# merge all the results
resFinal <- do.call(rbind, resFinal)

# merge with other result identifiers: podID and modelID (which is a bit redundant)
resFinal <- data.frame(podID = row.names(resFinal),
                       modelID = gsub("_.*$", "",row.names(resFinal)),
                       resFinal,
                       stringsAsFactors = FALSE)

# the original resFinal contains columns for DIV and IBD,
#   whose values are the posterior probabilities of each model for each pod
# simply rename them here
names(resFinal)[3:4] <- c("DIV", "IBD")

############################################
#### NAIVE CONFUSION MATRIX             ####      
############################################
# generate confusion matrix based on 
# who has higher posterior probability

# identify best model
resFinal$naiveSel <- ifelse((
  resFinal$IBD - resFinal$DIV) > 0 , 
  "IBD", "DIV")

# make an id for the true model and selected model
conMatrix <- paste(resFinal$modelID, 
                   resFinal$naiveSel, sep = "__")

# table(conMatrix)

# examine the distribution
conMatrix <- table(resFinal$modelID, resFinal$naiveSel)
#      DIV  IBD
# DIV 9989   11
# IBD    5 9995

# precision for IBD
conMatrix[2, 2]/sum(conMatrix[ , 2])
# 99.89007

# precission for DIV
conMatrix[1, 1]/sum(conMatrix[ , 1])
# 99.94997

# error/misclassification rate
(conMatrix[2, 1] + conMatrix[1, 2])/ sum(conMatrix[ , ])
# 8e-04

# accuracy
(conMatrix[1, 1] + conMatrix[2, 2])/ sum(conMatrix[ , ])
# 0.9992

############################################
#### EXAMINE POSTERIOR PROBABILITY      ####
#### BASED ON ROBUSTNESS                ####      
############################################
# instead of a naive classifier, here support
# will be determined by giving a precision (i.e., robustness) threshold
# i.e., the best supported model will be considered to be
#     good only if the support (posterior probability) is 
#     higher than the probability that will give a precision
#     of greater than 95%

# convert to data table
resFinalDT <- data.table::data.table(resFinal)

# mean probability of DIV given a model
resFinalDT[ , mean(DIV), by = .(modelID)]
# modelID           V1
# 1:     DIV 0.9988971036
# 2:     IBD 0.0007638627


# mean probability of IBD given a model
resFinalDT[ , mean(IBD), by = .(modelID)]
# modelID          V1
# 1:     DIV 0.001102896
# 2:     IBD 0.999236137


# make a sequence of precision threshold that will be used to
# validate a support; e.g., if 
thresSeq <- seq(0.01, 0.99, 0.001)

ibdRobust <- vector("numeric", length = length(thresSeq))

divRobust <- vector("numeric", length = length(thresSeq))

for(ind in 1:length(thresSeq)){
  
  threshold = thresSeq[ind]
  
  # number of DIV PODs that are greater than the threshold
  #   resFinalDT[modelID == "DIV" & DIV >= threshold, .N]
  # number of DIV PODs
  #   resFinalDT[modelID == "DIV", .N]
  # numDiv: proportion of DIV PODs that are greater than the threshold
  #         this is the probability of a DIV POD being identified as a DIV POD 
  #         this is the numerator in the Roux et al 2016 paper
  numDiv <- resFinalDT[naiveSel == "DIV" & 
                         modelID == "DIV" & 
                         DIV >= threshold, .N]/resFinalDT[modelID == "DIV", .N]
  
  # number of PODs (either DIV or IBD) that will be classified as DIV given a threshold
  # tempDen <- resFinalDT[DIV >= threshold, .N, by = .(modelID)]
  tempDen <- resFinalDT[naiveSel == "DIV" & 
                          DIV >= threshold, .N, by = .(modelID)]
  
  # number of PODs per model
  numSim <- resFinalDT[ , .N, by = .(modelID)]$N
  
  # # pull out the column
  # numSim <- numSim$N 
  
  # this
  # tempDen$N/numSim: probability that a POD under model i will be classified
  #                   as DIV given a certain threshold
  #                   this is the denominator in the Roux et al 2016 paper;
  denDiv <- sum(tempDen$N/numSim)
  
  # calculate robustness
  robDiv <- numDiv/denDiv
  
  # save robustness given a probability threshold
  divRobust[ind] <- robDiv
  
  #### THIS SECTION IS SAME AS THE ABOVE, BUT FOR IBD
  numIbd <- resFinalDT[naiveSel == "IBD" & 
                         modelID == "IBD" & 
                         IBD >= threshold, 
                       .N]/resFinalDT[modelID == "IBD", .N]
  
  tempDen <- resFinalDT[naiveSel == "IBD" & 
                          IBD >= threshold, 
                        .N, by = .(modelID)]
  
  numSim <- resFinalDT[ , .N, by = .(modelID)]$N
  
  denIbd <- sum(tempDen$N/numSim)
  
  robIbd <- numIbd/denIbd
  
  ibdRobust[ind] <- robIbd 
  # robRes[[ind]] = list(robDiv, robIbd)
  
}


############################################
#### EXAMINE THE ROBUSTNESS             ####
#### (OR PRECISION) FOR EACH MODEL      ####
#### UNDER VARIOUS THRESHOLD            ####
#############################################
# the goal really is just to define the minimum
# threshold X that can give robustness of greater than
# or equal to 0.95 for each model.
# this can be done by the following:

# minimum threshold that will give at least 0.95 robustness for IBD 
thresSeq[which(ibdRobust >= 0.95)[1]]

# minimum threshold that will give at least 0.95 robustness for IBD 
thresSeq[which(ibdRobust >= 0.95)[1]]

# make a plot of theshold vs robustness
plot(x = thresSeq,
     y = c(ibdRobust),
     ylim = c(min(c(ibdRobust, divRobust)), 
              max(c(ibdRobust, divRobust))),
     type = "n",
     xlab = "threshold",
     ylab = "robustness")

# add line for the divergence model
lines(x = thresSeq,
      y = divRobust,
      col = "red")

# add line for the ibd model
lines(x = thresSeq,
      y = ibdRobust,
      col = "blue")
