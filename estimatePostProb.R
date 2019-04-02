# NOTE on file names and structure:
#   1) For each model, there must be two table files (tab delimited) --
#             one for the summary statistics and one for the parameters;
#             this directory should also contain the empircal data summary statistics,
#             in the same format as the simulated summary statistics table (same header, or column names)
#   2) These files should be named so that the model can be identifed, and the
#             parameter file should have "Param" characters, the summary statistics file should have
#             "Sumstat", and the emperical should have "Empirical".
#             As an example, for the Acropora millepora data set, these are the files in the directory:
#             divParam_amilFinal      divergence parameter file
#             divSumstat_amilFinal    divergence summary statistics file
#             ibdParam_amilFinal      ibd paramater file
#             ibdSumstat_amilFinal    ibd summary statistics file
#             amilEmpirical           empirical sumary statistics


# full path to the directory containing the 
# files mentioned above
dataDir <- "./sample/amilSample/"

# this specifies the number of simulations to use for each model
# need to do this one because the original number of simulations
# was not reached for IBD -- hence, need to scale down to 500K
numSim <- 10000

# this correonds to the pattern of summary statistcs to be used
# K = number of alleles | R = range of allele size | Fst = differentiation
sumstatPat <- "K|R|FST"

# name to suffix the output; here it's amil for acropora millepora outputs
species = "amil"

#############################################
#### load data set                       ####
#############################################
# this section simply loads the data set and this
# is very specific to how the files are structured
# and how they are named (see description above)


# set the working directory to that directory
setwd(dataDir)

# get the date today
dataRun <- gsub("\\-" , "", Sys.Date())

# file name for the pca plot
pcaFile <- paste(dataDir, "/pcaPlot", sep = "")

# file name for the estimation of the
# psoterior probability
postprNeural <- paste(dataDir, "/", dataRun, "_",
                      "neuralModelsel", ".",
                      species, sep = "")

# list all the paramter files
paramFiles <- list.files(dataDir, 
                         pattern = "Param", 
                         full.names = TRUE)

# list all the summary statistics file
sumStatFiles <- list.files(dataDir, 
                           pattern = "Sum", 
                           full.names = TRUE)

# list the full path of the empirical file
empFile <- list.files(dataDir, 
                      pattern = "Emp", 
                      full.names = TRUE)

# load the parameter files
paramData <- lapply(paramFiles, 
                    function(x){
                      
                      read.table(x, sep = "\t", 
                                 header = TRUE)
                      
                    })

# rename the parameter list with 
repNames <- gsub("^.*/", "", paramFiles)

# replace the names by model
names(paramData) <- repNames

# remove the names object
rm(repNames)

# load the summary statistics file
sumStatData <- lapply(sumStatFiles, 
                      function(x) {
                        
                        read.table(x, sep = "\t", 
                                   header = TRUE)
                        
                      })

# get the names of the models
repNames <- gsub("^.*/", "", sumStatFiles)

# rename the summary statistics
names(sumStatData) <- repNames

# remove the model names object
rm(repNames)

# load the empirical data set
empData <- read.table(empFile, 
                      sep = "\t", 
                      header = TRUE)

# subset the parameters by the specified number of simulations
paramUse <- lapply(paramData, function(x) x[1:numSim, ])

# subset the summary statistics by the specified number of simulations
sumstatUse <- lapply(sumStatData, function(x) x[1:numSim, ])

# generate model label based on the name
modelLabel <- lapply(1:length(paramUse), 
                     function(x) {
                       
                       rep(names(paramUse)[x], 
                           nrow(paramUse[[x]]))
                       
                     })

# unlist the model indices
modelLabel <- unlist(modelLabel)

# make a data frame of the summary statistics
sumStat <- do.call("rbind", sumstatUse)

# get the indices of teh summary statistics to use
statUse <- grep(sumstatPat, names(empData))

#############################################
#### perform PCA                         ####
#############################################
# ideal to do this just to see if the empirical
# data falls within the cloud of simulated datasets
# NOTE: might be more ideal to simply include the empirical
#       data in generating the pca; in this example, only
#       the simulated data were included. The pc score of
#       the empirical data was determined afterwards.

# make an index for the color
colInd <- factor(modelLabel)

# make a pca
pcaRes <- prcomp(sumStat[, statUse], 
                 center = TRUE, 
                 scale. = TRUE)


# predict the empirical data based on the PCA of the
# other stuff
# OR an alternative is simply add the empirical data
# when generating the PCA
pcaEmp <- predict(pcaRes, empData[ ,statUse])

# make a device for the pca plot
jpeg(paste(pcaFile, ".jpg", sep = ""),
     height = 1800, width = 2400,
     res = 200, quality = 100)

# make a background plot
plot(x = pcaRes$x[ ,1],
     y = pcaRes$x[ ,2],
     col = c("black", "blue")[colInd],
     pch = 20,
     type = "n",
     xlab = "PC1",
     ylab = "PC2")

# shuffle the points
shuffle <- sample(1:nrow(pcaRes$x), 
                  nrow(pcaRes$x))

# plot the simulated data
points(x = pcaRes$x[shuffle ,1],
       y = pcaRes$x[shuffle ,2],
       col = c("black", "blue")[colInd][shuffle],
       cex = 0.75,
       pch = 20)

# plot the empirical data
points(x = pcaEmp[ ,1],
       y = pcaEmp[ ,2],
       col = "red",
       cex = 3,
       pch = 20)

# turn off the device
dev.off()

#### END SECTION


#############################################
#### estimate posterior support          ####
#############################################

# set seed for repeatability of the 
# posterior estimates of the empirical data
set.seed(666)

# perform abc anaylis with neural net method and 
# tolerance  of 0.01  
neuralPostpr <- abc::postpr(target = empData[statUse],
                            index = modelLabel,
                            sumstat = sumStat[statUse],
                            tol = 0.01,
                            method = "neuralnet",
                            maxit = 1500,
                            numnet = 20,
                            sizenet = 15,
                            MaxNWts = 2000)

# check results
summary(neuralPostpr)

# save the output as a compressed binary
saveRDS(neuralPostpr,
        file = postprNeural)

#############################################
#### check distribution of post samples  ####
#############################################
# there are instances of high support for a model
# if no weighing is applied, leading to a very large
# discrepancy between the posterior probability estimates of 
# neural network and raw proportion of posterior samples.
# Here, the weight of each posterior sample
# is examined. Note that the weight is proportional
# to the distance of a posterior sample from the empirical data 

# summary of simulated data set selected
postSamp <- as.character(neuralPostpr$values)

# remove other characters
postModel <- gsub("Param_amilFinal.*$", "", postSamp)

# generate a table summarizing model ID and weight
postWeight <- data.frame(model = postModel, 
                         weight = neuralPostpr$weights,
                         stringsAsFactors = FALSE)

# examine distributions of weights for the two models
boxplot(postWeight$weight ~ postWeight$model)

# make summary
summary(postWeight$weight[postWeight$model == "div"])
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000156 0.0172683 0.0392267 0.0517135 0.07323 91 0.3422224 

# make summary
summary(postWeight$weight[postWeight$model == "ibd"])
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001654 0.029186 0.057564 0.107234 0.191292 0.310447 

