# rm(list = ls())
# NOTE:
# This scripts links fastsimcoal and arlequin to simulate genetic data and summarize them in some statistics.
#   The script utilizes the .param file of fastsimcoal to pull out the parameters and the bounds of their distributions.
#   Using this information, random values are then drawn. These values are then substituted in the .tpl file of fastsimcoal.
#   After that, system calls for fastsimcoal and arlequin are made. This system calls can be edited to change the settings of
#   these programs (see fscCom and arlSumCom objects below).
#
# IMPORTANT: The filenames of the fastsimcoal and arelquin executables were hard-coded below;
#               see the fscDir and arlDir objects below to change them.
#            The arlequin command-line can be tricky. Make sure to select the proper summary statistics in
#               the ssdef.txt file, and make sure that the arl_run.ars is properly configured. The easiest
#               way to configure the arl_run.ars is to start an arp project in the arlequin gui, and configure
#               the project accordingly and save it. Then pull out the necessary files from it.
#               see instructions from arlequin http://cmpg.unibe.ch/software/arlequin35/
#
#
#### USAGE
# Rscripts fscArp_CL.R dirWork pathDir numSim modelSuf parRep fscCores
# dirWork:        directory containing the .tpl, .param, fastsimcoal executable, 
#                     arlequin executable and its associated file
# numSim:         number of simulations to perform
# modelSuf:       a label to prefix/suffix the outputs
# parRep:         number of replicates for each parameter combination (kept to 1)
# fscCores:       number of cores to pass to fastsimcoal

# Example:
# # Rscript fscArp_CL.R "./DIVrep1/" 50000 DIV 1 2


#### SOME FUNCTIONS

# function to draw parameters
drawParam <- function(parTemp, numRep){
  
  # get the values from the parameter specification
  parTemp <- unlist(strsplit(parTemp, " "))
  
  # get the flag for integer/float
  intID <- as.integer(parTemp[1])
  
  # define the bounds and assign an appropriate type
  # if 1 the bounds should be integer
  if(intID){
    
    # get the lower bound
    paramLB <- as.integer(parTemp[4])
    
    # get the upper bound
    paramUB <- as.integer(parTemp[5])
    
    # message("YES")
    
    # if 0 bounds should be float  
  }else{
    
    # get the lower bound
    paramLB <- as.double(parTemp[4])
    
    # get the upper bound
    paramUB <- as.double(parTemp[5])
    
    # message("NO")
  }
  
  # generate the parameters
  paramVal <- runif(numRep, paramLB, paramUB)
  
  # if the paramater is integer
  # round the values
  if(intID){
    
    paramVal <- round(paramVal)
    
    
  }
  
  # convert to a list
  paramVal <- list(paramVal)
  
  # use the paramter name as name of the list
  names(paramVal) <- parTemp[2]
  
  # return the list of 1 parameter
  return(paramVal)
  
}

#### END FUNCTIONS

##########################
##### INPUTS SECTION #####
##########################

# capture input
inputs <- commandArgs(trailingOnly = TRUE)

st <- Sys.time()
# number of parameter combinations to draw
# this is the basic number of simulations to perform
numSim = as.integer(inputs[2])

modelSuf = as.character(inputs[3])

# number of replication for each parameter combination
# this is to capture the variance for each parameter combination
# the effective number of simulations is equalt to numSim * parRep
parRep = as.integer(inputs[5])

# number of cores to utilize for the fastsimcoal
fscCores = as.integer(inputs[4])

# working directory
# this directory should include the following:
#   1) .est file which contains the parameter prior distribution
#   2) .tpl file which contains the model
#   3) fsc26 simulation program
#   4) arlsumstat and its associated files
dirWork = inputs[1]

# directory of the fastsimcoal
fscDir = paste(dirWork, "/fsc26", sep = "")

# directory of the arlsumstat
arlDir = paste(dirWork, "/arlsumstat3522_64bit", sep = "")

# directory of the summary statistics table
# this is the file where the sumstat will be saved
ssOut = paste(dirWork, "/ss", 
              modelSuf, sep = "")

# directroy of paramaters table
# this is the file where the paramaters will be saved
paramOut = paste(dirWork, "/param", 
                 modelSuf, sep = "")


##########################
##### INPUTS END     #####
##########################

# generate a template file for the
# summary statistics
sstempDir <- paste(dirWork, "/ssTemp", sep = "")

# set the working directroy to the dirWork
setwd(dirWork)

# find the .est file
estName <- list.files(dirWork, 
                      pattern = ".est$", 
                      full.names = TRUE)

# load the .est file
estFile <- readLines(estName)

# remove the comments
estFile <- estFile[-grep("//", estFile)]

# remove blank lines
estFile <- estFile[-grep("^$", estFile)]

# replace tabs with space
estFile <- gsub("\t", " ", estFile)

# remove multiple spaces
# https://stackoverflow.com/questions/25707647/merge-multiple-spaces-to-single-space-remove-trailing-leading-spaces
estFile <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", 
                "", 
                estFile, 
                perl=TRUE)

# extract the parameter lines
paramLines <- estFile[(grep("\\[PARAMETERS\\]", estFile) + 1):
                        (grep("\\[RULES\\]", estFile) - 1) ]


# extract rules line

if((grep("\\[COMPLEX PARAMETERS\\]", estFile) - 1) >= (grep("\\[RULES\\]", estFile) + 1)){
  
  rulesLines <- estFile[(grep("\\[RULES\\]", estFile) + 1):
                          (grep("\\[COMPLEX PARAMETERS\\]", estFile) - 1)]
  
  
}else {
  
  rulesLines <- estFile[0]
  
}

# draw values for the set of parameters using drawParam function
tempParam <- unlist(lapply(paramLines, 
                           function(x) 
                             drawParam(parTemp = x,
                                       numRep = numSim)),
                    recursive = FALSE)

if(length(rulesLines) >= 1){
  
  for(rule in rulesLines){
    
    # extract the rules
    ruleTemp <- unlist(strsplit(rule, " "))
    
    # get the the line of B
    parTemp <- unlist(strsplit(paramLines[grep(ruleTemp[3], 
                                               paramLines)], " "))
    
    # get the flag for integer/float
    intID <- as.integer(parTemp[1])
    
    if(ruleTemp[2] == "GREATER"){
      
      # if rule is A > B and it is false (B > A)
      # we would want to make B smaller
      # to do that, we resample B by setting the lower bound to 
      # its initial lower bound and its higher bound by A - 1
      
      # find parameters to modify 
      redrawParam <- which(!(tempParam[[ruleTemp[1]]] > 
                               tempParam[[ruleTemp[3]]]))

      # define the bounds and assign an appropriate type
      # if 1 the bounds should be integer
      if(intID){
        
        # get the upper bound
        paramLB <- as.integer(parTemp[4])
        
      # if 0 bounds should be float  
      }else{
        
        # get the upper bound
        paramLB <- as.double(parTemp[4])
        
      }
      
      # resample the parameters by applying the above rules
      paramVal <- mapply(LB = rep(paramLB, length(redrawParam)), 
                         UB = (tempParam[[ruleTemp[1]]][redrawParam] - 1), 
                         function(LB, UB){
                           
                           return(runif(1, LB, UB))
                           
                           })
      
      # define the bounds and assign an appropriate type
      # if 1 the bounds should be integer
      if(intID){
        
        # get the upper bound
        paramVal <- as.integer(paramVal)
        
        # message("YES")
        
        # if 0 bounds should be float  
      }else{
        
        # get the upper bound
        paramVal <- as.double(paramVal)
        
        # message("NO")
      }
      
      # replace the parameter values
      tempParam[[ruleTemp[3]]][redrawParam] <- paramVal

    }else if(ruleTemp[2] == "LESS"){
      
      # if rule is A < B and it is false (B < A)
      # we would want to make B larger
      # to do that, we resample B by setting the lower bound to 
      # A + 1 and its  its higher bound to its original higher bound
      
      # find parameters to modify 
      redrawParam <- which(!(tempParam[[ruleTemp[1]]] < tempParam[[ruleTemp[3]]]))
      
      # define the bounds and assign an appropriate type
      # if 1 the bounds should be integer
      if(intID){
        
        # get the upper bound
        paramUB <- as.integer(parTemp[5])
        
        # if 0 bounds should be float  
      }else{
        
        # get the upper bound
        paramUB <- as.double(parTemp[5])
        
      }
      
      # resample the parameters by applying the above rules
      paramVal <- mapply(LB = (tempParam[[ruleTemp[1]]][redrawParam] + 1), 
                         UB = rep(paramUB, length(redrawParam)), 
                         function(LB, UB){
                           
                           return(runif(1, LB, UB))
                           
                         })
      
      # define the bounds and assign an appropriate type
      # if 1 the bounds should be integer
      if(intID){
        
        # get the upper bound
        paramVal <- as.integer(paramVal)
        
      # if 0 bounds should be float  
      }else{
        
        # get the upper bound
        paramVal <- as.double(paramVal)
        
      }
      
      # replace the parameter values
      tempParam[[ruleTemp[3]]][redrawParam] <- paramVal
      
    }else{
      
      # should just be GREATER OR LESS
      stop("Check the rules!\nIt can only be GREATER or LESS")
      
    }

  }

}

# convert the result to a data.frame
tempParam <- as.data.frame(tempParam, stringsAsFactors = FALSE)

#### CHANGE PARAM

# pull out the template file
tplName <- list.files(dirWork, 
                      pattern = ".tpl$", 
                      full.names = TRUE)

# read the template file
tplFile <- readLines(tplName)

# list the parameters
paramNames <- names(tempParam)

# start the simulation
for(i in 1:numSim){
  
  # i = 1
  # copy the template to a parameter file
  parFile <- tplFile
  
  # replace the variable with the parameter draw
  for(param in paramNames){
    
    # replace the values
    parFile <- gsub(param, tempParam[i, param], parFile)
    
  }
  
  # create a file name for the parameter file
  parName <-  gsub(".tpl", ".par", tplName)
  
  # write the parameter file
  write(parFile, file = parName)
  
  # create a fsc26 call
  fscCom = paste(fscDir, "-i", parName, "-n", parRep, "-c", fscCores, "-g 1")
  
  # simulate using the parameter file
  system(fscCom)
  
  # get the directory of arp files
  resPath <- gsub(".tpl", "", tplName)
  
  # create a call to move the arp files
  mvCom <- paste("mv ", resPath, "/*.arp", " ", dirWork, sep = "")
  
  # move the arp file to the working directory
  system(mvCom)
  
  # list the arpfiles
  arpFiles <- list.files(dirWork, 
                         pattern = ".arp$", 
                         full.names = TRUE)
  
  # if there is replication per parameter files
  # loop over arp files
  if(length(arpFiles) > 1){
    
    # the number of arp files should correspond to parRep
    for(arpNum in 1:parRep){
      
      # if this is the first file
      # generate a header
      if(arpNum == 1){
        
        # generate call for arlsumstat
        arlsumCom <- paste(arlDir, arpFiles[arpNum], sstempDir, "0", "1")
        
        # if not the first file
        # just append
      }else{
        
        # generate call for arlsumstat
        arlsumCom <- paste(arlDir, arpFiles[arpNum], sstempDir, "1", "0")
        
      }
      
      # calculate summary statistics
      system(arlsumCom)
      
    }
    
    # if there is no replication, there should only be one arp file  
  }else{
    
    # generate call for arlsumstat
    arlsumCom <- paste(arlDir, arpFiles, sstempDir, "0", "1")
    
    # calculate summary statistics
    system(arlsumCom)
    
  }
  
  # this bit is to save the outputs
  # if this is the first simulation
  if(i == 1){
    
    # copy all sumstat to ssOut with header
    system(paste("cat", sstempDir, ">>", ssOut))
    
    # write the parameters with header
    write.table(tempParam[rep(i, parRep), ], 
                file = paramOut,
                sep = "\t",
                col.names = TRUE, row.names = FALSE)
    
  }else{
    
    # copy all sumstat to ssOut without header
    system(paste("tail -n +2", sstempDir,  ">>", ssOut))
    
    # write the parameters without header
    write.table(tempParam[rep(i, parRep), ], 
                file = paramOut, 
                sep = "\t",
                col.names = FALSE, row.names = FALSE,
                append = TRUE)
    
    
  }
  
  # remove all arp files
  system(paste("rm ", dirWork, "/*.arp", sep = ""))
  
  # remove the sstempOut file
  system(paste("rm ", sstempDir, sep = ""))
  
  # remove the fsc result directory
  system(paste("rm -r ", resPath, sep = ""))
  
  # remove arlsum folders
  system("rm -r *.res")
  
}

et <- Sys.time()

runTime <- et - st

write(print(runTime), file = "runTime")

# write(attr(runTime, "units"), file = "runTime", append = TRUE)

write(print(st), file = "runTime", append = TRUE)

write(print(et), file = "runTime", append = TRUE)
