# The approximate bayesian analysis for the two acropora species is divided into the following four components:
#	1) simulation of genetic data (generation of parameters and summary statistics)
#	2) estimation of posterior probability of the models
#	3) simulation of Pseudo-Observed-Data (POD) and estimating their posterior probabilities
#   4) examination of the distribution of the PODs posterior probabilties to determine
#			robustness

# These components are implemented in the scripts contained in this directory, namely:
#	1) fscArp_CL.R
#	2) estimatePostProb.R
#	3) crossValid_CL.R
#	4) examineRobust.R
#
#	NOTE: All of these scripts and their required programs were ran under a unix system (ubuntu).
#
# CONTENTS:
# Included in this folder are the following:
#		a) four Rscripts (described below)
#		b) an example working directory (inlcuding all necessary files) for fscArp_CL.R (DIVrep1);
#		c) an example working directory for running the crossValid_CL.R (amilSample), which contains the following:
#				- summary statistics table for two models (divSumstat_amilFinal and ibdSumstat_amilFinal) 
#					-- these files are needed by estimatePostProb.R and crossValid_CL.R;
#				- paramter table for the two models (divParam_amilFinal and ibdParam_amilFinal) -- not really needed, just model indices;
#				- directories for the pseudo-observed data for two models (amil_DIVpods and amil_IBDpods), each of which
#						contains a paramter and summary statistics table of the pods -- this is used by crossValid_CL.R;
#				- a directory containing the sample output of crossValid_CL.R -- the content of this directory is utilized by examineRobust.R;
#		d) an empirical data file (amilEmpirical), which is needed in estimatePostProb.R;
#
# IMPORTANT: These files are just examples and are subsets of the analysed data;
#
# BRIEF OVERVIEW:
# fscArp_CL.R utilizes fastsimcoal (http://cmpg.unibe.ch/software/fastsimcoal2/) and arlequin (http://cmpg.unibe.ch/software/arlequin35/)
#		to generate genetic data and compute summary statistics from it. It utilizes the .tpl and .param files of the fastsimcoal to
#		draw parameter values (which here is limited to uniform distribution). The fastsimcoal can do the generation of parameter values,
#		but it was implemented in the fscArp_CL.R to continuously feed the outputs to arlequin (NOTE: named pipes might
#		an alternative for this). Thus, models and parameter values will still be defined in the .tpl and .param files -- the construction of
#		which is will documented in the fastsimcoal manual. NOTE: setting up the command-line version of the arlequin can be tricky, and
#		requires careful reading of its manual (keep in mind that the arl_run.ars and ssdefs.txt are critical components of 
#		the arlequin command-line). This script is designed for command-line use.
#
#		Simulation of genetic data involving many individuals takes a lot of time; thus, running parallel instances of 
#		fasArp_CL.R is needed. The easiest way to do this is to place all the needed files in a directory, including the 
#		fastsimcoal and arlequin (and their associated file; see the contents of the DIVrep1 example). Making multiple copies
#		of that directory and running parallel fscArp_CL.R on them will hasten the simulations.
#		Example implentation:
#		# copy the whole directory	
#		cp -r DIVrep1 DIVrep2 
#
# 		# run the simulation in background (&) and disown (see fscArp_CL.R for details)
#		Rscript fscArp_CL.R /mnt/dinFiles/2018_CoralABC/sample/DIVrep2/ 100 DIV 1 2 & disown
#
# estimatePostProb.R performs the estimation of posterior probabilites for each model. This requires the summary statistics of the
#		empirical data set (which should be calculated in the same way as the simulated data set). If parallel processes were performed
#		to simulate genetic data, it will be ideal to merge all of them in a single parameter and summary statistics file for each model.
#		This script work with such structure: one paramter file and one summary statistics file for each model, and one empirical summary statistics
#		file (unless multiple empirical data points are availble). See estimatePostProb.R for the details. This is script is better to be used
#		in an R GUI.
#		
# crossValid_CL.R is the initial step in performing the cross-validation procedure. The main goal of the cross-validate is to give confidence on the
#		estimate of posterior probability of the supported model (i.e., somehow test that a model that is supported is not due to random chance). Performing
#		the cross-validation requires simulation of PODs, which can be done using the fscArp_CL.R. This script is designed for command-line use
#		For example, simulating 1000 PODs for the DIV model can be performed with this command:
#		
#		# use the same files 
#		Rscript fscArp_CL.R /mnt/dinFiles/2018_CoralABC/sample/DIVrep1/ 1000 DIV 1 2 & disown 
#		
#		# NOTE: the PODs should be simulated using the same approach as the simulated-data sets being utilized. In the R abc package, the cross-validation actually
#		# makes use of a certain proportion of the simulated genetic data as PODs
#
#		After simulating the PODs, simply get the output files and use it for the crossValid_CL.R, which then performs a similar analysis implemented in
#		estimatePostProb.R -- that is estimating the posterior probability of each POD. Similar to simulating genetic data set, estimating the posterior probability
#		can take some time, which prevents performance of cross-validation with high number of PODs. This is the rationale behind the crossValid_CL.R, which allows
#		parallel processes of the PODs (well, it simply splits the summary statistics and perform posterior probability estimation on those PODs).
#		See crossValid_CL.R for details.
#
# examineRobust.R gathers the results of the crossValid_CL.R (which are binary rds files) and identifies a threshold for posterior probability by which a support for
#		 model has a certain level of confidence X (which in this case is 95%). Please see Roux et al 2016 PloS Biology 14 e2000234 for details about the robustness, which 
#		 is the basis of the examineRobust.R (see the script for futher details). This is script is better to be used
#		 in an R GUI
#

			
			
