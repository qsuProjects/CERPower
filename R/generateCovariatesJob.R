

# number_cores <- 16
################################################
# Set up parallel backend
################################################

require(doSNOW)
require(plyr)

cl<-makeCluster(number_cores)
registerDoSNOW(cl)

number_cores
getDoParWorkers()

clusterApply(cl, seq(along=cl), function(.core_num, .segment_index, .log_write_directory, 
                                         .covariates_name, .data_write_directory) {
  set.seed(.segment_index * 2 ^ .core_num)
  
  library(devtools)
  library(MBESS)
  library(ICC)
  library(miscTools)
  library(corpcor)
  library(psych)
  library(ggplot2)
  library(stringr)
  
  WORKER.ID <<- paste0("job_", str_pad(.segment_index, 4, side = "left", pad = "0"),
                       "_core_", str_pad(.core_num, 2, side = "left", pad = "0"))
  
  LOG.FILE.PATH <<- paste0(.log_write_directory, "/", WORKER.ID, ".txt")
  # MODEL.OUTPUT.FILE.PATH <<- paste0(.results_write_directory, "/", WORKER.ID, "_", .sim_results_name)
  
  COVARIATE.OUTPUT.FILE.PATH <<- paste0(.data_write_directory, .covariates_name)
  
  # cat("var,coef,se,p,type,approach,beta,prop_obs_dropped,prop_subjects_dropped,prop_censored,prop_final_NA,rep,cov_file\n",
  #     file = MODEL.OUTPUT.FILE.PATH, append = FALSE)
  
  cat(paste0("################################################\n",
             "Starting covariate generation...\n",
             "################################################\n"), 
      file = LOG.FILE.PATH, append = FALSE)
  
}, segment_index, log_write_directory, covariates_name, data_write_directory)


l_ply(range_low:range_high, .parallel = TRUE, function(.set_number,
                                                       .n_subjects, 
                                                       .n_observations,
                                                       .n_drugs,
                                                       .source_location) {
  
  
  ######################### SET SIMULATION PARAMETERS #########################
  # if need to install paxkages
  #install.packages("plyr")
  
  cat(paste0(Sys.time(),
             " Loading functions for dataset ", .set_number,
             "\n"), 
      file = LOG.FILE.PATH, append = TRUE)
  
  source(paste0(.source_location, "jointly_generate_binary_normal_modified_v2.R"))
  source(paste0(.source_location, "load_functions.R"))
  
  # source code from GitHub
  # source_url("https://raw.githubusercontent.com/qsuProjects/PCORI/master/r/GENCOV/load_functions.R")
  
  # source_url("https://raw.githubusercontent.com/qsuProjects/PCORI/master/r/GENCOV/jointly_generate_binary_normal_modified_v2.R", 
  # local = FALSE)
  
  #############################################
  ##        Set Simulation Parameters        ##
  #############################################
  
  # name prefix for all datasets
  # name_prefix = "ex2"
  if(FALSE){
    # n (number of subjects)
    .n_subjects = 1000
    
    # obs (number of observations per subject)
    .n_observations = 150
    
    # n.Reps (number of datasets to generate)
    # n.Reps = 100
    .set_number <- 5
    # n.Drugs (number of drug variables)
    .n_drugs = 15
    .source_location <- "~/shexport/PCORI/power/src/r/"
  }
  
  #############################################
  ## Read in Correlation and Parameter Files ##
  #############################################
  cat(paste0(Sys.time(),
            " Loading parameters for dataset ", .set_number,
            "\n"), 
      file = LOG.FILE.PATH, append = TRUE)
  # within-subject correlation matrix
  .wcor = read.csv(paste0(.source_location, "ex2_wcor_time_vary.csv"), header=FALSE)[-1,-1]
  
  # population correlation matrix
  .pcor = read.csv(paste0(.source_location, "ex2_pcor_time_vary.csv"), header=TRUE)[,-1]
  
  cat(paste0(Sys.time(),
             " Parameters for dataset ", .set_number, " loaded, generating parameter data frame",
             "\n"), 
      file = LOG.FILE.PATH, append = TRUE)
  # read in and complete parameters dataframe
  .parameters_csv <- read.csv(paste0(.source_location, "power_parameters.csv"))
  .parameters = complete_parameters(.parameters_csv, 
                                    .n_subjects, .set_number = .set_number, .log_path = LOG.FILE.PATH)
  
  cat(paste0(Sys.time(),
             " Parameters data frame generated for dataset ", .set_number,
             "\n"), 
      file = LOG.FILE.PATH, append = TRUE)
  #cat.parameters = read.csv("C:/Users/npurington/Documents/PCORI/Maria Paper/ex1_categorical_parameters.csv")
  
  
  #############################################
  ##               Simulate!                 ##
  #############################################
  
  ### following two parameters defined here because we can't change source code
  ### in program make_one_dataset change line 
  # bins.prop = proportionize(bins, zero, one)  # proportionized version of the binaries
  ### to
  # bins.prop = proportionize(bins)  # proportionized version of the binaries
  ### unless we want to change the defaults for variables zero and one (which are the values below)
  
  
  # ptm0 <- proc.time()
  
  cat(paste0(Sys.time(),
             " Running repeat_sim for dataset ", .set_number,
             "\n"), 
      file = LOG.FILE.PATH, append = TRUE)
  # simulate results
  .gen_time <- system.time({ 
    sim = repeat_sim(n=.n_subjects, obs = .n_observations, parameters = .parameters, prop.target=NULL,
                     mean.target = NULL, n.Drugs = .n_drugs, 
                     pcor = .pcor, wcor = .wcor, n.Reps = 1,
                     # write.data=TRUE,
                     #name_prefix= paste( .name_prefix, WORKER.ID, sep="_" ),  # used with Sherlock
                     # name_prefix=name_prefix,
                     cat.parameters=NULL,
                     zero=0.00001,
                     one=0.999, 
                     .set_number = .set_number, .log_path = LOG.FILE.PATH)
  })

  cat(paste0(Sys.time(),
             " Running repeat_sim for dataset ", .set_number,
             " took ", round(.gen_time[3]/60,2), " minutes\n"), 
      file = LOG.FILE.PATH, append = TRUE)
  
  cat(paste0(Sys.time(),
             " Writing covariates for dataset ", .set_number,
             "\n"), 
      file = LOG.FILE.PATH, append = TRUE)
  
  write.table(x = sim[[1]], file = paste0(COVARIATE.OUTPUT.FILE.PATH, str_pad(.set_number, width = 4, side = "left", pad = "0"), ".csv"), 
              sep = ",", col.names = TRUE, row.names = FALSE)
  
  cat(paste0(Sys.time(),
             " Dataset ", .set_number, " complete",
             "\n"), 
      file = LOG.FILE.PATH, append = TRUE)
  
}, n_subjects, n_observations, n_drugs, source_location)

# tdf <- sim[[1]]
# summary(tdf)
# sim[[1]][1,]
# sim[[1]][1,]
# names(sim)
# # extract dataset
# d = sim$data
# head(d)
# 
# proc.time() - ptm0
# #~2.7 days
# 
# ######################### QUICK TOUR THROUGH THE SIMULATED DATA ########################
# 
# ##### Example of Normal Variable Clustered within a Subject #####
# temp = d[ d$id %in% 1:10, ]  # look at first 6 subjects
# 
# ### create time variable; everyone got 150 records
# temp$t <- rep(c(1:150),10)
# 
# ggplot( data=temp, aes(x=t, y=bmi) ) + geom_line() + facet_grid(id~.) 
# 
# ggplot( data=temp, aes(x=t, y=log_vln) ) + geom_line() + facet_grid(id~.) 
# 
# ggplot( data=temp, aes(x=t, y=bps) ) + geom_line() + facet_grid(id~.) 
# 
# table(temp$d_abac) 
# 
# ##### Correlation Matrix Across Subjects #####
# # extract observed across-subjects correlation matrix
# sim$corr
# 
# # how biased were the correlations?
# sim$corr.bias
# 
# 
