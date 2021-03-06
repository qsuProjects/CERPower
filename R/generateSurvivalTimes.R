# scp ~/shexport/PCORI/generateByScenario.R kikapp@sherlock:~/PCORI/lib
args <- commandArgs(trailingOnly = T)
# args <- vector("character", length = 11)
# print(args)
# args[1] <- "/share/PI/manishad/PCORI/power/src/r/generateSurvivalTimesJob.R"
# args[2] <- "/share/PI/manishad/PCORI/power/src/r/betasPower.R"
# args[3] <- "/scratch/PI/manishad/PCORI/power/covariates/"
# args[4] <- "8"
# args[5] <- "1"
# args[6] <- "1"
# args[7] <- "/scratch/PI/manishad/PCORI/power/modelFit"
# args[8] <- "survtime.csv"
# args[9] <- "/scratch/PI/manishad/PCORI/power/logSurvival/"
# args[10] <- "/scratch/PI/manishad/PCORI/power/covariatesSurvival/"
# args[11] <- "120"


print(args)



source_location <- args[1]
beta_location <- args[2]
covariate_data_path <- args[3]
number_cores <- as.numeric(args[4])
n_segments <- as.numeric(args[5])
segment_index <- as.numeric(args[6])
results_write_directory <- args[7]
sim_results_name <- args[8]
log_write_directory <- args[9]
data_write_directory <- args[10]
beta_set <- as.numeric(args[11])

#load betas
source(beta_location)
print(betas)

#determine which files to load
n_files <- length(list.files(covariate_data_path))
segment_endpoints <- round(quantile(1:n_files, probs = seq(0,1, length.out = (1 + n_segments))),0)
range_low <- segment_endpoints[segment_index]
range_high <- segment_endpoints[segment_index + 1]
cat(range_low, "\n")
cat(range_high, "\n")
# 
# segment_endpoints <- round(quantile(1:200, seq(0,1,1/20)))
# inde <- 20
# segment_endpoints[inde]
# segment_endpoints[inde+1]
# 
# 
# ?quantile

#initialize seed
the_seed <- segment_index

source(source_location)



#VVVV garbage VVVV#

# 
# source_location <- paste0("/share/PI/manishad/PCORI/power/src/r/survivalTimeGenerationAllFrailty.R")
# beta_location <- paste0("/share/PI/manishad/PCORI/power/src/r/betas.R")
# covariate_data_path <- paste0("/scratch/PI/manishad/PCORI/power/src/covariates/")
# n_cores <- 16
# n_segments <- 200
# segment_index <- 1
# results_write_directory <- paste0("/scratch/PI/manishad/PCORI/genSurv/output/genSurv/complete")
# sim_results_name <- "survtime.csv"
# log_write_directory <- paste0("/scratch/PI/manishad/PCORI/genSurv/log/genSurv/")
# data_write_directory <- paste0("/scratch/PI/manishad/PCORI/genSurv/covariates/timeCovariates/")
# 
# 
# source("~/shexport/PCORI/power/src/r/generateSurvivalTimesJob.R")
# source_location <- "~/shexport/PCORI/power/src/r/generateSurvivalTimesJob.R"
# beta_location <- "~/shexport/PCORI/power/src/r/betasPower.R"
# covariate_data_path <- "~/shexport/PCORI/power/covariates/"
# number_cores <- as.numeric("1")
# range_low <- as.numeric("1")
# range_high <- as.numeric("1")
# segment_index <- 1
# results_write_directory <- "~/shexport/PCORI/power/output/covariatesSurvival/"
# sim_results_name <- "testResults.csv"
# log_write_directory <- "~/shexport/PCORI/power/output/logSurvival/"

# source_location <- "/share/PI/manishad/PCORI/genSurv/r/generateTimeAndFit.R"
# beta_location <- "/share/PI/manishad/PCORI/power/src/r/betasPower.R"
# covariate_data_path <- "/scratch/PI/manishad/PCORI/genSurv/covariates/data_1000.csv"
# number_cores <- as.numeric("2")
# range_low <- as.numeric("1")
# range_high <- as.numeric("1000")
# results_write_directory <- "/scratch/PI/manishad/PCORI/genSurv/output/test"
# sim_results_name <- "testResults.csv"
# log_write_directory <- "/scratch/PI/manishad/PCORI/genSurv/log/test"
# data_write_directory <- "/scratch/PI/manishad/PCORI/genSurv/data/test"

# 
# if(is.na(output_location)) {output_location <- "/scratch/users/kikapp/dump"}
# if(is.na(n_cores)) {n_cores <- 16}
# if(is.na(range_low)) {range_low <- 0}
# if(is.na(range_high)) {range_high <- 5}
# if(is.na(source_location)) {source_location <- "~/PCORI/lib"}
# if(is.na(scenario_location)) {scenario_location <- "~/PCORI/lib"}
# if(is.na(scenario_row)) {scenario_row <- 1}
# 
# source_location <- "~/shexport/PCORI/genSurv/r/generateTimeAndFit.R"
# beta_location <- paste0("/share/PI/manishad/PCORI/genSurv/r/betas.R")
# covariate_data_path <-  paste0("/share/PI/manishad/sim_23_sher/")
# number_cores <- 8
# n_segments <- 50
# segment_index <- 28
# results_write_directory <- args[7]
# sim_results_name <- args[8]
# log_write_directory <- args[9]
# data_write_directory <- args[10]
