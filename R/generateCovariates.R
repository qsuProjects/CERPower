# .source_location <- source_location <- "~/shexport/PCORI/power/src/r/"
# LOG.FILE.PATH <- "~/shexport/PCORI/power/debugLog.txt"
# number_cores <- 16
# n_data_sets <- 100
# n_jobs <- 7
# segment_index <- 1
# covariates_name <- "powerCovariates_"
# log_write_directory <- "/scratch/PI/manishad/PCORI/power/log/"
# data_write_directory <- "/scratch/PI/manishad/PCORI/power/covariates/"

# scp ~/shexport/PCORI/bigSim/r/runImposeImpute.R kikapp@sherlock:~/PCORI/lib/bigSim/r/
args <- commandArgs(trailingOnly = T)
print(args)

source_location <- args[1]
number_cores <- as.numeric(args[2])
n_data_sets <-  as.numeric(args[3])
n_jobs <- as.numeric(args[4])
segment_index <- as.numeric(args[5])
covariates_name <- args[6]
log_write_directory <- args[7]
data_write_directory <- args[8]
n_subjects <- as.numeric(args[9] )
n_observations <- as.numeric(args[10])
n_drugs <- as.numeric(args[11])

segment_endpoints <- ceiling(quantile(1:n_data_sets, 
                                      seq(from = 0, to = 1, length.out = (n_jobs + 1)) ))


#figure out which data sets to generate
range_low <- segment_endpoints[segment_index]
range_high <- segment_endpoints[segment_index + 1]
range_low
range_high
# load impose missingness code
source(paste0(source_location, "generateCovariatesJob.R") )
