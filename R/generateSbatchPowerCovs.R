source("/Volumes/QSU/Datasets/sherlock/r/generateSbatchPSwitch.r")
library(stringr)

n_jobs <- 7
jobs <- expand.grid(segnum = 1:n_jobs,
                    jobtime = "5:00:00",
                    quality = "normal",
                    node_number = 1,
                    mem_per_node = 64000,
                    mailtype =  "ALL",
                    user_email = "kikapp@stanford.edu",
                    tasks_per_node = 1,
                    cpus_per_task = 12,
                    path_to_r_script = "/share/PI/manishad/PCORI/power/src/r/generateCovariates.R",
                    stringsAsFactors = F)

source_location <- "/share/PI/manishad/PCORI/power/src/r/"
number_cores <- 12
n_data_sets <-  100
n_jobs <- 7
segment_index <- jobs$segnum
covariates_name <- "powerCovariates_"
log_write_directory <- "/scratch/PI/manishad/PCORI/power/log/"
data_write_directory <- "/scratch/PI/manishad/PCORI/power/covariates/"
n_subjects <- 20000 
n_observations <- 150
n_drugs <- 15

jobs$jobname <- paste0("pow", str_pad(jobs$segnum, 4, side = "left", pad = "0"))
jobs$outfile <- paste0("opow", str_pad(jobs$segnum, 4, side = "left", pad = "0"), ".out")
jobs$errorfile <- paste0("epow", str_pad(jobs$segnum, 4, side = "left", pad = "0"), ".err")

jobs$args_to_r_script <- paste("--args ", source_location, number_cores,
                               n_data_sets, n_jobs, segment_index, covariates_name,
                               log_write_directory, data_write_directory, n_subjects,
                               n_observations, n_drugs, sep = " ")

jobs$write_path <- paste0("~/shexport/PCORI/power/src/sbatch/", 
                          str_pad(jobs$segnum, 4, side = "left", pad = "0"), 
                          "_power.sbatch")
jobs$server_sbatch_path =  paste0("/share/PI/manishad/PCORI/power/src/sbatch/", 
                                  str_pad(jobs$segnum, 4, side = "left", pad = "0"), 
                                  "_power.sbatch")

runfile_path <-  paste0("~/shexport/PCORI/power/src/r/runPowerCovs.R")
files <- generateSbatch(jobs, runfile_path = runfile_path, run_now = FALSE)


#transfer files to sherlock (must have kinit credentials)
system("scp ~/shexport/PCORI/power/src/sbatch/*.sbatch kikapp@sherlock:/share/PI/manishad/PCORI/power/src/sbatch/")
system("scp ~/shexport/PCORI/power/src/r/*.R kikapp@sherlock:/share/PI/manishad/PCORI/power/src/r/")
system("scp ~/shexport/PCORI/power/src/r/*.csv kikapp@sherlock:/share/PI/manishad/PCORI/power/src/r/")
# system("scp ~/shexport/PCORI/power/src/r/*.txt kikapp@sherlock:/share/PI/manishad/PCORI/power/src/r/")
