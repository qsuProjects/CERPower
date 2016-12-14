  source("/Volumes/QSU/Datasets/sherlock/r/generateSbatchPSwitch.r")
  library(stringr)
  # list.files("/Volumes/QSU-1")
  r_version <- "3.0.2"
  n_jobs <- 120
  jobs <- expand.grid(segnum = 1:n_jobs,
                                    jobtime = "8:00:00",
                                    quality = "normal",
                                    node_number = 1,
                                    mem_per_node = 64000,
                                    mailtype =  "ALL",
                                    user_email = "kikapp@stanford.edu",
                                    tasks_per_node = 1,
                                    cpus_per_task = 6,
                                    path_to_r_script = "/share/PI/manishad/PCORI/power/src/r/generateSurvivalTimes.R",
                                    stringsAsFactors = F)
  
  source_location <- paste0("/share/PI/manishad/PCORI/power/src/r/generateSurvivalTimesJob.R")
  beta_location <- paste0("/share/PI/manishad/PCORI/power/src/r/betasPower.R")
  covariate_data_path <- paste0("/scratch/PI/manishad/PCORI/power/covariates/")
  n_cores <- 6
  n_segments <- 1
  segment_index <- 1
  results_write_directory <- paste0("/scratch/PI/manishad/PCORI/power/modelFit")
  sim_results_name <- "survtime.csv"
  log_write_directory <- paste0("/scratch/PI/manishad/PCORI/power/logSurvival/")
  data_write_directory <- paste0("/scratch/PI/manishad/PCORI/power/covariatesSurvival/")
  # beta_set <- 1:120
  
  jobs$jobname <- paste0("pgs", str_pad(jobs$segnum, 4, side = "left", pad = "0"))
  jobs$outfile <- paste0("opgs", str_pad(jobs$segnum, 4, side = "left", pad = "0"), ".out")
  jobs$errorfile <- paste0("epgs", str_pad(jobs$segnum, 4, side = "left", pad = "0"), ".err")
  jobs$r_version <- paste0("/", r_version)
  
  jobs$args_to_r_script <- paste("--args ", source_location, beta_location, covariate_data_path,
                                 n_cores, n_segments, segment_index, results_write_directory,
                                 sim_results_name, log_write_directory, data_write_directory,
                                 jobs$segnum, sep = " ")
  
  jobs$write_path <- paste0("~/shexport/PCORI/power/src/sbatchSurvival/", 
                            str_pad(jobs$segnum, 4, side = "left", pad = "0"), 
                            "_powerSurv.sbatch")
  jobs$server_sbatch_path =  paste0("/share/PI/manishad/PCORI/power/src/sbatchSurvival/", 
                                    str_pad(jobs$segnum, 4, side = "left", pad = "0"), 
                                    "_powerSurv.sbatch")
  
  runfile_path <-  paste0("~/shexport/PCORI/power/src/r/runPowerSurv.R")
  files <- generateSbatch(jobs, runfile_path = runfile_path, run_now = F)
  
  
  #transfer files to sherlock (must have kinit credentials)
  system("scp ~/shexport/PCORI/power/src/sbatchSurvival/*.sbatch kikapp@sherlock:/share/PI/manishad/PCORI/power/src/sbatchSurvival/")
  system("scp ~/shexport/PCORI/power/src/r/*.R kikapp@sherlock:/share/PI/manishad/PCORI/power/src/r/")
  
