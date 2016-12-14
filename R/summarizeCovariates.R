

number_cores <- 16
data_read_directory <- "/scratch/PI/manishad/PCORI/power/covariates/"
data_write_directory <- "/scratch/PI/manishad/PCORI/power/covariatesSummary/"
covariates_name <- "covariateSummary_"
################################################
# Set up parallel backend
################################################

require(doSNOW)
require(plyr)

cl<-makeCluster(number_cores)
registerDoSNOW(cl)

number_cores
getDoParWorkers()

clusterApply(cl, seq(along=cl), function(.core_num, .data_write_directory, .covariates_name) {
  
  library(devtools)
  library(MBESS)
  library(ICC)
  library(miscTools)
  library(corpcor)
  library(psych)
  library(ggplot2)
  library(stringr)
  
  WORKER.ID <<- paste0("core_", str_pad(.core_num, 2, side = "left", pad = "0"))
  
}, data_write_directory, covariates_name)


l_ply(list.files(data_read_directory), .parallel = TRUE, function(.file_name, .data_read_directory,
                                                                  .data_write_directory) {
  
  .covariates <- read.table(file = paste0(.data_read_directory, "/", .file_name), 
                            sep = ",", header = T, stringsAsFactors = F)
  
  .covariates_sum <- ddply(.covariates, .(id), function(..subject_data) {
    ldply(names(..subject_data), function(...var, ...covariate_subject) {
      ...var_data <- ...covariate_subject[[...var]]
      data.frame(var = ...var,
                 mean = mean(...var_data, na.rm = TRUE),
                 sd = sd(...var_data, na.rm = TRUE),
                 n_NA = sum(is.na(...var_data)),
                 n = length(...var_data))
      
    }, ..subject_data )
  })
  
  write.table(x = .covariates_sum, file = paste0(data_write_directory, "summary_", .file_name), 
              sep = ",", col.names = TRUE, row.names = FALSE)
  
  
  .correlations <- data.frame(cor(.covariates))
  .correlations$var <- row.names(.correlations)
  write.table(x = .correlations, file = paste0(data_write_directory, "correlation_", .file_name), 
              sep = ",", col.names = TRUE, row.names = FALSE)
  
}, data_read_directory, data_write_directory)

