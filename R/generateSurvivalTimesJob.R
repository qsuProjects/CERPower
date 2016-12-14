if(FALSE) {
  rm(list = ls())
  source("~/.Rprofile")
}


###########################################
#this script loads simulated covariate data,
#uses those data to generate survival times,
#merges those survival times with covariate data,
#then fits various models on results

#It is written to be sourced from a script 
##executed on a multi-core system.

#Sourcing script should initialize:

#number_cores: Integer, number of cores to use
#the_seed: Integer, the seed
#covariate_data_path: String, path and filename for covariate data
#results_write_directory: String, where to write results
#beta: Data frame, one row, each column has a corresponding column
##in the covariate file. Names column names should be identical.
##values are the values to use to generate linear predictor
##for use in survival time generation algorithm
#range_low: Integer, first rep to do
#range_high: Integer, last rep to do
#sim_results_name: String, file name for simulation results
#log_write_directory: where to write log files
#data_write_directory: where to write data files


require(msm)
require(survival)
require(plyr)


#set up parallelization
require(doSNOW)
cl<-makeCluster(number_cores)
registerDoSNOW(cl)

number_cores
getDoParWorkers()

clusterApply(cl, seq(along=cl), function(.id, .results_write_directory, 
                                         .log_write_directory, .sim_results_name, 
                                         .segment_index, .beta_set) {
  
  
  set.seed(.segment_index*(2^.id))
  
  #load required packages
  
  require(msm)
  require(survival)
  require(plyr)
  require(coxme)
  require(geepack)
  require(gee)
  require(stringr)
  require(PermAlgo)
  require(combinat)
  
  WORKER.ID <<- paste0("job_", str_pad(.segment_index, 4, side = "left", pad = "0"),
                       "_beta_set_", str_pad(.beta_set, 4, side = "left", pad = "0"),
                       "_core_", str_pad(.id, 2, side = "left", pad = "0"))
  
  # }, results_write_directory, log_write_directory, sim_results_name, segment_index, beta_set)
  
  LOG.FILE.PATH <<- paste0(.log_write_directory, "/", WORKER.ID, ".txt")
  
  MODEL.OUTPUT.FILE.PATH <<- paste0(.results_write_directory, "/", WORKER.ID, "_", .sim_results_name)
  
  cat("var,coef,se,p,type,proportion_censored,rep,beta\n", file = MODEL.OUTPUT.FILE.PATH, append = F)
  
  cat(paste0("################################################\n",
             "Starting simulation...\n",
             "################################################\n"), file = LOG.FILE.PATH, append = F)
  
}, results_write_directory, log_write_directory, sim_results_name, segment_index, beta_set)
# list.files(log_write_directory)
# cat("test", file = paste0(log_write_directory, "test.csv"))
##main simulation loop
l_ply(range_low:range_high, .parallel = T, function(.rep, .covariate_data_path, .beta, .results_write_directory, 
                                                    .sim_results_name, .log_write_directory, .data_write_directory,
                                                    .beta_set) {
  # .rep <- 1
  # .covariate_data_path <- covariate_data_path
  # .beta <- betas
  # .results_write_directory <- results_write_directory
  # .sim_results_name <- "test" 
  # .log_write_directory <- log_write_directory
  # .data_write_directory <- data_write_directory
  # .beta_set <- 120
  # .coxme_object <- coxme(.coxme_formula, data = .data_none)
  # class(bdsmatrix::diag(.coxme_object$var))
  # as.vector(.coxme_object$var)
  
  # see coxme:::print.coxme
  coxmeCoefs <- function (.coxme_object) 
  {
    .tmp <- ""
    .beta <- .coxme_object$coefficients
    .nvar <- length(.beta)
    .nfrail <- nrow(.coxme_object$var) - .nvar
    .omit <- .coxme_object$na.action
    if (.nvar > 0) {
      .se <- sqrt(bdsmatrix::diag(.coxme_object$var)[.nfrail + 1:.nvar])
      .tmp <- cbind(.beta, 
                    .se, 
                    signif(1 - pchisq((.beta/.se)^2, 1), 2))
      dimnames(.tmp) <- list(names(.beta), c("coef", 
                                             "se(coef)", 
                                             "p"))
    }
    return(.tmp)
  }
  # .rep <- 1
  .rep_str <- str_pad(.rep, 4, pad = "0")
  # .covariate_data_path <- paste0("~/shexport/PCORI/power/covariates/")
  # .covariate_data_path <- covariate_data_path
  # .current_data_file <- "2015-04-09_job_580_worker_1_dataset_1"
  # .beta <- betas
  # .rep <- (range_low:range_high)[14]
  # .rep <- 1
  
  #load data
  .current_data_file <- list.files(.covariate_data_path)[.rep]
  cat(paste0("Reading covariates from ", .covariate_data_path, "/", .current_data_file, " JOB: ", WORKER.ID,
             " ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
  
  .covariates <- read.table(file = paste0(.covariate_data_path, "/", .current_data_file), 
                            sep = ",", header = T, stringsAsFactors = F)
  
  cat(paste0("Covariates read from ", .current_data_file, " JOB: ", WORKER.ID,
             " ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
  
  if(.beta_set == 1) {
    cat(paste0("Summarizing covariates for", " JOB: ", WORKER.ID,
               " ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    
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
    
    cat(paste0("Covariates summarized for", " JOB: ", WORKER.ID,
               " ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    
    cat(paste0("Writing covariates summary for", " JOB: ", WORKER.ID,
               " ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    
    write.table(x = .covariates_sum, file = gsub("survtime", "covariate_summary", MODEL.OUTPUT.FILE.PATH), 
                sep = ",", col.names = TRUE, row.names = FALSE)
    
    rm(list = ".covariates_sum")
    gc()
    
    cat(paste0("Covariate summary writing for", " JOB: ", WORKER.ID,
               " ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    
    cat(paste0("Generating correlations for", " JOB: ", WORKER.ID,
               " ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    
    .correlations <- data.frame(cor(.covariates))
    .correlations$var <- row.names(.correlations)
    
    cat(paste0("Correlations generated for", " JOB: ", WORKER.ID,
               " ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    
    cat(paste0("Writing correlations for", " JOB: ", WORKER.ID,
               " ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    write.table(x = .correlations, file = gsub("survtime", "covariate_correlations", MODEL.OUTPUT.FILE.PATH), 
                sep = ",", col.names = TRUE, row.names = FALSE)
    
    rm(list = ".correlations")
    gc()
    
    cat(paste0("Correlations written for", " JOB: ", WORKER.ID,
               " ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
  }
  
  
  #determine how many subjects are in the data
  .n_subjects <- length(unique(.covariates$id))
  
  #find number of observations per person
  .n_observations_per <- unique(table(.covariates$id))
  
  ### process covariates ####
  
  #viral load
  #split into quintiles
  .covariates$log_vln_quint <- cut(.covariates$log_vln, 
                                   breaks = quantile(.covariates$log_vln, seq(from = 0, to = 1, by = 0.2)),
                                   labels = c("first", "second", "third", "fourth", "fifth"))
  .covariates$log_vln_2 <- .covariates$log_vln_3 <- .covariates$log_vln_4 <- .covariates$log_vln_5 <- 0
  .covariates$log_vln_2[.covariates$log_vln_quint == "second"] <- 1
  .covariates$log_vln_3[.covariates$log_vln_quint == "third"] <- 1
  .covariates$log_vln_4[.covariates$log_vln_quint == "fourth"] <- 1
  .covariates$log_vln_5[.covariates$log_vln_quint == "fifth"] <- 1
  
  #cd4
  #use cutpoints from paper
  .covariates$cd4_cuts <- cut(.covariates$cd4, 
                              breaks = c(-1e6, 50, 100, 200, 350, 1e6),
                              labels = c("lt_50", "50_100", "100_200", "200_350", "350+"))
  .covariates$ind_cd4_50_100 <- 
    .covariates$ind_cd4_100_200 <-
    .covariates$ind_cd4_200_350 <- 
    .covariates$ind_cd4_350_500 <- 0
  
  .covariates$ind_cd4_50_100[.covariates$cd4_cuts == "50_100"] <- 1
  .covariates$ind_cd4_100_200[.covariates$cd4_cuts == "100_200"] <- 1
  .covariates$ind_cd4_200_350[.covariates$cd4_cuts == "200_350"] <- 1
  .covariates$ind_cd4_350_500[.covariates$cd4_cuts == "350+"] <- 1
  
  #bmi
  #use cutpoints from paper
  
  #use cutpoints from paper
  .covariates$bmi_cuts <- cut(.covariates$bmi, 
                              breaks = c(-1e6, 20, 25, 30, 1e6),
                              labels = c("lt_20", "20_25", "25_30", "gt30"))
  .covariates$ind_bmi_lt_20 <- 
    .covariates$ind_bmi_25_30 <-
    .covariates$ind_bmi_gt_30 <- 0
  
  .covariates$ind_bmi_lt_20[.covariates$bmi_cuts == "lt_20"] <- 1
  .covariates$ind_bmi_25_30[.covariates$bmi_cuts == "25_30"] <- 1
  .covariates$ind_bmi_gt_30[.covariates$bmi_cuts == "gt30"] <- 1
  
  #generate dummy variables for race
  .covariates$race_black <- 0
  .covariates$race_other <- 0
  # table(.covariates$race_black)
  .covariates$race_black[.covariates$racecatexp == 1] <- 1
  .covariates$race_other[.covariates$racecatexp == 2] <- 1
  
  #make age constant for each subject
  .first_ages <- .covariates$age[ ( (0:(.n_subjects - 1) )*.n_observations_per ) + 1]
  .covariates$age <- rep(.first_ages, each = .n_observations_per)
  
  #give everyone a frailty term
  .covariates$frailty <- rep(rnorm(.n_subjects, mean = 0, sd = 0.5 ), each = .n_observations_per)
  
  #give everyone a t0 for merging in with final covariate set
  .covariates$t0 <- rep(0:(.n_observations_per - 1), times = .n_subjects)
  
  # 
  # .covariates <- ddply(.covariates, .(id), function(.df) {
  #   
  #   #assign everyone a constant age
  #   if(length(unique(.covariates$age) ) == 1) {
  #     .df$age <- round(runif(1,min = 40, max = 70),0)
  #   } else {
  #     .df$age <- .df$age[1]
  #   }
  #   
  #   #give everyone a frailty term
  #   .df$frailty <- rnorm(n = 1, mean = 0, sd = 0.5 )
  #   
  #   #give everyone a t0 for merging in with final covariate set
  #   .df$t0 <- 0:(nrow(.df) - 1)
  #   
  #   return(.df)
  #   
  # })
  # .covariates$t0[1:201]
  
  .beta$frailty <- 1
  
  cat(paste0("Covariates processed for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
  
  #generate job-specific betas
  .beta_pool <- ldply(permn(c(log(1.05), log(1.15), log(1.2), log(1.5), log(2))))
  names(.beta_pool) <- c("d_abac", "d_ataz", "d_dida", "d_efav", "d_emtr")
  
  .beta[ ,names(.beta_pool)] <- .beta_pool[.beta_set, names(.beta_pool)]     
  
  use_PermAlgo <- T
  if (use_PermAlgo) {
    
    cat(paste0("Generating survival and censoring times for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    
    #load empirical survival times
    # source("/share/PI/manishad/PCORI/genSurv/r/empiricalSurvivalTimes.R")
    # source("~/shexport/PCORI/genSurv/r/empiricalSurvivalTimes.R")
    
    #rescale and bin into integers empirical survival times
    # .stime_scale_factor <- ( max(obs_period_days)/.n_observations_per )
    # .binned <- ceiling( x = obs_period_days / .stime_scale_factor )
    .survival_median <- 60
    .survival_shape <- 1.5
    .survival_scale <- .survival_median / (log(2) ^ (1/.survival_shape))
    # median(rweibull(10000, shape = .shape, scale = .scale))
    # .survival_times <- rweibull(10000, shape = .shape, scale = .scale)
    .survival_times <- ceiling(rweibull(.n_subjects, shape = .survival_shape, scale = .survival_scale))
    # hist(.survival_times)
    # quantile(.survival_times, probs = 0.05)
    # .p_censor <- 0.95
    # .max_censor_time <- qweibull(1 - .p_censor, shape = .survival_shape, scale = .survival_scale)
    # .mean_censor_time <- quantile(.survival_times, (1 - .p_censor ) )
    .censor_shape <- 1
    .censor_scale <- 10
    .censor_times <- ceiling(rweibull(.n_subjects, shape = .censor_shape, scale = .censor_scale))
    # hist(.censor_times)
    # mean(.survival_times > .censor_times)
    # hist(.survival_times[.survival_times < .censor_times])
    #sample a survival time for each person using distribution of empirical survival times
    # .survival_times <- sample(1:.n_observations_per, .n_subjects, replace = TRUE, 
    # prob = table(.binned)/length(.binned))
    cat(paste0("Survival and censoring times generated for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    cat(paste0("Running permalgo for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    
    #use permalgorithm to generate survival times
    .data_none <- permalgorithm(numSubjects = .n_subjects, maxTime = .n_observations_per, 
                                Xmat = as.matrix(.covariates[ , names(.beta)]), XmatNames = names(.beta),
                                beta = t(.beta), censorRandom = .censor_times,
                                eventRandom = .survival_times )
    
    cat(paste0("Permalgo complete for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    
    names(.data_none)[c(1,2,4,5)] <- c("id", "d", "t0", "t")
    
    # hist(rweibull(.n_subjects, shape = 2, scale = 75))
    
    .data_to_merge <- .covariates[ ,c("id", "t0", setdiff(names(.covariates), names(.data_none)))]
    
    .data_none <- merge(x = .data_none,
                        y = .data_to_merge,
                        by = c("id", "t0") )
    
    
    #rename variables to fit with previously existing code
    # .data_none[1:10,1:5]
    1 - sum(.data_none$d) / .n_subjects
    .data_none$proportion_censored <- 1 - sum(.data_none$d) / .n_subjects
    .data_none$source_file <- .current_data_file
    # .data_none$survtime_scale_factor <- .stime_scale_factor
    
    #     events <- ddply(.data_none, .(id), function(.df) {
    #       return(.df[nrow(.df), 1:4])
    #     })
    #     
    #     mean(events$d)
    #     hist(events$Fup)
    #     summary(events$Fup)
    #     .data_none[50:60,1:10]
    #     .data_none[1,]
    #     
    .data_none <- .data_none[order(.data_none$id, .data_none$t0), ]
    
    # .survival_times <- ldply(.covars_list, .gen_y, .g_inv_t,  .g_inv_t_min, .g_inv_t_max)
    cat(paste0("Survival time data generated for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    write.table(.data_none, file = paste0(.data_write_directory, "/SURV_", str_pad(.beta_set, 4, "left", "0"), 
                                          "_", .current_data_file), sep = ",", row.names = F, col.names = T)
    cat(paste0("Survival time data written for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    
    # cat(paste0("Survival times generated w/PermAlgo for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    
  }
  use_hendry <- F
  if(use_hendry) {
    #generate time-dependent lambda
    # .covariates <- .covariates[ , c(names(.beta), "id")]
    .covariates$linpred <- as.matrix(.covariates[ , names(.beta)]) %*% t(.beta)
    .covariates$frailty <- rep(rnorm(n = .n_subjects, mean = 0, sd = 0.5 ), each = .n_observations_per)
    
    .covariates$xB <-  exp(.covariates$linpred + .covariates$frailty)
    
    #put data frame in list form by id
    .covars_list <- dlply(.covariates, .(id))
    
    #define simulation parameters
    .ZBbar <- mean(.covariates$linpred)
    .nu <- 2.0
    .mediangoal <- 75
    #was 400 2015/06/22
    #nu = 10 and mediangoal = 50 work
    .lambda <- (log(2)/exp(.ZBbar))*.mediangoal^(-.nu)
    
    # the g function is defined as the inverse of the baseline cummulative hazard from
    ## a Weibull with shape nu and scale lambda defined above
    .g <- function(x){
      ((1/.lambda)*x)^(1/.nu)
    }
    .g_inv <- function(x){
      .lambda*(x^.nu)
    }  
    
    .t <- 0:(.n_observations_per-1)
    .t_diff <- (.t[-1] - .t[1:(length(.t) - 1)])[-(length(.t) - 1)]
    .g_inv_t <- .g_inv(.t)
    .g_inv_t_diff <- (.g_inv(.t[-1]) - .g_inv(.t[1:(length(.t) - 1)]))[-(length(.t) - 1)]
    
    #CREATING THE BOUNDS OF TRUNCATION
    .t_max <- .n_observations_per
    .t_min <- 1
    
    .g_inv_t_max <- .g_inv(.t_max)
    .g_inv_t_min <- .g_inv(.t_min)
    
    
    #K function applies ACCEPT-REJECT algorithm
    .k <- function(..x, ..m, ..M, ..rates, ..t){
      ifelse(..x <= ..m | ..x >= ..M, 0, dpexp(..x, ..rates, ..t))
    }
    
    #define survival time generation function
    .gen_y <- function(.x, .g_inv_t,  .g_inv_t_min, .g_inv_t_max) {
      .x1 <- .x$xB
      .d <- ppexp(.g_inv_t_max, .x1, .g_inv_t) - ppexp(.g_inv_t_min, .x1, .g_inv_t)
      .M <- 1 / .d
      .r <- 60
      .count<-0
      #counter of times repeat is run
      while (.count<1000) {
        .count <- .count+1
        .y <- rpexp(.r, .x1, .g_inv_t)
        .u <- runif(.r)
        .t <- .M * (.k(.y, .g_inv_t_min, .g_inv_t_max, .x1, .g_inv_t) / .d / dpexp(.y, .x1, .g_inv_t))
        .y <- .y[.u <= .t][1]
        if (!is.na(.y)) {break}
      }
      # print(.count)
      .y
    }
    
    cat(paste0("Covariates processed for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    .survival_times <- ldply(.covars_list, .gen_y, .g_inv_t,  .g_inv_t_min, .g_inv_t_max)
    cat(paste0("Survival times generated for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    .survival_times$g_y <- .g(.survival_times[ ,2])
    # summary(.survival_times$g_y/0.03921569)
    hist(.survival_times$g_y)
    
    if (sum(is.na(.survival_times$V1)) == 0) {
      
      #create uncensored model-ready dataset
      .data_none <-  ldply(1:.n_subjects, function(..subject, ..survival_times, ..covars_list) {
        ..survival_time <- ceiling(..survival_times$g_y[..subject])
        ..to_return <- ..covars_list[[..subject]][1:..survival_time, ]
        ..to_return$id <- ..subject
        ..to_return$t <- c(1:..survival_time)
        ..to_return$t0 <- 0:(..survival_time - 1)
        ..to_return$d <- c( rep(0, ..survival_time - 1), 1)
        ..to_return$proportion_censored = 0
        
        return(..to_return)
        
      }, .survival_times, .covars_list)
      
      .data_none$source_file <- .current_data_file
      
      cat(paste0("Uncensored data set generated for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
      write.table(.data_none, file = paste0(.data_write_directory, "/SURV_", .current_data_file), sep = ",", row.names = F, col.names = T)
      cat(paste0("Uncensored data written for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    }
  }
  
  #initialize model container list
  .uncensored_results <- list()
  
  #fit some models
  .pooled_cox_coefs <- names(.beta)[ !grepl( "race|frailty", names(.beta))]
  .pooled_cox_formula <- as.formula(paste0("Surv(t0, t, d) ~ ", paste0(.pooled_cox_coefs, collapse = " + "), " + cluster(id) "))
  if (ncol(.beta) == 1) {
    .pooled_cox_results <- data.frame(matrix(summary(coxph(.pooled_cox_formula, data = .data_none))$coef[ ,c(1,4,6)], ncol = 3))
    names(.pooled_cox_results) <- c("coef", "se", "p")
    .pooled_cox_results$var <- names(.beta)[1]
  } else {
    .pooled_cox_results <- data.frame(summary(coxph(.pooled_cox_formula, data = .data_none))$coef[ ,c(1,4,6)])
    names(.pooled_cox_results) <- c("coef", "se", "p")
    .pooled_cox_results$var <- row.names(.pooled_cox_results)
  }
  
  .pooled_cox_results$type <- "Pooled Cox"
  .pooled_cox_results$proportion_censored <- 1 - sum(.data_none$d) / .n_subjects
  
  .uncensored_results[[1]] <- .pooled_cox_results
  
  cat(paste0("Pooled Cox fit for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
  ################################      
  # .na_coef_vars <- .pooled_cox_results$var[is.na(.pooled_cox_results$coef)]
  .coxme_coefs <- names(.beta)[ !grepl( "race|frailty", names(.beta))]
  # .coxme_coefs <- .pooled_cox_coefs[!(.pooled_cox_coefs %in% .na_coef_vars)]
  
  .coxme_formula <- as.formula(paste0("Surv(t0, t, d) ~ ", paste0(.coxme_coefs, collapse = " + "), " + (1|id) "))
  .coxme_results <- data.frame(coxmeCoefs(coxme(.coxme_formula, data = .data_none)))
  names(.coxme_results) <- c("coef", "se", "p")
  .coxme_results$var <- row.names(.coxme_results)
  .coxme_results$type <- "Frailty Cox"
  .coxme_results$proportion_censored <- 1 - sum(.data_none$d) / .n_subjects
  
  #add in NA vars
  # .na_rows <- .coxme_results[1:length(.na_coef_vars), ]
  # .na_rows$var <- .na_coef_vars
  # .na_rows[ ,c("coef", "se", "p")] <- NA
  # .coxme_results <- rbind(.coxme_results, .na_rows)
  .uncensored_results[[2]] <- .coxme_results
  
  cat(paste0("Frailty Cox fit for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
  ################################  
  if(FALSE) {
    .pooled_poisson_formula <- as.formula(paste0("d ~ ", paste0(names(.beta), collapse = " + ")))
    .pooled_poisson_results <- data.frame( summary(geeglm(data = .data_none, formula = .pooled_poisson_formula, id = id, 
                                                          corstr = "independence", family = "poisson"))$coef[-1,c(1,2,4)])
    names(.pooled_poisson_results) <- c("coef", "se", "p")
    .pooled_poisson_results$var <- row.names(.pooled_poisson_results)
    .pooled_poisson_results$type <- "Pooled Poisson"
    .pooled_poisson_results$proportion_censored <- 1 - sum(.data_none$d) / .n_subjects
    
    .uncensored_results[[3]] <- .pooled_poisson_results
    
    cat(paste0("Pooled Poisson fit for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
    
    .pooled_logistic_formula <- as.formula(paste0("d ~ ", paste0(names(.beta), collapse = " + ")))
    .pooled_logistic_results <- data.frame( summary(geeglm(data = .data_none, formula = .pooled_logistic_formula, id = id, 
                                                           corstr = "independence", family = "binomial"))$coef[-1,c(1,2,4)])
    names(.pooled_logistic_results) <- c("coef", "se", "p")
    .pooled_logistic_results$var <- row.names(.pooled_logistic_results)
    .pooled_logistic_results$type <- "Pooled Logistic"
    .pooled_logistic_results$proportion_censored <- 1 - sum(.data_none$d) / .n_subjects
    
    .uncensored_results[[4]] <- .pooled_logistic_results
    
    cat(paste0("Pooled Logistic fit for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
  }
  
  .uncensored_results <- ldply(.uncensored_results)
  
  cat(paste0("Uncensored models fit for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
  
  
  .to_write <-  .uncensored_results
  .to_write$rep <- .rep
  
  
  .beta_merge <- data.frame(beta = t(.beta))
  .beta_merge$var <- rownames(.beta_merge)                   
  
  .to_write <-  merge(x = .to_write,
                      y = .beta_merge,
                      by = "var",
                      all.x = T)
  
  write.table(.to_write, file = MODEL.OUTPUT.FILE.PATH, 
              append = T, row.names = F, col.names = F, sep = ",")
  # paste0(.results_write_directory, "/", WORKER.ID, "_", .sim_results_name), append = F)
  cat(paste0("Results written for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
  cat(paste0("-------------------END-REP------------------------\n"), file = LOG.FILE.PATH, append = T)
  
  # } else {
  # cat(paste0("Convergence fail for rep ", .rep, " at ", Sys.time(), "\n"), file = LOG.FILE.PATH, append = T)
  # }
  gc()
}, covariate_data_path, betas, results_write_directory, sim_results_name, log_write_directory, data_write_directory, beta_set)
