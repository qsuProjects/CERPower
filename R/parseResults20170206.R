
openFilesInDirectory <- function(directory,
                                 match_string,
                                 merge = FALSE,
                                 sep = ",",
                                 na.strings = "NA",
                                 header = T,
                                 fill=T,
                                 skip = 0) {
  require(plyr)
  file_array <-  paste0(directory, "/", list.files(directory)[grep(pattern=match_string, list.files(directory))])
  
  data_list <- llply(file_array, function(file_path, delim_str) {
    cat(file_path, "\n")
    to_return <- read.table(file = file_path, header = header, sep = sep, stringsAsFactors = FALSE, fill=fill, quote="\"", na.strings = na.strings, skip = skip )
    
    if(nrow(to_return) > 0 ) {
      to_return["loaded_file_name"] <- tail(strsplit(file_path, "/")[[1]],1)
      return(to_return)
    } else {
      cat(tail(strsplit(file_path, "/")[[1]],1), " was empty\n")
      warning(tail(strsplit(file_path, "/")[[1]],1), " was empty\n")
      to_return["loaded_file_name"] <- NULL
      return(to_return)
      
    } 
  }, delim_str)
  
  if(merge |  length(file_array) == 1) {
    data_list <- ldply(data_list, identity)
  }
  
  return(data_list)
}


library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

#the number of drugs you're testing
alpha_mod <- 15

#load model results files
all_results[1,]
results_path <- "path to your model fit result files"

#open and combine files
all_results <- openFilesInDirectory(directory = results_path, 
                                    sep = ",", match_string = "survtime", 
                                    merge = TRUE, header = TRUE)

#make sure results are numeric
all_results$coef <- as.numeric(all_results$coef)
all_results$se <- as.numeric(all_results$se)
all_results$p <- as.numeric(all_results$p)
all_results$beta <- as.numeric(all_results$beta)

#define beta sets (eg unique combinations of beta x drug
all_results$beta_set <- as.numeric(gsub(".{1,}beta_set_|_core.{1,}", "", all_results$loaded_file_name))

#define confidence intervals and an indicator for coverage
all_results$lcl <- all_results$coef - 1.96 * all_results$se
all_results$ucl <- all_results$coef + 1.96 * all_results$se
all_results$covered <- as.numeric((all_results$lcl <= all_results$beta) & 
                                    (all_results$ucl >= all_results$beta))


#look only at drugs (which are variables whose names start with "d_"
drugs <- filter(all_results, grepl("^d_", var))

#also, grab drugs with nonzero betas (not sure if this is still used)
sig_drugs <- filter(drugs, beta != 0)

#summaryize results by rep
rep_summary <- ddply(drugs, .(rep, beta_set, type), function(.rep_data, .alpha_mod) {
  
  .sigs <- filter(.rep_data, beta != 0)
  .nonsigs <- filter(.rep_data, beta == 0)
  
  data.frame(all_sig_detected = as.numeric(mean(.sigs$p <= 0.05/.alpha_mod) == 1),
             any_four_sig = as.numeric(sum( .sigs$p <= 0.05/.alpha_mod) == 4),
             any_three_sig = as.numeric(sum( .sigs$p <= 0.05/.alpha_mod) == 3),
             any_two_sig = as.numeric(sum( .sigs$p <= 0.05/.alpha_mod) == 2),
             any_one_sig = as.numeric(sum( .sigs$p <= 0.05/.alpha_mod) >= 1),
             at_least_two_sig = as.numeric(sum( .sigs$p <= 0.05/.alpha_mod) >= 2),
             prop_sig_detected = mean(.sigs$p <= 0.05/.alpha_mod),
             n_sig_detected = sum(.sigs$p <= 0.05/.alpha_mod),
             any_nonsig_detected = as.numeric( sum(.nonsigs$p <= 0.05/.alpha_mod) > 0),
             n_nonsig_detected = as.numeric(sum(.nonsigs$p <= 0.05/.alpha_mod)),
             prop_censored = .rep_data$proportion_censored[1])
  
}, alpha_mod)

summary(rep_summary$n_nonsig_detected)

write.table(x = rep_summary, file = paste0("output path for summary data/repSummary.csv"), 
            sep = ",", col.names = TRUE, row.names = FALSE)



prev_summary <-  read.table(file = paste0("output/summaryData/prevSummary.csv"), 
                            sep = ",", header = TRUE, stringsAsFactors = FALSE)
all_results[1,]


cov_sums_agg <- read.table(file = paste0("output/summaryData/covariateSummary.csv"), 
                           sep = ",", header = TRUE, stringsAsFactors = FALSE)


all_results <- merge(x = filter(all_results, type == "Frailty Cox"),
                     y = cov_sums_agg[ ,c(1,2)],
                     by.x = "var",
                     by.y = "var",
                     all.x = TRUE)

all_results$prevalence <- printDec(all_results$mean_mean)
all_results$d_prevalence <- paste0("Drug ", match(all_results$prevalence, 
                                                   c("0.05", "0.10", "0.20", "0.30", "0.40")),
                                    "\nPrevalence: ", all_results$prevalence)


drugs <- filter(all_results, grepl("^d_", var))

drugs[1,]

table(drugs$rep)

write.table(x = drugs, file = paste0("output/summaryData/drugResults.csv"), 
            sep = ",", col.names = TRUE, row.names = FALSE)

