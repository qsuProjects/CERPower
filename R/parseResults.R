if(FALSE) {
  rm(list = ls())
  source("~/.Rprofile")
}

system("scp kikapp@sherlock:/scratch/PI/manishad/PCORI/power/modelFit/* ~/shexport/PCORI/power/modelFit/", 
       intern = TRUE, ignore.stdout = FALSE, ignore.stderr = FALSE)

setwd("~/shexport/PCORI/power/")
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

alpha_mod <- 15
#load model results files
all_results[1,]
all_results <- openFilesInDirectory(directory = paste0("modelFit/"), sep = ",", match_string = "survtime", 
                     merge = TRUE, header = TRUE)
all_results$coef <- as.numeric(all_results$coef)
all_results$se <- as.numeric(all_results$se)
all_results$p <- as.numeric(all_results$p)
all_results$proportion_censored <- as.numeric(all_results$proportion_censored)
all_results$beta <- as.numeric(all_results$beta)

all_results$beta_set <- as.numeric(gsub(".{1,}beta_set_|_core.{1,}", "", all_results$loaded_file_name))

colClasses((all_results))
all_results$lcl <- all_results$coef - 1.96 * all_results$se
all_results$ucl <- all_results$coef + 1.96 * all_results$se
all_results$covered <- as.numeric((all_results$lcl <= all_results$beta) & 
                                    (all_results$ucl >= all_results$beta))


drugs <- filter(all_results, grepl("^d_", var))
sig_drugs <- filter(drugs, beta != 0)

# .rep_data <- filter(drugs, rep == 1 & beta_set == 1 & type == "Frailty Cox")
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

write.table(x = rep_summary, file = paste0("output/summaryData/repSummary.csv"), 
            sep = ",", col.names = TRUE, row.names = FALSE)

################################################
# summarize by prevalence
################################################
drugs[1,]

# .rep_data <- filter(sig_drugs, var == "d_abac" & type == "Frailty Cox")
prev_summary <- ddply(sig_drugs, .(var, type), function(.rep_data, .alpha_mod) {
  
  .sigs <- filter(.rep_data, beta != 0)
  # .nonsigs <- filter(.rep_data, beta == 0)
  # .sigs[1:10,]
  data.frame(#all_sig_detected = as.numeric(mean(.sigs$p <= 0.05/.alpha_mod) == 1),
             # any_four_sig = as.numeric(sum( .sigs$p <= 0.05/.alpha_mod) == 4),
             # any_three_sig = as.numeric(sum( .sigs$p <= 0.05/.alpha_mod) == 3),
             # any_two_sig = as.numeric(sum( .sigs$p <= 0.05/.alpha_mod) == 2),
             # any_one_sig = as.numeric(sum( .sigs$p <= 0.05/.alpha_mod) >= 1),
             # at_least_two_sig = as.numeric(sum( .sigs$p <= 0.05/.alpha_mod) >= 2),
             prop_sig_detected = mean(.sigs$p <= 0.05/.alpha_mod),
             n_sig_detected = sum(.sigs$p <= 0.05/.alpha_mod))#,
             # any_nonsig_detected = as.numeric( sum(.nonsigs$p <= 0.05/.alpha_mod) > 0),
             # n_nonsig_detected = as.numeric(sum(.nonsigs$p <= 0.05/.alpha_mod)),
             # prop_censored = .rep_data$proportion_censored[1])
  
}, alpha_mod)

# summary(rep_summary$n_nonsig_detected)

write.table(x = prev_summary, file = paste0("output/summaryData/prevSummary.csv"), 
            sep = ",", col.names = TRUE, row.names = FALSE)

#summarize results by var
# .var_data <- filter(sig_drugs, var == "d_abac")
.var_data[1,]
var_summary <- ddply(sig_drugs, .(var, beta_set, type), function(.var_data) {
  data.frame(prop_detected = mean(.var_data$p <= 0.05),
             n_detected = sum(.var_data$p <= 0.05),
             coverage_p = mean(.var_data$covered),
             mean_coef = mean(.var_data$coef),
             mean_se = mean(.var_data$se),
             beta = .var_data$beta[1],
             mean_abs_bias = mean(.var_data$coef - .var_data$beta),
             n = nrow(.var_data)) 
})

var_summary[1,]

write.table(x = var_summary, file = paste0("output/summaryData/varSummary.csv"), 
            sep = ",", col.names = TRUE, row.names = FALSE)


betaprev_summary <- ddply(sig_drugs, .(var, beta, type), function(.var_data) {
  data.frame(prop_detected = mean(.var_data$p <= 0.05),
             n_detected = sum(.var_data$p <= 0.05),
             coverage_p = mean(.var_data$covered),
             mean_coef = mean(.var_data$coef),
             mean_se = mean(.var_data$se),
             beta = .var_data$beta[1],
             mean_abs_bias = mean(.var_data$coef - .var_data$beta)) 
})

write.table(x = betaprev_summary, file = paste0("output/summaryData/betaPrevSummary.csv"), 
            sep = ",", col.names = TRUE, row.names = FALSE)
#summarize covariates values

cov_sums <- openFilesInDirectory(directory = paste0("modelFit/"), sep = ",", match_string = "summary", 
                                    merge = TRUE, header = TRUE)

summary(cov_sums$mean[cov_sums$var == "ever_abac"])
cov_sums[1,]
cov_sums_agg <- ddply(cov_sums, .(var), function(.var_data) {
  data.frame(mean_mean = printDec(mean(.var_data$mean), 4),
             median_mean = printDec(mean(.var_data$mean), 4),
             sd_mean = printDec(sd(.var_data$mean), 4)) 
})

write.table(x = cov_sums_agg, file = paste0("output/summaryData/covariateSummary.csv"), 
            sep = ",", col.names = TRUE, row.names = FALSE)

