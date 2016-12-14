
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

