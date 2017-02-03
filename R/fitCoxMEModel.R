
library(survival)
library(coxme)

#returns coefficients from coxme object
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

#initialize model results file
#we just write the column names and then for each rep
#we append model results
cat("var,coef,se,p,type,rep,beta\n", file = "path/to/where/you're/storing/your/results", append = F)


################################################
# the following block gets run once on each
# simulated data set after survival times have
# been generated
# note that if you're running on Sherlock
# you have to use R/3.0.2 because coxme is
# very unstable on the default version of R on Sherlock
################################################

{
#betas should contain the same beta that was used to generate survival times
coxme_coefs <- names(betas)

#create formula object using betas for use in coxme()
coxme_formula <- as.formula(paste0("Surv(t0, t, d) ~ ", paste0(coxme_coefs, collapse = " + "), " + (1|id) "))

#run a coxme model using survival data contained in survival_data and model specified
#by coxme_formula, extract its coefficients, convert them to a data frame and store in 
#coxme_results
coxme_results <- data.frame(coxmeCoefs(coxme(coxme_formula, data = survival_data)))

#rename columns
names(coxme_results) <- c("coef", "se", "p")

#add column with variable names
coxme_results$var <- row.names(coxme_results)
coxme_results$type <- "Frailty Cox"

#get results ready to write
to_write <-  coxme_results
to_write$rep <- rep #rep corresponds to the replication number of this data set in the power simulation

#merge in beta values
beta_merge <- data.frame(beta = t(beta))
beta_merge$var <- rownames(beta_merge)                   

to_write <-  merge(x = to_write,
                    y = beta_merge,
                    by = "var",
                    all.x = T)

write.table(to_write, file = "path/to/where/you're/storing/your/results", 
            append = T, row.names = F, col.names = F, sep = ",")

}