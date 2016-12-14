
######################### SET SIMULATION PARAMETERS #########################
# if need to install paxkages
#install.packages("plyr")

library(devtools)
library(MBESS)
library(ICC)
library(miscTools)
library(corpcor)
library(psych)
library(ggplot2)


# source code from GitHub
source_url("https://raw.githubusercontent.com/qsuProjects/PCORI/master/r/GENCOV/load_functions.R")

source("https://raw.githubusercontent.com/qsuProjects/PCORI/master/r/GENCOV/jointly_generate_binary_normal_modified_v2.R", 
       local=TRUE)

#############################################
##        Set Simulation Parameters        ##
#############################################

# name prefix for all datasets
name_prefix = "ex2"

# n (number of subjects)
n.Subj = 10000

# obs (number of observations per subject)
obs = 150

# n.Reps (number of datasets to generate)
n.Reps = 100

# n.Drugs (number of drug variables)
n.Drugs = 15


#############################################
## Read in Correlation and Parameter Files ##
#############################################

# within-subject correlation matrix
wcor = read.csv("/Users/mmrath/Documents/Projects/pcori/paper2/power analysis/PCORI_Simulated_Data_Set no slopes/ex2_wcor_time_vary.csv", header=FALSE)[-1,-1]

# population correlation matrix
pcor = read.csv("/Users/mmrath/Documents/Projects/pcori/paper2/power analysis/PCORI_Simulated_Data_Set no slopes/ex2_pcor_time_vary.csv", header=TRUE)[,-1]

# read in and complete parameters dataframe
parameters = complete_parameters(read.csv("/Users/mmrath/Documents/Projects/pcori/paper2/power analysis/PCORI_Simulated_Data_Set no slopes/ex2_parameters_time_vary.csv"), n.Subj)

#cat.parameters = read.csv("C:/Users/npurington/Documents/PCORI/Maria Paper/ex1_categorical_parameters.csv")


#############################################
##               Simulate!                 ##
#############################################

### following two parameters defined here because we can't change source code
### in program make_one_dataset change line 
# bins.prop = proportionize(bins, zero, one)  # proportionized version of the binaries
### to
# bins.prop = proportionize(bins)  # proportionized version of the binaries
### unless we want to change the defaults for variables zero and one (which are the values below)
zero=0.00001
one=0.999

ptm0 <- proc.time()

# simulate results
sim = repeat_sim(n=n.Subj, obs=obs, parameters=parameters, prop.target=NULL,
                       mean.target=NULL, n.Drugs=n.Drugs, 
                       pcor=pcor, wcor=wcor, n.Reps=n.Reps,
                       write.data=TRUE,
                       #name_prefix= paste( .name_prefix, WORKER.ID, sep="_" ),  # used with Sherlock
                       name_prefix=name_prefix,
                       cat.parameters=NULL)

# extract dataset
d = sim$data
head(d)

proc.time() - ptm0
#~2.7 days

######################### QUICK TOUR THROUGH THE SIMULATED DATA ########################

##### Example of Normal Variable Clustered within a Subject #####
temp = d[ d$id %in% 1:10, ]  # look at first 6 subjects

### create time variable; everyone got 150 records
temp$t <- rep(c(1:150),10)

ggplot( data=temp, aes(x=t, y=bmi) ) + geom_line() + facet_grid(id~.) 

ggplot( data=temp, aes(x=t, y=log_vln) ) + geom_line() + facet_grid(id~.) 

ggplot( data=temp, aes(x=t, y=bps) ) + geom_line() + facet_grid(id~.) 

table(temp$d_abac) 

##### Correlation Matrix Across Subjects #####
# extract observed across-subjects correlation matrix
sim$corr

# how biased were the correlations?
sim$corr.bias


