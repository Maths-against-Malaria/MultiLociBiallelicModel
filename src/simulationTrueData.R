# Title        : Simulation study using haplotype frequencies from real-world data
# Objective    : Implement the EM-algorithm on simulated data and save the estimates
# Created by   : Christian Tsoungui Obama
# Created on   : 03.04.21
# Last modified: 22.12.21

# Loading library
library(xlsx)
library(dplyr)

# Relative path
#path <- "/Volumes/GoogleDrive-117934057836063832284/My Drive/Maths against Malaria/Christian/Models/MultiLociBiallelicModel/"
path <- '/Volumes/GoogleDrive/My Drive/Maths against Malaria/Christian/Models/MultiLociBiallelicModel/'

# Loading external ressources
source(paste0(path, "src/nbiallelicModel.R"))          ## Loading Model
source(paste0(path, "src/dataGenerator.R"))            ## Loading the data generaor for model

df <- read.xlsx(paste0(path,'dataset/MutFreq2005-2017_CT.xlsx'), 1)

# Data origin
name <- 'Kenya'

## True frequencies adjustments (to compensate for the removal of the )
df[17, "est2010"] <- "0.081"

## Data preprocessing
df1 <- df %>%
        select(c(1:11, 14, 17)) %>%
        mutate_at(1:10, as.factor)

for (i in 11:13){
    sel <- df1[,i] == "-"
    df1[,i][sel] = "0"
}

df1 <- df1 %>%
            mutate_at(11:13, as.numeric)

## Defining the matrix of true haplotype frequencies
est_years <- c(2005, 2010, 2017)

# True Poisson parameter
lbdavec <- c(0.1, 0.25, 0.5, 1, 1.5, 2, 2.5)
NLbd <- length(lbdavec)

# Number of loci considered
NumbLoci <- c(10)
Nn <- length(NumbLoci)
nvec <- matrix(NumbLoci, nrow = 1, ncol = Nn)

# Number of possible haplotypes
Hvec <- 2^nvec
NH <- length(Hvec)
Hvecpo <- Hvec + 1

# Sample sizes considered
Nvec <- c(50, 100, 150, 200, 500)
NN <- length(Nvec)

# Number of estimates generated in the simulation
NEst <- 10000

# Number of distributions of true frequencies for each number of loci (we remove 2017 estimates)
NFreq <- length(est_years)-1

# Extra parameters
ParExtra <- list(NLbd, Nn, Hvec, NN, NEst, NFreq)

# True haplotype frequencies (without estimates for 2017)
Pvec <- vector(mode="list", length=Nn)
for (i in 1:Nn){
  Pvec[[i]] <- array(0, c(1, Hvec[i]))
  tmp <- rep(0, Hvec[i])
  for (j in 1:NFreq){
    tmp <-  c(round(df1[,paste0("est", est_years[j])], 10), rep(0, (Hvec[i]-nrow(df1))))
    Pvec[[i]] <- rbind(Pvec[[i]], tmp)
  }
  Pvec[[i]] <- Pvec[[i]][-1,]
}

# True parameter
True_param <- list(Pvec, lbdavec, Nvec)

# Simulation
out <- vector(mode = "list", length = Nn)

for (i in 1:Nn){
  sizelist <- vector(mode = "list", length = NN)
  for (j in 1:NN){                                                                                ## For each value of the sample size
    lbdalist <- vector(mode = "list", length = NLbd)
    for (k in 1:NLbd){                                                                            ## For each value of the lambda parameter
      Estim <- array(0, dim = c(Hvecpo[i], NEst, NFreq))
      for (cnt in 1:NFreq){
        for (l in 1:NEst){
          infct <-  sampNew(unlist(Pvec[[i]][cnt,]) ,unlist(lbdavec[k]) ,Nvec[j], nvec[,i])       ## Generating data for the simulation
          Estim[,l,cnt] <- unlist(nbialModel(infct[[2]], infct[[1]]))                             ## Evaluating and saving the Estimates
        }
      }
      lbdalist[[k]] <- Estim
    }
    sizelist[[j]] <- lbdalist
  }
  out[[i]] <- sizelist
}

# Saving the list for post-processing
saveRDS(out, file = paste0(path, "dataset/modelEstimates", name, ".rds"))

# Saving the true parameters
saveRDS(True_param, file = paste0(path, "dataset/trueParameters", name, ".rds"))

# Saving the extra parameters
saveRDS(ParExtra, file = paste0(path, "dataset/extraParameters", name, ".rds"))

print("Simulation finished, check your data.")
