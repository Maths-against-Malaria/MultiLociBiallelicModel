# Title        : Simulation study
# Objective    : Compute the bias and coefficient of variation of the estimates,
#                and the estimates of prevalence
# Created by   : Christian Tsoungui Obama
# Created on   : 03.04.21
# Last modified: 31.01.22

# Importing libraries
library(wordspace)

# Functions
psi <- function (inp){
  inp / (1 - exp(-inp))
}

hapl <- function(n){
  H <- array(0,c(2^n,n))
  H[1:2,1] <- c(0,1)
  for(k in 2:n){
    H[(2^(k-1)+1):2^k,1:(k-1)] <- H[1:2^(k-1),1:(k-1)]
    H[(2^(k-1)+1):2^k,k] <- 1
  }
  H <- H[,n:1]
  H
}

gen_func <- function(x, lambd){
  (exp(x*lambd)-1)/(exp(lambd) - 1)
}

bias <- function(df_Estim, df_Param, name){
  # This function implements bias of mean MOI and haplotype frequencies as defined in the manuscript
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence 
  # alongside multiplicity of infection from SNPs data"

  moibiasloc <- vector(mode = "list", length = Nn)
  freqbiasloc <- vector(mode = "list", length = Nn)

  for (l in 1:numbloci){   # For each number of locus
    freqbiasSamp <- vector(mode = "list", length = NN)
    moibiasSamp  <- vector(mode = "list", length = NN)

    # Number of haplotypes
    numH   <- Hvec[l]
    numHpo <- numH + 1

    for (k in 1:NN){  # For each sample size
      freqbias_lamb <- vector(mode = "list", length = NFreq)
      moibias_lamb  <- vector(mode = "list", length = NFreq)

      for (j in 1:NLbd){ # For each true Lambda
        freq_bias <- vector(mode = "list", length = NFreq)
        moi_bias  <- vector(mode = "list", length = NFreq)

        for (i in 1:NFreq){ # For each choice of true frequency distribution
          # Access each of the 10000 estimates
          tmp1 <- df_Estim[[l]][[k]][[j]][, , i]

          # True frequencies
          tmp2 <- t(df_Param[[1]][[l]])
          tmp3 <- tmp2[, i]

          # Relative bias of haplotype freq. in percent
          freqbias <- rowMeans((tmp1[2:numHpo,]/tmp3 - 1), na.rm = TRUE) * 100

          # Remove the estimates of Lambda that yield an "Infinite mean MOI"
          mean_moi_estim <- psi(tmp1[1, ])
          
          # Relative bias of MOI in percent
          bias_moi <- (mean_moi_estim/mean_moi_true[j] - 1)
          bias_moi <- mean(bias_moi[!is.infinite(bias_moi)], na.rm = TRUE) * 100

          # Save frequency and MOI bias in lists
          freq_bias[[i]] <- freqbias
          moi_bias[[i]]  <- bias_moi
        }
        freqbias_lamb[[j]] <- freq_bias
        moibias_lamb[[j]]  <- moi_bias
      }
      freqbiasSamp[[k]] <- freqbias_lamb
      moibiasSamp[[k]]  <- moibias_lamb
    }
    freqbiasloc[[l]] <- freqbiasSamp
    moibiasloc[[l]]  <- moibiasSamp
  }

  # Saving the estimates
  saveRDS(freqbiasloc, file = paste0(path, "dataset/freqbias",name ,".rds"))
  saveRDS(moibiasloc, file = paste0(path, "dataset/moibias",name ,".rds"))
}

coefvar <- function(df_Estim, df_Param, name){
  # This function implements coefficient of variation of mean MOI as defined in the manuscript
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence 
  # alongside multiplicity of infection from SNPs data"

  moicvloc <- vector(mode = "list", length = Nn)

  for (l in 1:numbloci){   # For each number of locus
    moicvSamp  <- vector(mode = "list", length = NN)

    # Number of haplotypes
    numH   <- Hvec[l]
    numHpo <- numH + 1

    for (k in 1:NN){  # For each sample size
      moicv_lamb  <- vector(mode = "list", length = NFreq)

      for (j in 1:NLbd){ # For each true Lambda
        moi_cv  <- vector(mode = "list", length = NFreq)

        for (i in 1:NFreq){ # For each choice of true frequency distribution
          # Access each of the 10000 estimates
          tmp1 <- df_Estim[[l]][[k]][[j]][, , i]

          # Remove the estimates of Lambda that yield an "Infinite mean MOI"
          mean_moi_estim <- psi(tmp1[1, ])
          mean_moi_estim <- mean_moi_estim[!is.infinite(mean_moi_estim)]
          
          # Dimensionless coefficient of variation of MOI
          moicv <- sd(mean_moi_estim, na.rm = TRUE) / mean_moi_true[j]

          # Save moi coefficient of variation in list
          moi_cv[[i]]  <- moicv
        }
        moicv_lamb[[j]]  <- moi_cv
      }
      moicvSamp[[k]]  <- moicv_lamb
    }
    moicvloc[[l]]  <- moicvSamp
  }

  # Saving the estimates
  saveRDS(moicvloc, file = paste0(path, "dataset/moicv", name, ".rds"))
}

true_amb_prevalence <- function(true_Par, name){
  # This function implements the true ambiguous prevalence as defined in the manuscript
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence 
  # alongside multiplicity of infection from SNPs data"

  qh_loc <- vector(mode = "list", length = Nn)
  for (l in 1:numbloci){   # For each number of locus
  qh_loc[[l]] <- vector(mode = "list", length = NFreq)

    # Number of haplotypes
    numH <- Hvec[l]
    numHpo <- Hvec[l] + 1

        for (i in 1:NFreq){ # For each choice of true frequency distribution
          amb_prevalence <- array(0, dim = c(numH, NLbd))

          # Ambiguous prevalence
          ## Access each of the 10000 estimates
          tmp1 <- true_Par[[l]][[i]]

          ## For each set of estimates, compute prevalence
          for (j in 1:NLbd){
            amb_prevalence[,j] <- (exp(tmp1[1,j]) - exp(1-tmp1[2:numHpo,j])^tmp1[1,j])/(exp(tmp1[1,j])-1)
          }

          qh_loc[[l]][[i]] <- amb_prevalence
        }
  }
  # Saving the estimates
  saveRDS(qh_loc, file = paste0(path, "dataset/TrueAmbPrevalence", name, ".rds"))
  qh_loc
}

true_unamb_prevalence <- function(true_Par, ParTru, name){
 # This function implements the unambiguous prevalence as defined in the manuscript
 # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"
 
  qh_loc <- vector(mode = "list", length = Nn)
  for (l in 1:numbloci){   # For each number of locus
    qh_loc[[l]] <- vector(mode = "list", length = NFreq)
    tru_freq <- ParTru[[1]][[l]]

      # Number of haplotypes
      numH <- Hvec[l]

      # Number of loci
      nLoci <- NLoci[l]

      # Table of all possible haplotypes
      Hapl <- hapl(nLoci)

      ## For each haplotype in the table, build the set of observation Uh
      nLociUh <- nLoci+1

        for (i in 1:NFreq){ # For each choice of true frequency distribution
          rh <- array(0, dim = c(numH, NLbd))

          ## For each combination of true lanmbda and haplotype frequencies values
          tmp1 <- true_Par[[l]][[i]]
          tmp2 <- tmp1[2:(numH+1),]

          ## For each haplotype build the set Uh for all l
          trufreq_vec <- tru_freq[i,]

          for (idx in 1:numH){ 
            if(trufreq_vec[idx] != 0){
              uh <- t(array(rep(Hapl[idx,], nLoci), dim=c(nLoci, nLociUh)))
              uh[2:nLociUh,] <- (uh[2:nLociUh,]+diag(nLoci))%%2

              ## Pick the right frequencies estimates 
              pickh <- which(colSums(uh[1,] == t(Hapl))==nLoci)
              GPh <- gen_func(tmp2[pickh,], tmp1[1,])

              GPartFreq <- rep(0, NLbd)
              GFreq <- rep(0, NLbd)

              pick1 <- rep(0, nLociUh)
              pick2 <- rep(0, nLociUh)

              for(idxUh in 1:nLociUh){ 
                pick1[idxUh] <- which(colSums(uh[idxUh,] == t(Hapl))==nLoci)
              }

              pick2 <- pick1
              pick2[1] <- 0

              for(idxUh in 2:nLociUh){ 
                GPartFreq <- gen_func(tmp2[pick2[idxUh],], tmp1[1,])
                GFreq <- gen_func(colSums(tmp2[pick1[c(1,idxUh)],]), tmp1[1,])
                rh[idx,] <- rh[idx,] + GFreq - GPartFreq
              }
                rh[idx,] <- rh[idx,] - (nLoci-1)*GPh
            }else{
              rh[idx,] <- 0
            }
          }
          ## Save the prevalence in a list
          qh_loc[[l]][[i]] <- rh
        }
  }

  # Saving the estimates
  saveRDS(qh_loc, file = paste0(path,  "dataset/TrueUnambPrevalence", name, ".rds"))

  qh_loc
}

true_prevalence <- function(true_Par, ParTru, name){
 # This function implements the unambiguous prevalence as defined in the manuscript
 # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"
 
  qh_loc <- vector(mode = "list", length = Nn)
  for (l in 1:numbloci){   # For each number of locus
    qh_loc[[l]] <- vector(mode = "list", length = NFreq)
    tru_freq <- ParTru[[1]][[l]]

    #for (k in 1:NN){  # For each true sample size
     # qh_lamb <- vector(mode = "list", length = NFreq)

      # Number of haplotypes
      numH <- Hvec[l]

      # Number of loci
      nLoci <- NLoci[l]

      # Table of all possible haplotypes
      Hapl <- hapl(nLoci)

      ## For each haplotype in the table, build the set of observation Uh
      nLociUh <- nLoci+1

      #for (j in 1:NLbd){ # For each true Lambda
       # qh_freq <- vector(mode = "list", length = NFreq)

        for (i in 1:NFreq){ # For each choice of true frequency distribution
          rh <- array(0, dim = c(numH, NLbd))

          ## For each combination of true lanmbda and haplotype frequencies values
          tmp1 <- true_Par[[l]][[i]]
          tmp2 <- tmp1[2:(numH+1),]

          ## For each haplotype build the set Uh for all l
          trufreq_vec <- tru_freq[i,]

          for (idx in 1:numH){ 
            if(trufreq_vec[idx] != 0){
              uh <- t(array(rep(Hapl[idx,], nLoci), dim=c(nLoci, nLociUh)))
              uh[2:nLociUh,] <- (uh[2:nLociUh,]+diag(nLoci))%%2

              ## Pick the right frequencies estimates 
              pickh <- which(colSums(uh[1,] == t(Hapl))==nLoci)
              GPh <- gen_func(tmp2[pickh,], tmp1[1,])

              GPartFreq <- rep(0, NLbd)
              GFreq <- rep(0, NLbd)

              pick1 <- rep(0, nLociUh)
              pick2 <- rep(0, nLociUh)

              for(idxUh in 1:nLociUh){ 
                pick1[idxUh] <- which(colSums(uh[idxUh,] == t(Hapl))==nLoci)
              }

              pick2 <- pick1
              pick2[1] <- 0

              for(idxUh in 2:nLociUh){ 
                GPartFreq <- gen_func(tmp2[pick2[idxUh],], tmp1[1,])
                GFreq <- gen_func(colSums(tmp2[pick1[c(1,idxUh)],]), tmp1[1,])
                rh[idx,] <- rh[idx,] + GFreq - GPartFreq
              }
                rh[idx,] <- rh[idx,] - (nLoci-1)*GPh
            }else{
              rh[idx,] <- 0
            }
          }
          ## Save the prevalence in a list
          for (q in 1:NLbd){
            rh[,q] <- rh[,q]/sum(rh[,q])
          }
          qh_loc[[l]][[i]] <- rh
        }
  }
  # Saving the estimates
  saveRDS(qh_loc, file = paste0(path,  "dataset/TruePrevalence", name, ".rds"))

  qh_loc
}

perf_amb_prevalence <- function(df_Estim, true_prev, name){
  # This function implements the ambiguous prevalence as defined in the manuscript
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence 
  # alongside multiplicity of infection from SNPs data"

  qh_loc <- vector(mode = "list", length = Nn)
  biash_loc <- vector(mode = "list", length = Nn)
  coefvar_loc <- vector(mode = "list", length = Nn)
  for (l in 1:numbloci){   # For each number of locus
    qh_Samp <- vector(mode = "list", length = NN)
    biash_Samp <- vector(mode = "list", length = NN)
    coefvar_Samp <- vector(mode = "list", length = NN)

    # Number of haplotypes
    numH <- Hvec[l]
    numHpo <- numH + 1

    for (k in 1:NN){  # For each true sample size
      qh_lamb <- vector(mode = "list", length = NFreq)
      biash_lamb <- vector(mode = "list", length = NFreq)
      coefvar_lamb <- vector(mode = "list", length = NFreq)

      for (j in 1:NLbd){ # For each true Lambda
        qh_freq <- vector(mode = "list", length = NFreq)
        biash_prev <- vector(mode = "list", length = NFreq)
        coefvar_prev <- vector(mode = "list", length = NFreq)

        for (i in 1:NFreq){ # For each choice of true frequency distribution
          amb_prevalence <- array(0, dim = c(numH, NEst))
          bias_prevalence <- array(0, dim = c(numH, NEst))
          coefvar_prevalence <- rep(0, numH)

          # Ambiguous prevalence
          ## True ambiguous prevalence
          tru_prev <- true_prev[[l]][[i]]

          ## Access each of the 10000 estimates
          tmp1 <- df_Estim[[l]][[k]][[j]][, , i]

          ## For each set of estimates, compute prevalence
          amb_prevalence <- (exp(tmp1[1,]) - exp(1-tmp1[2:numHpo,])^tmp1[1,])/(exp(tmp1[1,])-1)

          ## Replace entries with NAN values by 0
          amb_prevalence[is.na(amb_prevalence)] <- 0.0

          # Bias
          bias_prevalence <- (amb_prevalence/tru_prev[,j] - 1)*100
          bias_prevalence[is.infinite(bias_prevalence)] <- 0.0

          qh <- rowMeans(amb_prevalence, na.rm = TRUE)
          bias_prev <- rowMeans(bias_prevalence, na.rm = TRUE)

          # coefficient of variation
          for (q in 1:numH){
            coefvar_prevalence[q] <- sd(amb_prevalence[q,], na.rm = TRUE)/tru_prev[q,j]
          }

          ## Save the prevalence in a list
          qh_freq[[i]] <- qh
          biash_prev[[i]] <- bias_prev
          coefvar_prev[[i]] <- coefvar_prevalence
        }
        qh_lamb[[j]] <- qh_freq
        biash_lamb[[j]] <- biash_prev
        coefvar_lamb[[j]] <- coefvar_prev
      }
      qh_Samp[[k]] <- qh_lamb
      biash_Samp[[k]] <- biash_lamb
      coefvar_Samp[[k]] <- coefvar_lamb
    }
    qh_loc[[l]] <- qh_Samp
    biash_loc[[l]] <- biash_Samp
    coefvar_loc[[l]] <- coefvar_Samp
  }
  # Saving the estimates
  saveRDS(qh_loc, file = paste0(path, "dataset/ambPrevalenceEstimates", name, ".rds"))
  saveRDS(biash_loc, file = paste0(path, "dataset/biasAmbPrevalence", name, ".rds"))
  saveRDS(coefvar_loc, file = paste0(path, "dataset/coefvarAmbPrevalence", name, ".rds"))
}

perf_unamb_prevalence <- function(df_Estim, ParTru, true_prev, name){
 # This function implements the unambiguous prevalence as defined in the manuscript
 # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"
 
  qh_loc <- vector(mode = "list", length = Nn)
  biash_loc <- vector(mode = "list", length = Nn)
  coefvar_loc <- vector(mode = "list", length = Nn)
  for (l in 1:numbloci){   # For each number of locus
    qh_Samp <- vector(mode = "list", length = NN)
    biash_Samp <- vector(mode = "list", length = NN)
    coefvar_Samp <- vector(mode = "list", length = NN)
    tru_freq <- ParTru[[1]][[l]]

    for (k in 1:NN){  # For each true sample size
      qh_lamb <- vector(mode = "list", length = NFreq)
      biash_lamb <- vector(mode = "list", length = NFreq)
      coefvar_lamb <- vector(mode = "list", length = NFreq)

      # Number of haplotypes
      numH <- Hvec[l]

      # Number of loci
      nLoci <- NLoci[l]

      # Table of all possible haplotypes
      Hapl <- hapl(nLoci)

      ## For each haplotype in the table, build the set of observation Uh
      nLociUh <- nLoci+1

      for (j in 1:NLbd){ # For each true Lambda
        qh_freq <- vector(mode = "list", length = NFreq)
        biash_prev <- vector(mode = "list", length = NFreq)
        coefvar_prev <- vector(mode = "list", length = NFreq)

        for (i in 1:NFreq){ # For each choice of true frequency distribution
          rh <- rep(0, numH)
          bias_prev <- rep(0, numH)
          coefvar_prevalence <- rep(0, numH)
          tmp_prev <- array(0, dim = c(numH, NEst))

          # True prevalence
          tru_prev <- true_prev[[l]][[i]]

          ## Access each of the 10000 estimates
          tmp1 <- df_Estim[[l]][[k]][[j]][, , i]
          tmp2 <- tmp1[2:(numH+1),]

          ## For each haplotype build the set Uh for all l
          trufreq_vec <- tru_freq[i,]

          for (idx in 1:numH){ 
            if(trufreq_vec[idx] != 0){
              uh <- t(array(rep(Hapl[idx,], nLoci), dim = c(nLoci, nLociUh)))
              uh[2:nLociUh,] <- (uh[2:nLociUh,]+diag(nLoci))%%2 

              ## Pick the right frequencies estimates 
              pickh <- which(colSums(uh[1,] == t(Hapl))==nLoci)
              GPh <- gen_func(tmp2[pickh,], tmp1[1,])

              GPartFreq <- rep(0, NEst)
              GFreq <- rep(0, NEst)

              pick1 <- rep(0, nLociUh)
              pick2 <- rep(0, nLociUh)

              for(idxUh in 1:nLociUh){ 
                pick1[idxUh] <- which(colSums(uh[idxUh,] == t(Hapl))==nLoci)
              }
              pick2 <- pick1
              pick2[1] <- 0

              for(idxUh in 2:nLociUh){ 
                GPartFreq      <- gen_func(tmp2[pick2[idxUh],], tmp1[1,])
                GFreq          <- gen_func(colSums(tmp2[pick1[c(1,idxUh)],]), tmp1[1,])
                tmp_prev[idx,] <- tmp_prev[idx,] + GFreq - GPartFreq
              }
              estim_prev     <- tmp_prev[idx,] - (nLoci-1)*GPh
              rh[idx]        <- mean(estim_prev, na.rm = TRUE)

              # Bias
              bias           <- estim_prev/tru_prev[idx,j] - 1
              bias_prev[idx] <- mean(bias, na.rm = TRUE)*100

              # Coefficient of variation
              coefvar_prevalence[idx] <- sd(estim_prev, na.rm = TRUE)/tru_prev[idx,j]

            }else{
              rh[idx]        <- 0
              bias_prev[idx] <- 0
              coefvar_prevalence[idx] <- 0
            }
          }
          ## Save the prevalence in a list
          qh_freq[[i]]    <- rh
          biash_prev[[i]] <- bias_prev
          coefvar_prev[[i]] <- coefvar_prevalence
        }
        qh_lamb[[j]]      <- qh_freq
        biash_lamb[[j]]   <- biash_prev
        coefvar_lamb[[j]]   <- coefvar_prev
      }
      qh_Samp[[k]]        <- qh_lamb
      biash_Samp[[k]]     <- biash_lamb
      coefvar_Samp[[k]]     <- coefvar_lamb
    }
    qh_loc[[l]]    <- qh_Samp
    biash_loc[[l]] <- biash_Samp
    coefvar_loc[[l]] <- coefvar_Samp
  }

  # Saving the estimates
  saveRDS(qh_loc, file = paste0(path,  "dataset/unambPrevalenceEstimates", name, ".rds"))
  saveRDS(biash_loc, file = paste0(path,  "dataset/biasUnambPrevalence", name, ".rds"))
  saveRDS(coefvar_loc, file = paste0(path,  "dataset/coefvarUnambPrevalence", name, ".rds"))

  coefvar_loc
}

perf_prevalence <- function(df_Estim, ParTru, true_prev, name){
 # This function implements the unambiguous prevalence as defined in the manuscript
 # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"
 
  qh_loc <- vector(mode = "list", length = Nn)
  biash_loc <- vector(mode = "list", length = Nn)
  coefvar_loc <- vector(mode = "list", length = Nn)
  for (l in 1:numbloci){   # For each number of locus
    qh_Samp <- vector(mode = "list", length = NN)
    biash_Samp <- vector(mode = "list", length = NN)
    coefvar_Samp <- vector(mode = "list", length = NN)
    tru_freq <- ParTru[[1]][[l]]

    for (k in 1:NN){  # For each true sample size
      qh_lamb <- vector(mode = "list", length = NFreq)
      biash_lamb <- vector(mode = "list", length = NFreq)
      coefvar_lamb <- vector(mode = "list", length = NFreq)

      # Number of haplotypes
      numH <- Hvec[l]

      # Number of loci
      nLoci <- NLoci[l]

      # Table of all possible haplotypes
      Hapl <- hapl(nLoci)

      ## For each haplotype in the table, build the set of observation Uh
      nLociUh <- nLoci+1

      for (j in 1:NLbd){ # For each true Lambda
        qh_freq <- vector(mode = "list", length = NFreq)
        biash_prev <- vector(mode = "list", length = NFreq)
        coefvar_prev <- vector(mode = "list", length = NFreq)


        for (i in 1:NFreq){ # For each choice of true frequency distribution
          rh <- rep(0, numH)
          bias_prev <- rep(0, numH)
          coefvar_prevalence <- rep(0, numH)
          tmp_prevalence <- rep(0, NEst)
          tmp_prev <- array(0, dim = c(numH, NEst))

          # True prevalence
          tru_prev <- true_prev[[l]][[i]]

          ## Access each of the 10000 estimates
          tmp1 <- df_Estim[[l]][[k]][[j]][, , i]
          tmp2 <- tmp1[2:(numH+1),]

          ## For each haplotype build the set Uh for all l
          trufreq_vec <- tru_freq[i,]

          for (idx in 1:numH){ 
            if(trufreq_vec[idx] != 0){
              uh <- t(array(rep(Hapl[idx,], nLoci), dim = c(nLoci, nLociUh)))
              uh[2:nLociUh,] <- (uh[2:nLociUh,]+diag(nLoci))%%2 

              ## Pick the right frequencies estimates 
              pickh <- which(colSums(uh[1,] == t(Hapl))==nLoci)
              GPh <- gen_func(tmp2[pickh,], tmp1[1,])

              GPartFreq <- rep(0, NEst)
              GFreq <- rep(0, NEst)

              pick1 <- rep(0, nLociUh)
              pick2 <- rep(0, nLociUh)

              for(idxUh in 1:nLociUh){ 
                pick1[idxUh] <- which(colSums(uh[idxUh,] == t(Hapl))==nLoci)
              }
              pick2 <- pick1
              pick2[1] <- 0

              for(idxUh in 2:nLociUh){ 
                GPartFreq      <- gen_func(tmp2[pick2[idxUh],], tmp1[1,])
                GFreq          <- gen_func(colSums(tmp2[pick1[c(1,idxUh)],]), tmp1[1,])
                tmp_prev[idx,] <- tmp_prev[idx,] + GFreq - GPartFreq
              }
              estim_prev     <- tmp_prev - (nLoci-1)*GPh
              for (q in 1:NEst){
                tmp_prevalence[q] <- estim_prev[,q]/sum(estim_prev[,q])
              }
              rh[idx]        <- mean(tmp_prevalence, na.rm = TRUE)

              # Bias
              bias <- tmp_prevalence/tru_prev[idx,j] - 1
              bias_prev[idx] <- mean(bias, na.rm = TRUE)*100

              # Coefficient of variation
              coefvar_prevalence[idx] <- sd(tmp_prevalence, na.rm = TRUE)/tru_prev[idx,j]
            }else{
              rh[idx]        <- 0
              bias[idx]      <- 0
              coefvar_prevalence[idx]  <- 0
            }
          }

          ## Save the prevalence in a list
          qh_freq[[i]]    <- rh
          biash_prev[[i]] <- bias
          coefvar_prev[[i]] <- coefvar_prevalence
        }
        qh_lamb[[j]]      <- qh_freq
        biash_lamb[[j]]   <- biash_prev
        coefvar_lamb[[j]]   <- coefvar_prev
      }
      qh_Samp[[k]]        <- qh_lamb
      biash_Samp[[k]]     <- biash_lamb
      coefvar_Samp[[k]]     <- coefvar_lamb
    }
    qh_loc[[l]]    <- qh_Samp
    biash_loc[[l]] <- biash_Samp
    coefvar_loc[[l]] <- coefvar_Samp
  }

  # Saving the estimates
  #saveRDS(qh_loc, file = paste0(path,  "dataset/PrevalenceEstimates", name, ".rds"))
  #saveRDS(biash_loc, file = paste0(path,  "dataset/biasPrevalence", name, ".rds"))
  #saveRDS(coefvar_loc, file = paste0(path,  "dataset/coefvarPrevalence", name, ".rds"))
  coefvar_loc
}

main <- function(df_Param, true_Par, name){
  # Loading estimated haplotype frequencies and MOI
  df_Estim <- readRDS(paste0(path, "dataset/modelEstimates", name, ".rds"))

  # Bias of frequencies and MOI
  bias(df_Estim, df_Param, name)

  # Coefficient of variation of MOI
  coefvar(df_Estim, df_Param, name)

  # True ambiguous prevalence
  tru_AmbPrev   <- true_amb_prevalence(true_Param, name)

  # True unambiguous prevalence
  tru_UnambPrev <- true_unamb_prevalence(true_Param, dfParam, name)

  # True prevalence
  tru_Prev      <- true_prevalence(true_Param, dfParam, name)

  # Estimated ambiguous prevalence
  perf_amb_prevalence(df_Estim, tru_AmbPrev, name)

  # Estimated unambiguous prevalence
  perf_unamb_prevalence(df_Estim, df_Param, tru_UnambPrev, name)

  # Estimated prevalence
  perf_prevalence(df_Estim, df_Param, tru_Prev, name)

}

path <- "/Volumes/GoogleDrive-117934057836063832284/My Drive/Maths against Malaria/Christian/Models/MultiLociBiallelicModel/"

# Define data origin
name <- ''

# Loading true haplotype frequencies and MOI
dfParam <- readRDS(paste0(path, "dataset/trueParameters", name, ".rds"))

# Loading extra parameters
parExtr  <- readRDS(paste0(path, "dataset/extraParameters", name, ".rds"))

# Variables initialization
NLbd          <- parExtr[[1]]
Nn            <- parExtr[[2]]
Hvec          <- parExtr[[3]]
NN            <- parExtr[[4]]
NEst          <- parExtr[[5]]
NFreq         <- parExtr[[6]]
NLoci         <- log2(Hvec)
mean_moi_true <- psi(dfParam[[2]])

numbloci <- length(Hvec)

# Reformatting true parameters to compute true prevalence
true_Param <- vector(mode='list', length=Nn)

for (i in 1:Nn){
  tot_row <- Hvec[i]+1
  true_Param[[i]] <- vector(mode='list', length=NFreq)
  for (j in 1:NFreq){
    true_Param[[i]][[j]] <- array(0, c(tot_row, NLbd))
    true_Param[[i]][[j]][1,] <- dfParam[[2]]
    true_Param[[i]][[j]][2:tot_row,] <- dfParam[[1]][[i]][j,]
  }
}

# Running the performance checker ('' <- simulated data, 'Kenya' <- kenyan data)
main(dfParam, true_Param, name)