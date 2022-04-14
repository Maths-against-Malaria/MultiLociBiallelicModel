# Title        : Simulation study
# Objective    : Compute the bias and coefficient of variation of the estimates,
#                and the estimates of prevalence
# Created by   : Christian Tsoungui Obama
# Created on   : 03.04.21
# Last modified: 06.03.22

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

bias     <- function(estim_Param, sim_Param, name){
  # This function implements bias of mean MOI and haplotype frequencies as defined in the manuscript
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence 
  # alongside multiplicity of infection from SNPs data"

  moibiasloc <- vector(mode = "list", length = n_Sim_Loci)
  freqbiasloc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    freqbiasSamp <- vector(mode = "list", length = n_Sampl)
    moibiasSamp  <- vector(mode = "list", length = n_Sampl)

    # Number of haplotypes
    num_Hapl   <- n_Hapl[l]
    num_Hapl_PlusOne <- num_Hapl + 1

    for (k in 1:n_Sampl){  # For each sample size
      freqbias_lamb <- vector(mode = "list", length = n_Freq_Distr)
      moibias_lamb  <- vector(mode = "list", length = n_Freq_Distr)

      for (j in 1:n_Lbda){ # For each true Lambda
        freq_bias <- vector(mode = "list", length = n_Freq_Distr)
        moi_bias  <- vector(mode = "list", length = n_Freq_Distr)

        for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
          # Access each of the 10000 estimates
          tmp1 <- estim_Param[[l]][[k]][[j]][, , i]

          # True frequencies
          tmp2 <- t(sim_Param[[1]][[l]])
          tmp3 <- tmp2[, i]

          # Relative bias of haplotype freq. in percent
          freqbias <- rowMeans((tmp1[2:num_Hapl_PlusOne,]/tmp3 - 1), na.rm = TRUE) * 100

          # Remove the estimates of Lambda that yield an "Infinite mean MOI"
          mean_moi_estim <- psi(tmp1[1, ])
          
          # Relative bias of MOI in percent
          bias_moi <- (mean_moi_estim/true_Mean_MOI[j] - 1)
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

coefvar  <- function(estim_Param, sim_Param, name){
  # This function implements coefficient of variation of mean MOI as defined in the manuscript
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence 
  # alongside multiplicity of infection from SNPs data"

  moicvloc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    moicvSamp  <- vector(mode = "list", length = n_Sampl)

    # Number of haplotypes
    num_Hapl   <- n_Hapl[l]
    num_Hapl_PlusOne <- num_Hapl + 1

    for (k in 1:n_Sampl){  # For each sample size
      moicv_lamb  <- vector(mode = "list", length = n_Freq_Distr)

      for (j in 1:n_Lbda){ # For each true Lambda
        moi_cv  <- vector(mode = "list", length = n_Freq_Distr)

        for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
          # Access each of the 10000 estimates
          tmp1 <- estim_Param[[l]][[k]][[j]][, , i]

          # Remove the estimates of Lambda that yield an "Infinite mean MOI"
          mean_moi_estim <- psi(tmp1[1, ])
          mean_moi_estim <- mean_moi_estim[!is.infinite(mean_moi_estim)]
          
          # Dimensionless coefficient of variation of MOI
          moicv <- sd(mean_moi_estim, na.rm = TRUE) / true_Mean_MOI[j]

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

true_amb_prevalence         <- function(reshap_Sim_Param, name){
  # This function implements the true ambiguous prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"

  qh_loc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){
    qh_loc[[l]] <- vector(mode = "list", length = n_Freq_Distr)

    # Number of haplotypes
    num_Hapl         <- n_Hapl[l]
    num_Hapl_PlusOne <- n_Hapl[l] + 1

    for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
      amb_prevalence <- array(0, dim = c(num_Hapl, n_Lbda))

      # Access each of the 10.000 parameters estimates
      tmp1 <- reshap_Sim_Param[[l]][[i]]

      # For each set of parameters estimates
      for (j in 1:n_Lbda){
        amb_prevalence[,j] <- (exp(tmp1[1,j]) - exp(1-tmp1[2:num_Hapl_PlusOne,j])^tmp1[1,j])/(exp(tmp1[1,j])-1)
      }

      qh_loc[[l]][[i]] <- amb_prevalence
    }
  }

  # Save the ambiguous prevalence estimates
  saveRDS(qh_loc, file = paste0(path, "dataset/true_Amb_Prevalence", name, ".rds"))
  qh_loc
}

true_unamb_prevalence       <- function(reshap_Sim_Param, sim_Param, name){
  # This function implements the true unambiguous prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"
 
  qh_loc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    qh_loc[[l]] <- vector(mode = "list", length = n_Freq_Distr)
    true_freq   <- sim_Param[[1]][[l]]

    # Number of haplotypes
    num_Hapl <- n_Hapl[l]

    # Number of loci
    numb_Loci <- n_Loci[l]

    # Table of all possible haplotypes
    Hapl <- hapl(numb_Loci)

    ## For each haplotype in the table, build the set of observation Uh
    numb_Hapl_Uh <- numb_Loci + 1

    for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
      rh <- array(0, dim = c(num_Hapl, n_Lbda))

      ## For each combination of true lanmbda and haplotype frequencies values
      tmp1 <- reshap_Sim_Param[[l]][[i]]
      tmp2 <- tmp1[2:(num_Hapl+1),]

      ## For each haplotype build the set Uh for all l
      trufreq_vec <- true_freq[i,]
      pickhap     <- which(trufreq_vec != 0)

      for (idx in pickhap){ 
        uh                  <- t(array(rep(Hapl[idx,], numb_Loci), dim=c(numb_Loci, numb_Hapl_Uh)))
        uh[2:numb_Hapl_Uh,] <- (uh[2:numb_Hapl_Uh,]+diag(numb_Loci))%%2

        ## Pick the right frequencies estimates 
        pickh <- which(colSums(uh[1,] == t(Hapl))==numb_Loci)
        GPh   <- gen_func(tmp2[pickh,], tmp1[1,])

        GPartFreq <- rep(0, n_Lbda)
        GFreq    <- rep(0, n_Lbda)

        pick1 <- rep(0, numb_Hapl_Uh)
        pick2 <- rep(0, numb_Hapl_Uh)

        for(idxUh in 1:numb_Hapl_Uh){ 
          pick1[idxUh] <- which(colSums(uh[idxUh,] == t(Hapl))==numb_Loci)
        }

        pick2    <- pick1
        pick2[1] <- 0

        for(idxUh in 2:numb_Hapl_Uh){ 
          GPartFreq <- gen_func(tmp2[pick2[idxUh],], tmp1[1,])
          GFreq     <- gen_func(colSums(tmp2[pick1[c(1,idxUh)],]), tmp1[1,])
          rh[idx,]  <- rh[idx,] + GFreq - GPartFreq
        }
        rh[idx,] <- rh[idx,] - (numb_Loci - 1)*GPh
      }
      qh_loc[[l]][[i]] <- rh
    }
  }
  # Save the unambiguous prevalence estimates
  saveRDS(qh_loc, file = paste0(path,  "dataset/true_Unamb_Prevalence", name, ".rds"))
  qh_loc
}

true_conditional_prevalence <- function(reshap_Sim_Param, sim_Param, name){
  # This function implements the true conditional prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"
 
  qh_loc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    qh_loc[[l]] <- vector(mode = "list", length = n_Freq_Distr)
    true_freq <- sim_Param[[1]][[l]]

    # Number of haplotypes
    num_Hapl <- n_Hapl[l]

    # Number of loci
    numb_Loci <- n_Loci[l]

    # Table of all possible haplotypes
    Hapl <- hapl(numb_Loci)

    ## For each haplotype in the table, build the set of observation Uh
    numb_Hapl_Uh <- numb_Loci + 1

    for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
      numh <- array(0, dim = c(num_Hapl, n_Lbda))
      denh <- array(0, dim = c(num_Hapl, n_Lbda))
      rh   <- array(0, dim = c(num_Hapl, n_Lbda))

      ## For each combination of true lanmbda and haplotype frequencies values
      tmp1 <- reshap_Sim_Param[[l]][[i]]
      tmp2 <- tmp1[2:(num_Hapl+1),]

      ## For eachobserved haplotype build the set Uh for all l
      trufreq_vec <- true_freq[i,]
      pickhap     <- which(trufreq_vec != 0)

      for (idx in pickhap){ 
        uh                  <- t(array(rep(Hapl[idx,], numb_Loci), dim=c(numb_Loci, numb_Hapl_Uh)))
        uh[2:numb_Hapl_Uh,] <- (uh[2:numb_Hapl_Uh,]+diag(numb_Loci))%%2

        ## Pick the right frequencies estimates 
        pickh <- which(colSums(uh[1,] == t(Hapl))==numb_Loci)
        GPh   <- gen_func(tmp2[pickh,], tmp1[1,])

        GPartFreq <- rep(0, n_Lbda)
        GFreq     <- rep(0, n_Lbda)

        pick1 <- rep(0, numb_Hapl_Uh)
        pick2 <- rep(0, numb_Hapl_Uh)

        for(idxUh in 1:numb_Hapl_Uh){ 
          pick1[idxUh] <- which(colSums(uh[idxUh,] == t(Hapl))==numb_Loci)
        }

        pick2    <- pick1
        pick2[1] <- 0

        for(idxUh in 2:numb_Hapl_Uh){ 
          GPartFreq  <- gen_func(tmp2[pick2[idxUh],], tmp1[1,])
          GFreq      <- gen_func(colSums(tmp2[pick1[c(1,idxUh)],]), tmp1[1,])
          tmp        <- GFreq - GPartFreq
          numh[idx,] <- numh[idx,] + tmp
          denh[idx,] <- denh[idx,] + tmp/2
        }
          numh[idx,] <- numh[idx,] - (numb_Loci - 1)*GPh
          denh[idx,] <- denh[idx,] - (numb_Loci/2 - 1)*GPh
      }

      den <- colSums(denh, na.rm = TRUE)

      for(q in 1:n_Lbda){
        rh[,q] <- numh[,q]/den[q]
      }
      qh_loc[[l]][[i]] <- rh
    }
  }  
  # Saving the cunditional prevalence estimates
  saveRDS(qh_loc, file = paste0(path,  "dataset/true_Cond_Prevalence", name, ".rds"))
  qh_loc
}

true_relative_prevalence    <- function(reshap_Sim_Param, sim_Param, name){
  # This function implements the true relative prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"
 
  qh_loc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    qh_loc[[l]] <- vector(mode = "list", length = n_Freq_Distr)
    true_freq   <- sim_Param[[1]][[l]]

    # Number of haplotypes
    num_Hapl <- n_Hapl[l]

    # Number of loci
    numb_Loci <- n_Loci[l]

    # Table of all possible haplotypes
    Hapl <- hapl(numb_Loci)

    ## For each haplotype in the table, build the set of observation Uh
    numb_Hapl_Uh <- numb_Loci + 1

      for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
        rh <- array(0, dim = c(num_Hapl, n_Lbda))

        ## For each combination of true lanmbda and haplotype frequencies values
        tmp1 <- reshap_Sim_Param[[l]][[i]]
        tmp2 <- tmp1[2:(num_Hapl+1),]

        ## For each haplotype build the set Uh for all l
        trufreq_vec <- true_freq[i,]
        pickhap <- which(trufreq_vec != 0)

        for (idx in pickhap){ 
          uh                  <- t(array(rep(Hapl[idx,], numb_Loci), dim=c(numb_Loci, numb_Hapl_Uh)))
          uh[2:numb_Hapl_Uh,] <- (uh[2:numb_Hapl_Uh,]+diag(numb_Loci))%%2

          ## Pick the right frequencies estimates 
          pickh <- which(colSums(uh[1,] == t(Hapl))==numb_Loci)
          GPh   <- gen_func(tmp2[pickh,], tmp1[1,])

          GPartFreq <- rep(0, n_Lbda)
          GFreq     <- rep(0, n_Lbda)

          pick1 <- rep(0, numb_Hapl_Uh)
          pick2 <- rep(0, numb_Hapl_Uh)

          for(idxUh in 1:numb_Hapl_Uh){ 
            pick1[idxUh] <- which(colSums(uh[idxUh,] == t(Hapl))==numb_Loci)
          }

          pick2    <- pick1
          pick2[1] <- 0

          for(idxUh in 2:numb_Hapl_Uh){ 
            GPartFreq <- gen_func(tmp2[pick2[idxUh],], tmp1[1,])
            GFreq     <- gen_func(colSums(tmp2[pick1[c(1,idxUh)],]), tmp1[1,])
            rh[idx,]  <- rh[idx,] + GFreq - GPartFreq
          }
            rh[idx,]  <- rh[idx,] - (numb_Loci - 1)*GPh
        }
      for (q in 1:n_Lbda){
        rh[,q] <- rh[,q]/sum(rh[,q], na.rm = TRUE)
      }
      qh_loc[[l]][[i]] <- rh
    }
  }
  # Saving the relative prevalence estimates
  saveRDS(qh_loc, file = paste0(path,  "dataset/true_Rel_Prevalence", name, ".rds"))
  qh_loc
}


estim_amb_prevalence <- function(estim_Param, true_prev, name){
  # This function estimates the ambiguous prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"

  qh_loc <- vector(mode = "list", length = n_Sim_Loci)
  #biash_loc <- vector(mode = "list", length = n_Sim_Loci)
  #coefvar_loc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    qh_Samp      <- vector(mode = "list", length = n_Sampl)
    #biash_Samp   <- vector(mode = "list", length = n_Sampl)
    #coefvar_Samp <- vector(mode = "list", length = n_Sampl)

    # Number of haplotypes
    num_Hapl         <- n_Hapl[l]
    num_Hapl_PlusOne <- num_Hapl + 1

    for (k in 1:n_Sampl){  # For each true sample size
      qh_lamb <- vector(mode = "list", length = n_Freq_Distr)
      #biash_lamb <- vector(mode = "list", length = n_Freq_Distr)
      #coefvar_lamb <- vector(mode = "list", length = n_Freq_Distr)

      for (j in 1:n_Lbda){ # For each true Lambda
        qh_freq <- vector(mode = "list", length = n_Freq_Distr)
        #biash_prev <- vector(mode = "list", length = n_Freq_Distr)
        #coefvar_prev <- vector(mode = "list", length = n_Freq_Distr)

        for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
          amb_prevalence <- array(0, dim = c(num_Hapl, n_Sampl_Gen))
          #bias_prevalence <- array(0, dim = c(num_Hapl, n_Sampl_Gen))
          #coefvar_prevalence <- rep(0, num_Hapl)

          # Ambiguous prevalence
          ## True ambiguous prevalence
          #true_prev <- true_prev[[l]][[i]]

          ## Access each of the 10000 estimates
          tmp1 <- estim_Param[[l]][[k]][[j]][, , i]

          ## For each set of estimates, compute prevalence
          amb_prevalence <- (exp(tmp1[1,]) - exp(1-tmp1[2:num_Hapl_PlusOne,])^tmp1[1,])/(exp(tmp1[1,])-1)

          ## Replace entries with NAN values by 0
          amb_prevalence[is.na(amb_prevalence)] <- 0.0

          # Bias
          #bias_prevalence <- (amb_prevalence/true_prev[,j] - 1)*100
          #bias_prevalence[is.infinite(bias_prevalence)] <- 0.0

          qh <- rowMeans(amb_prevalence, na.rm = TRUE)
          #bias_prev <- rowMeans(bias_prevalence, na.rm = TRUE)

          # coefficient of variation
          #for (q in 1:num_Hapl){
          #  coefvar_prevalence[q] <- sd(amb_prevalence[q,], na.rm = TRUE)/true_prev[q,j]
          #}

          ## Save the prevalence in a list
          qh_freq[[i]] <- qh
          #biash_prev[[i]] <- bias_prev
          #coefvar_prev[[i]] <- coefvar_prevalence
        }
        qh_lamb[[j]] <- qh_freq
        #biash_lamb[[j]] <- biash_prev
        #coefvar_lamb[[j]] <- coefvar_prev
      }
      qh_Samp[[k]] <- qh_lamb
      #biash_Samp[[k]] <- biash_lamb
      #coefvar_Samp[[k]] <- coefvar_lamb
    }
    qh_loc[[l]] <- qh_Samp
    #biash_loc[[l]] <- biash_Samp
    #coefvar_loc[[l]] <- coefvar_Samp
  }
  # Saving the estimates
  saveRDS(qh_loc, file = paste0(path, "dataset/estim_Amb_Prevalence", name, ".rds"))
  #saveRDS(biash_loc, file = paste0(path, "dataset/biasAmbPrevalence", name, ".rds"))
  #saveRDS(coefvar_loc, file = paste0(path, "dataset/coefvarAmbPrevalence", name, ".rds"))
  qh_loc
}

estim_unamb_prevalence <- function(estim_Param, sim_Param, true_prev, name){
  # This function estimates the unambiguous prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"

  qh_loc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    qh_Samp <- vector(mode = "list", length = n_Sampl)

    # True haplotype frequencies
    true_freq <- sim_Param[[1]][[l]]

    for (k in 1:n_Sampl){  # For each true sample size
      qh_lamb <- vector(mode = "list", length = n_Freq_Distr)

      # Number of haplotypes
      num_Hapl <- n_Hapl[l]

      # Number of loci
      numb_Loci <- n_Loci[l]

      # Table of all possible haplotypes
      Hapl <- hapl(numb_Loci)

      ## For each haplotype in the table, build the set of observation Uh
      numb_Hapl_Uh <- numb_Loci + 1

      for (j in 1:n_Lbda){ # For each true Lambda
        qh_freq <- vector(mode = "list", length = n_Freq_Distr)
        #biash_prev <- vector(mode = "list", length = n_Freq_Distr)
        #coefvar_prev <- vector(mode = "list", length = n_Freq_Distr)

        for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
          rh <- rep(0, num_Hapl)
          #bias_prev <- rep(0, num_Hapl)
          #coefvar_prevalence <- rep(0, num_Hapl)
          tmp_prev <- array(0, dim = c(num_Hapl, n_Sampl_Gen))

          # True prevalence
          #true_prev <- true_prev[[l]][[i]]

          ## Access each of the 10000 estimates
          tmp1 <- estim_Param[[l]][[k]][[j]][, , i]
          tmp2 <- tmp1[2:(num_Hapl+1),]

          ## For each haplotype build the set Uh for all l
          trufreq_vec <- true_freq[i,]
          pickhap     <- which(trufreq_vec != 0)

          for (idx in pickhap){ 
            uh                  <- t(array(rep(Hapl[idx,], numb_Loci), dim = c(numb_Loci, numb_Hapl_Uh)))
            uh[2:numb_Hapl_Uh,] <- (uh[2:numb_Hapl_Uh,]+diag(numb_Loci))%%2 

            ## Pick the right frequencies estimates 
            pickh <- which(colSums(uh[1,] == t(Hapl))==numb_Loci)
            GPh   <- gen_func(tmp2[pickh,], tmp1[1,])

            GPartFreq <- rep(0, n_Sampl_Gen)
            GFreq     <- rep(0, n_Sampl_Gen)

            pick1 <- rep(0, numb_Hapl_Uh)
            pick2 <- rep(0, numb_Hapl_Uh)

            for(idxUh in 1:numb_Hapl_Uh){ 
              pick1[idxUh] <- which(colSums(uh[idxUh,] == t(Hapl))==numb_Loci)
            }

            pick2    <- pick1
            pick2[1] <- 0

            for(idxUh in 2:numb_Hapl_Uh){ 
              GPartFreq      <- gen_func(tmp2[pick2[idxUh],], tmp1[1,])
              GFreq          <- gen_func(colSums(tmp2[pick1[c(1,idxUh)],]), tmp1[1,])
              tmp_prev[idx,] <- tmp_prev[idx,] + GFreq - GPartFreq
            }
            estim_prev     <- tmp_prev[idx,] - (numb_Loci - 1)*GPh
            rh[idx]        <- mean(estim_prev, na.rm = TRUE)

            # Bias
            #bias           <- estim_prev/true_prev[idx,j] - 1
            #bias_prev[idx] <- mean(bias, na.rm = TRUE)*100

            # Coefficient of variation
            #coefvar_prevalence[idx] <- sd(estim_prev, na.rm = TRUE)/true_prev[idx,j]
          }
          ## Save the prevalence in a list
          qh_freq[[i]]    <- rh
          #biash_prev[[i]] <- bias_prev
          #coefvar_prev[[i]] <- coefvar_prevalence
        }
        qh_lamb[[j]]      <- qh_freq
        #biash_lamb[[j]]   <- biash_prev
        #coefvar_lamb[[j]]   <- coefvar_prev
      }
      qh_Samp[[k]]        <- qh_lamb
      #biash_Samp[[k]]     <- biash_lamb
      #coefvar_Samp[[k]]     <- coefvar_lamb
    }
    qh_loc[[l]]    <- qh_Samp
    #biash_loc[[l]] <- biash_Samp
    #coefvar_loc[[l]] <- coefvar_Samp
  }
  # Saving the estimates of unambiguous prevalence
  saveRDS(qh_loc, file = paste0(path,  "dataset/estim_Unamb_Prevalence", name, ".rds"))
  #saveRDS(biash_loc, file = paste0(path,  "dataset/biasUnambPrevalence", name, ".rds"))
  #saveRDS(coefvar_loc, file = paste0(path,  "dataset/coefvarUnambPrevalence", name, ".rds"))
  qh_loc
}

estim_conditional_prevalence <- function(estim_Param, sim_Param, true_prev, name){
  # This function estimates the unambiguous prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"
 
  qh_loc <- vector(mode = "list", length = n_Sim_Loci)
  #biash_loc   <- vector(mode = "list", length = n_Sim_Loci)
  #coefvar_loc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    qh_Samp <- vector(mode = "list", length = n_Sampl)
    #biash_Samp   <- vector(mode = "list", length = n_Sampl)
    #coefvar_Samp <- vector(mode = "list", length = n_Sampl)

    # True haplotype frequencies
    true_freq <- sim_Param[[1]][[l]]

    for (k in 1:n_Sampl){  # For each true sample size
      qh_lamb <- vector(mode = "list", length = n_Freq_Distr)
      #biash_lamb   <- vector(mode = "list", length = n_Freq_Distr)
      #coefvar_lamb <- vector(mode = "list", length = n_Freq_Distr)

      # Number of haplotypes
      num_Hapl <- n_Hapl[l]

      # Number of loci
      numb_Loci <- n_Loci[l]

      # Table of all possible haplotypes
      Hapl <- hapl(numb_Loci)

      ## For each haplotype in the table, build the set of observation Uh
      numb_Hapl_Uh <- numb_Loci + 1

      for (j in 1:n_Lbda){ # For each true Lambda
        qh_freq      <- vector(mode = "list", length = n_Freq_Distr)
        #bias_prev    <- vector(mode = "list", length = n_Freq_Distr)
        #coefvar_prev <- vector(mode = "list", length = n_Freq_Distr)

        for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
          rh   <- array(0, dim = c(num_Hapl, n_Sampl_Gen))
          numh <- array(0, dim = c(num_Hapl, n_Sampl_Gen))
          denh <- array(0, dim = c(num_Hapl, n_Sampl_Gen))
          prev <- array(0, dim = c(num_Hapl, n_Sampl_Gen))

          prevalence <- rep(0, num_Hapl)
          #bias_prevalence    <- rep(0, num_Hapl)
          #coefvar_prevalence <- rep(0, num_Hapl)

          ## Access each of the 10000 estimates
          tmp1 <- estim_Param[[l]][[k]][[j]][, , i]
          tmp2 <- tmp1[2:(num_Hapl+1),]

          ## For each haplotype build the set Uh for all l
          trufreq_vec <- true_freq[i,]
          pickhap <- which(trufreq_vec != 0)

          # Find ambiguous prevalence for each haplotype
          for (idx in pickhap){ 
            uh                  <- t(array(rep(Hapl[idx,], numb_Loci), dim=c(numb_Loci, numb_Hapl_Uh)))
            uh[2:numb_Hapl_Uh,] <- (uh[2:numb_Hapl_Uh,]+diag(numb_Loci))%%2

            ## Pick the right frequencies estimates 
            pickh <- which(colSums(uh[1,] == t(Hapl))==numb_Loci)
            GPh   <- gen_func(tmp2[pickh,], tmp1[1,])

            GPartFreq <- rep(0, n_Lbda)
            GFreq     <- rep(0, n_Lbda)

            pick1 <- rep(0, numb_Hapl_Uh)
            pick2 <- rep(0, numb_Hapl_Uh)

            for(idxUh in 1:numb_Hapl_Uh){ 
              pick1[idxUh] <- which(colSums(uh[idxUh,] == t(Hapl))==numb_Loci)
            }

            pick2 <- pick1
            pick2[1] <- 0

            for(idxUh in 2:numb_Hapl_Uh){ 
              GPartFreq  <- gen_func(tmp2[pick2[idxUh],], tmp1[1,])
              GFreq      <- gen_func(colSums(tmp2[pick1[c(1,idxUh)],]), tmp1[1,])
              tmp        <- GFreq - GPartFreq
              numh[idx,] <- numh[idx,] + tmp
              denh[idx,] <- denh[idx,] + tmp/2
            }
              numh[idx,] <- numh[idx,] - (numb_Loci - 1)*GPh
              denh[idx,] <- denh[idx,] - (numb_Loci/2 - 1)*GPh
          }

          den <- colSums(denh, na.rm = TRUE)
          for (q in 1:n_Sampl_Gen){
            prev[,q] <- numh[,q]/den[q]
          }

          prevalence <- rowMeans(prev, na.rm = TRUE)
          
          #for (q in 1:num_Hapl){
            # Conditional prevalence
            # prevalence[q] <- mean(prev[q,], na.rm = TRUE)

            # Bias
            # bias <- prev[q,]/true_prev[q,j] - 1
            # bias_prevalence[q] <- mean(bias, na.rm = TRUE)*100

            # Coefficient of variation
            # coefvar_prevalence[q] <- sd(prev[q,], na.rm = TRUE)/true_prev[q,j]
          #}
  
          ## Save the prevalence in a list
          qh_freq[[i]]  <- prevalence
         # bias_prev[[i]]    <- bias_prevalence
         # coefvar_prev[[i]] <- coefvar_prevalence
        }
        qh_lamb[[j]] <- qh_freq
        #biash_lamb[[j]]   <- bias_prev
        #coefvar_lamb[[j]] <- coefvar_prev
      }
      qh_Samp[[k]]   <- qh_lamb
      #biash_Samp[[k]]     <- biash_lamb
      #coefvar_Samp[[k]]   <- coefvar_lamb
    }
    qh_loc[[l]]      <- qh_Samp
    #biash_loc[[l]]   <- biash_Samp
    #coefvar_loc[[l]] <- coefvar_Samp
  }
  # Saving the estimates
  saveRDS(qh_loc, file = paste0(path, "dataset/estim_Cond_Prevalence", name, ".rds"))
  #saveRDS(biash_loc, file = paste0(path,  "dataset/biasConditionalPrevalence", name, ".rds"))
  #saveRDS(coefvar_loc, file = paste0(path,  "dataset/coefvarConditionalPrevalence", name, ".rds"))
  qh_loc
}

estim_relative_prevalence <- function(estim, name){
  # This function estimates the relative prevalence as defined in the manuscript of tsoungui et.al, titled
  # "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data"

  qh_loc <- vector(mode = "list", length = n_Sim_Loci)
  #biash_loc <- vector(mode = "list", length = n_Sim_Loci)
  #coefvar_loc <- vector(mode = "list", length = n_Sim_Loci)

  for (l in 1:n_Sim_Loci){   # For each number of locus
    qh_Samp <- vector(mode = "list", length = n_Sampl)
    #biash_Samp <- vector(mode = "list", length = n_Sampl)
    #coefvar_Samp <- vector(mode = "list", length = n_Sampl)

    # True haplotype frequencies
    true_freq <- sim_Param[[1]][[l]]

    for (k in 1:n_Sampl){  # For each true sample size
      qh_lamb <- vector(mode = "list", length = n_Freq_Distr)
      #biash_lamb <- vector(mode = "list", length = n_Freq_Distr)
      #coefvar_lamb <- vector(mode = "list", length = n_Freq_Distr)

      # Number of haplotypes
      num_Hapl <- n_Hapl[l]

      # Number of loci
      numb_Loci <- n_Loci[l]

      # Table of all possible haplotypes
      Hapl <- hapl(numb_Loci)

      ## For each haplotype in the table, build the set of observation Uh
      numb_Hapl_Uh <- numb_Loci + 1

      for (j in 1:n_Lbda){ # For each true Lambda
        qh_freq <- vector(mode = "list", length = n_Freq_Distr)
        #bias_prev <- vector(mode = "list", length = n_Freq_Distr)
        #coefvar_prev <- vector(mode = "list", length = n_Freq_Distr)

        for (i in 1:n_Freq_Distr){ # For each choice of true frequency distribution
          prev       <- array(0, dim = c(num_Hapl, n_Sampl_Gen))
          prevalence <- rep(0, num_Hapl)
          #bias_prevalence <- rep(0, num_Hapl)
          #coefvar_prevalence <- rep(0, num_Hapl)
          #tmp_prev <- array(0, dim = c(num_Hapl, n_Sampl_Gen))

          # True prevalence
          #true_prev <- true_prev[[l]][[i]]

          ## Access each of the 10000 estimates
          prev <- estim[[l]][[k]][[j]][, , i]
          prevalence <- rowMeans(prev, na.rm = TRUE)
          
          #for (q in 1:num_Hapl){
            # Relative prevalence
            #prevalence[q] <- mean(prev[q,], na.rm = TRUE)

            # Bias
            #bias <- prev[q,]/true_prev[q,j] - 1
            #bias_prevalence[q] <- mean(bias, na.rm = TRUE)*100

            # Coefficient of variation
            #coefvar_prevalence[q] <- sd(prev[q,], na.rm = TRUE)/true_prev[q,j]
          #}
  
          ## Save the prevalence in a list
          qh_freq[[i]]      <- prevalence
          #bias_prev[[i]]    <- bias_prevalence
          #coefvar_prev[[i]] <- coefvar_prevalence
        }
        qh_lamb[[j]]      <- qh_freq
        #biash_lamb[[j]]   <- bias_prev
        #coefvar_lamb[[j]] <- coefvar_prev
      }
      qh_Samp[[k]]        <- qh_lamb
      #biash_Samp[[k]]     <- biash_lamb
      #coefvar_Samp[[k]]   <- coefvar_lamb
    }
    qh_loc[[l]]      <- qh_Samp
    #biash_loc[[l]]   <- biash_Samp
    #coefvar_loc[[l]] <- coefvar_Samp
  }

  # Saving the estimates
  saveRDS(qh_loc, file = paste0(path,  "dataset/estim_Rel_Prevalence", name, ".rds"))
  #saveRDS(biash_loc, file = paste0(path,  "dataset/biasRelativePrevalence", name, ".rds"))
  #saveRDS(coefvar_loc, file = paste0(path,  "dataset/coefvarRelativePrevalence", name, ".rds"))
  qh_loc
}

main <- function(sim_Param, reshap_Sim_Param, name){
  # Loading estimated haplotype frequencies and MOI
  estim_Param       <- readRDS(paste0(path, "dataset/modelEstimates", name, ".rds"))
  adhoc_estim_Param <- readRDS(paste0(path, "dataset/adhocModelEstimates", name, ".rds"))

  # Bias of frequencies and MOI
  bias(estim_Param, sim_Param, name)

  # Coefficient of variation of MOI
  coefvar(estim_Param, sim_Param, name)

  # True ambiguous prevalence
  true_Amb_Prev          <- true_amb_prevalence(reshap_Sim_Param, name)

  # True unambiguous prevalence
  true_Unamb_Prev        <- true_unamb_prevalence(reshap_Sim_Param, sim_Param, name)

  # True conditional prevalence
  true_Conditional_Prev  <- true_conditional_prevalence(reshap_Sim_Param, sim_Param, name)

  # True relative prevalence
  true_Relative_Prev     <- true_relative_prevalence(reshap_Sim_Param, sim_Param, name)


  # Estimated ambiguous prevalence
  estim_amb_prevalence(estim_Param, true_Amb_Prev, name)

  # Estimated unambiguous prevalence
  estim_unamb_prevalence(estim_Param, sim_Param, true_Unamb_Prev, name)

  # Estimated conditional prevalence
  estim_conditional_prevalence(estim_Param, sim_Param, true_Conditional_Prev, name)

  # Estimated relative prevalence (adhoc Model)
  estim_relative_prevalence(adhoc_estim_Param, name)
}
 
path <- "/Volumes/GoogleDrive-117934057836063832284/My Drive/Maths against Malaria/Christian/Models/MultiLociBiallelicModel/"

# Define data origin ('' <- simulated data, 'Kenya' <- kenyan data)
namelist <- c('', 'Kenya')

for (name in namelist){
  print(paste0('Ongoing simulation for: ', name, ' data!'))

  # Loading true haplotype parameters for the simulation
  sim_Param <- readRDS(paste0(path, "dataset/true_Parameters", name, ".rds"))

  # Loading extra parameters
  extra_Sim_Param  <- readRDS(paste0(path, "dataset/extra_Parameters", name, ".rds"))

  # Simulation parameters
  n_Lbda        <- extra_Sim_Param[[1]]
  n_Sim_Loci    <- extra_Sim_Param[[2]]
  n_Hapl        <- extra_Sim_Param[[3]]
  n_Sampl       <- extra_Sim_Param[[4]]
  n_Sampl_Gen   <- extra_Sim_Param[[5]]
  n_Freq_Distr  <- extra_Sim_Param[[6]]

  n_Loci        <- log2(n_Hapl)
  true_Mean_MOI <- psi(sim_Param[[2]])

  # Reformatting true parameters to compute true prevalence
  reshap_Sim_Param <- vector(mode='list', length=n_Sim_Loci)

  for (i in 1:n_Sim_Loci){
    numb_Row <- n_Hapl[i]+1
    reshap_Sim_Param[[i]] <- vector(mode='list', length=n_Freq_Distr)

    for (j in 1:n_Freq_Distr){
      reshap_Sim_Param[[i]][[j]]             <- array(0, c(numb_Row, n_Lbda))
      reshap_Sim_Param[[i]][[j]][1,]         <- sim_Param[[2]]
      reshap_Sim_Param[[i]][[j]][2:numb_Row,] <- sim_Param[[1]][[i]][j,]
    }
  }

  main(sim_Param, reshap_Sim_Param, name)
}
