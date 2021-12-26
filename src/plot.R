# Title        : Plotting the hapl. freq. and MOI Bias & CV, and also prevalence
# Objective    : Plot the bias and coefficient of variation of the estimates
# Created by   : christian Tsoungui Obama
# Created on   : 03.04.21
# Last modified: 25.12.21

# Importing libraries
library(dplyr)
library(ggplot2)

# Functions
beautify <- function (p, legende, pos, colpal, linety, name){
  p <- p + theme(panel.grid.minor = element_blank(),panel.grid=element_blank())
  p <- p + theme(panel.grid.major = element_blank(),panel.grid=element_blank())
  p <- p + theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
  p <- p + theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank())
  p <- p + theme(panel.background = element_rect(colour='black',fill="transparent"))
  p <- p + theme(plot.background = element_rect(fill = "transparent", colour = NA))

  # Title
  p <- p + theme(plot.title =element_text(face ="italic",size=rel(2.2),hjust = 0.5))

  # Legend
  p <- p + theme(legend.key=element_blank(), legend.text.align = 0, legend.background = element_blank())
  p <- p + theme(legend.key=element_blank(), legend.text.align = 0, legend.position = pos, legend.background = element_blank())
  p <- p + scale_linetype_manual(values=linety, labels=Ilab <-  legende, guide = guide_legend(title = "none"))
  p <- p + scale_colour_manual(values=colpal, labels=Ilab <-  legende, guide = guide_legend(title = name))
  p <- p + theme(legend.text = element_text(size = rel(2.2)))
  p <- p + theme(legend.title = element_text(size = rel(1.5),face="italic"), legend.margin = margin(t = 10, b = 0.1))
  p <- p + theme(legend.key.width = unit(9.5,"mm"))

  # Axis
  p <- p + theme(axis.text = element_text(colour='black'))
  p <- p + theme(axis.text = element_text(size = rel(2.1)))
  p <- p + theme(axis.title = element_text(size = rel(2.2)))
  p <- p + theme(axis.title.x = element_text(face="plain"))
  p <- p + theme(axis.title.y = element_text(face="plain"))
  p <- p + theme(axis.ticks = element_line(color = "black"))

  p
}

psi <- function (inp){
  inp / (1 - exp(-inp))
}

dataframe_builder_prev <- function(prev_estim, name, locNumb){
  NRow <- length(NSamp)*length(lbdavec)*NFreq*Hvec[locNumb] 

  cnames <- c('prev', 'sample', 'freq', 'shape')
  df1 <- array(0, dim = c(NRow, length(cnames)))
  df1 <- as.data.frame(df1)
  colnames(df1) <- cnames
  samp_vec <- rep(NSamp, each=NRow/(length(NSamp)*NFreq))
  df1[,'sample'] <- as.factor(rep(samp_vec, NFreq))
  df1[,'freq']   <- as.factor(rep(1:Hvec[locNumb], NRow/(length(lbdavec)*Hvec[locNumb])))
  df1[,'shape']  <- as.factor(rep(c("sym", "asym"), each=NRow/NFreq))

  exp_prev <- c()
  for (l in 1:NFreq){
    for (k in 1:length(NSamp)) {
      for (j in 1:NLbd){
        for (i in 1:Hvec[locNumb]){
          exp_prev <- c(exp_prev, prev_estim[[locNumb]][[k]][[j]][[l]][i])
        }
      }
    }
  }

  df1[,'prev'] <- exp_prev
  df1$type <- as.factor(name)
  df1
}

dataframe_builder_Freqperf <- function(perf_estim, locNumb){
  NRow <- length(NSamp)*length(lbdavec)*NFreq*Hvec[locNumb]

  # frequencies
  cnames <- c('bias', 'sample', 'freq', 'shape')
  df1 <- array(0, dim = c(NRow, length(cnames)))
  df1 <- as.data.frame(df1)
  colnames(df1) <- cnames

  samp_vec <- rep(NSamp, each=NRow/(length(NSamp)*NFreq))
  df1[,'sample'] <- as.factor(rep(samp_vec, NFreq))
  df1[,'freq']   <- as.factor(rep(1:Hvec[locNumb], NRow/(length(lbdavec)*Hvec[locNumb])))
  df1[,'shape']  <- as.factor(rep(c("sym", "asym"), each=NRow/NFreq))

  exp_perf <- c()
  for (l in 1:NFreq){
    for (k in 1:length(NSamp)) {
      for (j in 1:NLbd){
          exp_perf <- c(exp_perf, perf_estim[[locNumb]][[k]][[j]][[l]])
      }
    }
  }
  
  df1[,'bias'] <- exp_perf
  df1$type <- as.factor('freq')
  df1
}

dataframe_builder_Moiperf <- function(perf_estim, locNumb){
  NRow <- length(NSamp)*length(lbdavec)*NFreq

  # frequencies
  cnames <- c('bias', 'sample', 'shape')
  df1 <- array(0, dim = c(NRow, length(cnames)))
  df1 <- as.data.frame(df1)
  colnames(df1) <- cnames

  samp_vec <- rep(NSamp, each=NRow/(length(NSamp)*NFreq))
  df1[,'sample'] <- as.factor(rep(samp_vec, NFreq))
  df1[,'shape']  <- as.factor(rep(c("sym", "asym"), each=NRow/NFreq))

  exp_perf <- c()
  for (l in 1:NFreq){
    for (k in 1:length(NSamp)) {
      for (j in 1:NLbd){
          exp_perf <- c(exp_perf, perf_estim[[locNumb]][[k]][[j]][[l]])
      }
    }
  }
  
  df1[,'bias'] <- exp_perf
  df1$type <- as.factor('moi')
  df1
}

main <- function(ParTru, name){
  # Plots parameters
  shape_typ <- c('sym', 'asym')

  if(name=="Kenya"){
    dir        <- 'DD'
    shape_typ2 <- ''
  }else{
    dir        <- 'SD'
  }
  
  # Color palette (color-blind friendly) for the plots
  cbPalette <- c(rgb(0,0,0), rgb(.35, .70, .90), rgb(.90,.60,0), rgb(0,.60,.50), rgb(0,.45,.70), rgb(.80,.40,0), rgb(.5, .5, .5))
  lty <- c("dashed", "solid")

  if(1==1){ # Plotting prevalence
    # Importing the data to plot
    amb_prev   <- readRDS(paste0(path, "dataset/ambPrevalenceEstimates", name, ".rds"))
    unamb_prev <- readRDS(paste0(path, "dataset/unambPrevalenceEstimates", name, ".rds"))
    prev       <- readRDS(paste0(path, "dataset/prevalenceEstimates", name, ".rds"))

    # Plots parameters
    legende1 <- c('ambiguous', 'unambiguous', 'relative')

    for(l in 1:numbloci){ # 2 or 5 loci
      # Building the prevalence dataframe
      df_ambprev   <- dataframe_builder_prev(amb_prev, 'amb_prev', l)
      df_unambprev <- dataframe_builder_prev(unamb_prev, 'unamb_prev', l)
      df_prev      <- dataframe_builder_prev(prev, 'prev', l)

      df <- rbind(df_ambprev,df_unambprev,df_prev)
      tru_freq <- ParTru[[1]][[l]]

      for(k in 1:NFreq){  # sym or asym
        trufreq_vec <- tru_freq[k,]
        for(i in 1:Hvec[l]){
          if(trufreq_vec[i] != 0){
              for(j in NSamp){
                df1 <- df %>%
                      filter(sample == j, freq == i, shape == shape_typ[k]) %>%
                      droplevels()
                df1$lbd <- lbdavec
                p <- ggplot(data = df1, aes(x=lbd))
                p <- p + geom_line(aes(y = df1[,'prev'], color = type), size=1.)
                p <- beautify(p, legende1, c(0.20, 0.25), cbPalette, lty, 'Prevalence')
                p <- p + labs(x=expression(lambda), y="Prevalence", title=paste0("P = ", round(trufreq_vec[i], 3), ", N = ", j))
                p <- p + expand_limits(y=0)
                if(name == 'Kenya'){
                    outfile <- paste0(path,"plots/Prev_plots_", dir, "/prev_freq_", i, "_SSize_", j, "_year_", est_years[k], "_", name, ".pdf")
                }else{
                    outfile <- paste0(path,"plots/Prev_plots_", dir, "/prev_", shape_typ[k], "_freq_", i, "_SSize_", j, "_nloci_", NLoci[l], "_", name, ".pdf")
                }
                pdf(outfile, height=5, width=8)
                print(p)
                dev.off()
              }
          }
        }
      }
    }
  }


  legende1  <- NSamp

  if(1==1){ # Plotting bias for haplotype frequencies
    # Importing the data to plot
    freqbias <- readRDS(paste0(path, "dataset/freqbias", name, ".rds"))

      for(l in 1:numbloci){ # 2 or 5 loci
        # Building the frequencies bias dataframe
        df <- dataframe_builder_Freqperf(freqbias, l)
        tru_freq <- ParTru[[1]][[l]]

        for(k in 1:NFreq){  # sym or asym
          trufreq_vec <- tru_freq[k,]
          for(i in 1:Hvec[l]){
            if(trufreq_vec[i] != 0){
              for(j in NSamp){
                df1 <- df %>%
                    filter(freq == i, shape == shape_typ[k]) %>%
                    droplevels()
                df1$lbd <- lbdavec

                p <- ggplot(data = df1, aes(x=lbd))
                p <- p + geom_line(aes(y = df1[,'bias'], color = sample), size=1.)
                p <- beautify(p, legende1, NULL, cbPalette, lty, 'Sample')
                p <- p + expand_limits(y=0)
                p <- p + labs(x=expression(frac(lambda, 1 - e^-lambda)), y=paste0('Bias frequencies in %'), title = paste0('p = ', round(trufreq_vec[i], 3)))
                if(name == 'Kenya'){
                  outfile <- paste0(path,"plots/Bias_cv_plots_", dir, "/frequency/", "bias_freq_", i, "_year_", est_years[k], "_", name, ".pdf")
                }else{
                  outfile <- paste0(path,"plots/Bias_cv_plots_", dir, "/frequency/", "bias_", shape_typ[k], "_freq_", i, "_nloci_", NLoci[l], "_", name, ".pdf")
                }
                pdf(outfile, height=5, width=8)
                print(p)
                dev.off()
              }
            }
          }
        }
    }
  }

  if(1==1){ # Plotting bias and coefficient of variation for MOI
    # Importing the data to plot
    moibias  <- readRDS(paste0(path, "dataset/moibias", name, ".rds"))
    moicv    <- readRDS(paste0(path, "dataset/moicv", name, ".rds"))

      for(l in 1:numbloci){ # 2 or 5 loci
        # Building the prevalence dataframe
        df_cv   <- dataframe_builder_Moiperf(moicv, l)
        df_bias <- dataframe_builder_Moiperf(moibias, l)

        tru_freq <- ParTru[[1]][[l]]

        for(k in 1:NFreq){  # sym or asym
          trufreq_vec <- tru_freq[k,]
          df1 <- df_bias %>%
                filter(shape == shape_typ[k]) %>%
                droplevels()
          df1$lbd <- lbdavec

          p <- ggplot(data = df1, aes(x=lbd))
          p <- p + geom_line(aes(y = df1[,'bias'], color = sample), size=1.)
          p <- beautify(p, legende1, NULL, cbPalette, lty, 'Sample')
          p <- p + expand_limits(y=0)
          p <- p + labs(x=expression(frac(lambda, 1 - e^-lambda)), y=paste0('Bias MOI in %'))
          if(name == 'Kenya'){
              outfile <- paste0(path,"plots/Bias_cv_plots_", dir, "/moi/", "bias_moi_year_", est_years[k], "_", name, ".pdf")
          }else{
              outfile <- paste0(path,"plots/Bias_cv_plots_", dir, "/moi/", "bias_", shape_typ[k], "_moi_nloci_", NLoci[l], "_", name, ".pdf")
          }
          pdf(outfile, height=5, width=8)
          print(p)
          dev.off()

          df2 <- df_cv %>%
                filter(shape == shape_typ[k]) %>%
                droplevels()
          df2$lbd <- lbdavec

          p <- ggplot(data = df2, aes(x=lbd))
          p <- p + geom_line(aes(y = df2[,'bias'], color = sample), size=1.)
          p <- beautify(p, legende1, NULL, cbPalette, lty, 'Sample')
          p <- p + expand_limits(y=0)
          p <- p + labs(x=expression(frac(lambda, 1 - e^-lambda)), y=paste0('Coef. variation MOI'))
           if(name == 'Kenya'){
              outfile <- paste0(path,"plots/Bias_cv_plots_", dir, "/moi/", "coefvar_moi_year_", est_years[k], "_", name, ".pdf")
          }else{
              outfile <- paste0(path,"plots/Bias_cv_plots_", dir, "/moi/", "coefvar_", shape_typ[k], "_moi_nloci_", NLoci[l], "_", name, ".pdf")
          }
          pdf(outfile, height=5, width=8)
          print(p)
          dev.off()
        }
      }
  }
}

path <- "/Volumes/GoogleDrive-117934057836063832284/My Drive/Maths against Malaria/Christian/Models/MultiLociBiallelicModel/"

# Define data origin
name <- ''

# Importing extra parameters
parExtr  <- readRDS(paste0(path, "dataset/extraParameters", name, ".rds"))
parTrue  <- readRDS(paste0(path, "dataset/trueParameters", name, ".rds"))

# Running the plots ('' <- simualted data, kenya <- kenyan data)
NLbd      <- parExtr[[1]]
Nn        <- parExtr[[2]]
Hvec      <- parExtr[[3]]
NN        <- parExtr[[4]]
NEst      <- parExtr[[5]]
NFreq     <- parExtr[[6]]
NLoci     <- log2(Hvec)
lbdavec   <- parTrue[[2]]
NSamp     <- parTrue[[3]]
lbda      <- psi(lbdavec)
numbloci  <- length(Hvec)
est_years <- c(2005, 2010)

# Plotting
main(parTrue, name)