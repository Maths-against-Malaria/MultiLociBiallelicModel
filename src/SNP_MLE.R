# Title        : Template script for the MLE of haplotype frequencies, MOI, and prevalence
#                from example dataset.
# Objective    : Estimate haplotype frequencies, MOI, and prevalence
#                from genomic/molecular data
# Created by   : Christian Tsoungui Obama
# Created on   : 25.04.22
# Last modified: 26.04.22

# Install the necessary packages if necessary
i#nstall.packages('xlsx')   # Comment this line if xlsx installed

# Loading libraries
library(xlsx)

# Importing the reformatted data as '.xlsx' file
# path <- "/home/janedoe/Documents/"
# DATA <- read.xlsx('/home/janedoe/Documents/example.xlsx', 1, header = TRUE)
path1 <- "/Volumes/GoogleDrive-117934057836063832284/My Drive/Maths against Malaria/Christian/dataset/"
path2 <- "/Volumes/GoogleDrive-117934057836063832284/My Drive/Maths against Malaria/Christian/Models/MultiLociBiallelicModel/"

DATA1 <- read.xlsx(paste0(path1,'CameroonMcCollum2008-SNP.xlsx'), 5, header = TRUE)
pick1 <- rowSums(is.na(DATA1))<1
DATA1 <- DATA1[pick1,]

DATA2 <- read.xlsx(paste0(path1,'CameroonMcCollum2008-SNP.xlsx'), 6, header = TRUE)
pick2 <- rowSums(is.na(DATA2))<1
DATA2 <- DATA2[pick2,]

# Loading external resources
# source("/home/janedoe/Documents/SNPModel.R")
source(paste0(path2,"src/SNPModel.R"))

# Finding the MLEs
est1 <- mle(DATA1, id=TRUE)
est2 <- mle(DATA2, id=TRUE)

# Estimating prevalence
## Unobservable prevalence
unobsprev1 <- estunobsprev(est1)
unobsprev2 <- estunobsprev(est2)

## Conditional prevalence
condprev1 <- estcondprev(est1)
condprev2 <- estcondprev(est2)

## Relative prevalence
relprev1 <- estrelprev(DATA1, id=TRUE)
relprev2 <- estrelprev(DATA2, id=TRUE)
