# Title        : MLE of lineage frequencies, MOI, and prevalence
# Objective    : Estimate haplotype frequencies, MOI, and prevalence
#                from genomic/molecular data
# Created by   : Christian Tsoungui Obama
# Created on   : 25.04.22
# Last modified: 26.04.22

# Install the necessary packages if necessary
#install.packages('xlsx')   # Comment this line if xlsx installed

# Loading libraries
library(xlsx)

# Importing the reformatted data as '.xlsx' file
#path <- "/home/john_doe/Documents/"
path <- "/Volumes/GoogleDrive-117934057836063832284/My Drive/Maths against Malaria/Christian/Models/MultiLociBiallelicModel/"
df <- read.xlsx(paste0(path,'dataset/exampleData.xlsx'),1)

# Loading external ressources
source(paste0(path, "src/SNPModel.R")) 

# Finding the MLEs
est <- mle(df, id=TRUE)

# Estimating prevalence
## Unobservable prevalence
unobsprev <- estunobsprev(est)

## Conditional prevalence
condprev <- estcondprev(est)

## Relative prevalence
relprev <- estrelprev(df, id=TRUE)