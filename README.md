# MultiLociBiallelicModel
This repository contains the implementation of the model described in the paper "A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNPs data".

The script "SNPModel.R" (in folder "src") contains the implementation of the method to obtain the maximum-likelihood estimates (MLE) of haplotype frequencies and multiplicity of infection (MOI). Moreover, the script contains the method to estimates to estimate haplotypes prevalence using the MLEs of haplotype frequencies and MOI as plug-in estimates. The script contains additional functions necessary for the functionning of the method (see "User manual"). 

The script "SNP_MLE.R" is a template that exemplifies the use of the method to obtain the MLEs, and the estimates of prevalence.
