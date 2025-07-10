# MultiLociBiallelicModel
This repository contains the implementation of the model described in the publication:

Tsoungui Obama HCJ and Schneider KA (2022) A maximum-likelihood method to estimate haplotype frequencies and prevalence alongside multiplicity of infection from SNP data. Front. Epidemiol. 2:943625. doi: 10.3389/fepid.2022.943625

The script "SNPModel.R" (in folder "src") contains the implementation of the method to obtain the maximum-likelihood estimates (MLE) of haplotype frequencies and multiplicity of infection (MOI). Moreover, the script contains the method to estimates to estimate haplotypes prevalence using the MLEs of haplotype frequencies and MOI as plug-in estimates. The script contains additional functions necessary for the functionning of the method (see "User manual"). 

The script "SNP_MLE.R" is a template that exemplifies the use of the method to obtain the MLEs, and the estimates of prevalence.
