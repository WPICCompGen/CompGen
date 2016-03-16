current_script_dir <- dirname(sys.frame(1)$ofile)
setwd(current_script_dir)
source("source_DAWN.R")  ##source file

## Load gene expression data
expdata <- read.table("sample_input_expression.txt",header=T)

## Load GLA gene names and p-values
GLA <- read.table('sample_input_GLA.txt',header=T)  ##load the GLA p-values
GLA_genenames <- GLA[, 'GLA_genenames']
GLA_pvalue <- GLA[, 'GLA_pvalue']

## Load genes to have fixed hidden states
FG <- read.delim("sample_input_fixed_genes.txt", sep='\t', header=F)
fixedGene <- as.character(FG[,1])

## Load genes in the additional covariate of HMRF analysis
CG <- read.delim("sample_input_additional_covariates.txt", sep='\t', header=F)
covGene <- as.character(CG[,1])

## Adjust parameters
pthres_pns <- 0.1  #threshold for p-value in PNS algorithm, default = 0.1
corthres <- 0.7  #threshold for pairwise-correlation in PNS algorithm, default = 0.7
lambda <- 0.24  #tuning parameter for lasso, default = 0.24; see user manual for more details
pthres_hmrf <- 0.05  #threshold for p-value in HMRF analysis, default = 0.05
iter <- 100  #number of iterations in HMRF analysis, default= 100
trimthres <- 5  #threshold for Z score trimming, default = 5; see user manual for more details
file_out_pns <- "pns.network"  #PNS network output file directory, default = 'pns.network'
file_out_hmrf <- "DAWN_result.csv"  #HMRF analysis output file directory, default = 'DAWN_result.csv'

##Function for performing a complete DAWN analysis
identify_risk_genes(GLA_pvalue = GLA_pvalue, GLA_genenames = GLA_genenames, expdata = expdata)
