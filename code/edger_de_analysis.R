#!/usr/bin/env Rscript

wants <- c('tidyverse')
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has], repos="http://ftp.acc.umu.se/mirror/CRAN/")

########BIOCONDUCTOR DEPENDENCIES######################3

########BIOCONDUCTOR DEPENDENCIES######################3


#load all the packages
lapply(wants, require, character.only = TRUE)


#find the "code" directory, this should be one folder away from the project root
getCurrentFileLocation <-  function()
{
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}
folder_vec <- unlist(strsplit(getCurrentFileLocation(), "/"))
code_folder_location <- tail(which(folder_vec == "code"), n=1)
code_folder_location <- paste(folder_vec[1:code_folder_location], collapse = '/')
setwd(code_folder_location)

#go one folder up and set as project root directory
setwd("../")
PROJ <- getwd()
setwd(paste0(PROJ, '/data/'))

#set folder paths relative to the project root
CURRENT <- getwd()
DATA <- getwd()
RAW_EXTERNAL <- paste0(DATA, '/raw_external/')
RAW_INTERNAL <- paste0(DATA, '/raw_internal/')
INTERMEDIATE = paste0(DATA, '/intermediate/')
FINAL = paste0(DATA, '/final/')

RESULTS <- paste0(PROJ, '/results/')
FIGURES <- paste0(RESULTS, 'figures/')
PICTURES <- paste0(RESULTS, 'pictures/')


####################################################################################
# Differential expression study starting from  count matrices 
####################################################################################

cli = length(args) == 0

####################################################################################
# pre-processing (should in future be done in python and should only translate to matrix)
####################################################################################

# read in raw data from file
count_matrix_df = read.csv(paste0(RAW_INTERNAL, "proteomics/protein_values.csv"))

# remove repeated value
count_matrix_df = count_matrix_df[!duplicated(count_matrix_df[,"Uniprot.Accession"]),]


start_column = grep("\\Glucose\\b", colnames(count_matrix_df))
end_column = grep("\\Fructose\\b", colnames(count_matrix_df))

start_column_unc = grep("\\Glucose.1\\b", colnames(count_matrix_df))
end_column_unc = grep("\\Fructose.1\\b", colnames(count_matrix_df))



# extract uncertainty values from data
unc_matrix = data.matrix(count_matrix_df[,start_column_unc[1]:end_column_unc[1]])


####################################################################################
# exampledata for trial of the whatever shit  
####################################################################################
library(edgeR)
library(limma)
library(Glimma)
library(ggplot2)
#library(org.Mm.eg.db)

# read in file provided on the command line
args = commandArgs(trailingOnly=TRUE)
print(args)
data = args[1]

results_path = args[2]


if(cli == TRUE){seqdata <- read.csv(data)}


## Calculate the Counts Per Million measure
# myCPM <- cpm(countdata)


# Subset the rows of countdata without nan values
counts.keep <- count_matrix_df[,start_column[1]:end_column[1]]
# set rownames to uniprot id
rownames(counts.keep) <- count_matrix_df[, "Uniprot.Accession"]
unc_matrix <- unc_matrix[complete.cases(counts.keep),]
counts.keep <- counts.keep[complete.cases(counts.keep),]
col_group <- rep(1:3, each=22)

# fake counts 
fake_counts = matrix(nrow=nrow(counts.keep), ncol = ncol(counts.keep)*2)
for(i in 0:1){
  for(j in 1:nrow(fake_counts)){
    for(k in 1:ncol(counts.keep)){
      fake_counts[j,k+(i*ncol(counts.keep))] = round(rnorm(1, counts.keep[j,k],unc_matrix[j,k]))
    }
  }
}

counts.keep <- cbind(counts.keep, as.data.frame(fake_counts))
colnames(counts.keep) <- sapply(c(1:dim(counts.keep)[2]), function(x) paste(colnames(counts.keep)[x], col_group[x], sep = ""))

## Convert to an edgeR object
group <- rep(1:22, 3)
dgeObj <- DGEList(counts.keep, group = group)

## Perform TMM normalisation
dgeObj <- calcNormFactors(dgeObj)




####################################################################################
# create design matrix from data matrix 
####################################################################################
# assuming rownames as treatments
treatments = colnames(count_matrix_df[,start_column[1]:end_column[1]])

design <- rbind(model.matrix(~ -1 + treatments), model.matrix(~ -1 + treatments), model.matrix(~ -1 + treatments))

# all against all combinations
# combs = c(paste(combn(colnames(countdata), 2)[1,], "-",  combn(colnames(countdata), 2)[2,], sep = ""))

# all against glucose combinations 
combs = sapply(treatments[2:length(treatments)], function(x){paste(x, "-Glucose", sep = "")}, USE.NAMES = FALSE)

# stress analysis against glucose
# combs = sapply(treatments[15:17], function(x){paste(x, "-Glucose", sep = "")}, USE.NAMES = FALSE)

# all stationary phase 


# glucose vs complex (LB)


# all chemostat


contr = makeContrasts(contrasts = combs, levels = factor(treatments))

####################################################################################
# estimate dispersion of the data  
####################################################################################

dgeObj <- estimateCommonDisp(dgeObj)
dgeObj <- estimateGLMTrendedDisp(dgeObj)
dgeObj <- estimateTagwiseDisp(dgeObj)

if(!cli){
  plotBCV(dgeObj)
  plotMDS(dgeObj)
}
####################################################################################
# testing for differential expression
####################################################################################
# fit a linear model (generalized)
names(fit)
head(coef(fit))

# create test results for all contrasts 
# results = as.data.frame(matrix(0, ncol = 0, nrow = dim(dgeObj$counts)[1]))
results = data.frame(matrix(nrow=nrow(counts.keep), ncol=0))
rownames(results) = rownames(counts.keep)
for(i in 1:dim(contr)[2]){
  fit <- glmFit(dgeObj, design)
  lrt.BvsL <- glmLRT(fit, contrast = contr[,i])
  # order the PValues by the names of the proteins, so that the rownames show the same name as the pvalue
  # results[colnames(contr)[i]] = topTags(lrt.BvsL, n = dim(dgeObj$counts)[1])$table$PValue[order(rownames(topTags(lrt.BvsL, n = dim(dgeObj$counts)[1])$table))]
  results[colnames(contr)[i]] = decideTests(lrt.BvsL)
  #rownames(results) = rownames(topTags(lrt.BvsL, n = dim(dgeObj$counts)[1])$table)[order(rownames(topTags(lrt.BvsL, n = dim(dgeObj$counts)[1])$table))]
  result = topTags(lrt.BvsL, n = dim(dgeObj$counts)[1])$table
}


# save test output
write.csv(results,paste0(FINAL, "de_edger_results_stress.csv"))
