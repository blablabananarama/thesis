
cli = length(args) == 0

####################################################################################
# pre-processing (should in future be done in python and should only translate to matrix)
####################################################################################
if (!cli){
  
  count_matrix_df = read.csv("~/Documents/masters/jupyter_notebooks/data/proteomics/protein_values.csv")
  
  start_column = grep("\\Glucose\\b", colnames(count_matrix_df))
  end_column = grep("\\Fructose\\b", colnames(count_matrix_df))
  
  start_column_unc = grep("\\Glucose.1\\b", colnames(count_matrix_df))
  end_column_unc = grep("\\Fructose.1\\b", colnames(count_matrix_df))
  
  
  
  # 
  seqdata = count_matrix_df[,c(1,start_column[1]:end_column[1])]
  unc_matrix = data.matrix(count_matrix_df[,start_column_unc[1]:end_column_unc[1]])
  
  # remove duplicates (to be done in python)
  seqdata = seqdata[!duplicated(seqdata[,1]),]
  
}


####################################################################################
# create synthetic dataset by bootstrap sampling
####################################################################################
library("gdata")

# 
# pseudo_count_matrix = matrix(nrow = nrow(count_matrix)*3, ncol = ncol(count_matrix))
# 
# for(k in 0:2){
#   factor = dim(count_matrix)[1] * k
#     for(i in 1:dim(count_matrix)[1]){
#       for(j in 1:dim(count_matrix)[2]){
#         pseudo_count_matrix[i+factor, j] = round(rnorm(1, mean=count_matrix[i,j],sd=unc_matrix[i,j]))
#     }
#   }
# }


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
countdata <- seqdata[,-1]

# Store ProteinID as rownames
rownames(countdata) <- seqdata[,1]
countdata[rownames(countdata) %in% rownames(countdata)[2], ]
countdata

intensities = log2(data.matrix(countdata))
intensities[intensities == -Inf] = NA

combs = c(paste(combn(colnames(countdata), 2)[1,], "-",  combn(colnames(countdata), 2)[2,], sep = ""))


####################################################################################
# create design matrix from data matrix 
####################################################################################
treatments = colnames(countdata)

design <- diag(length(treatments))
colnames(design) = treatments

design_matrix = model.matrix(~ 0 + treatments)

contr = makeContrasts(contrasts = combs, levels = treatments)

####################################################################################
# testing for differential expression
####################################################################################
# fit a linear model 
fit = lmFit(intensities, design_matrix)

# 
fit2 = contrasts.fit(fit, contr) 
fit2 = eBayes(fit2)

res = decideTests(fit2, lfc=1, p.value = 0.05 )

write.csv(results$table,paste(results_path, "/limma_de_analysis.csv", sep = ""))

