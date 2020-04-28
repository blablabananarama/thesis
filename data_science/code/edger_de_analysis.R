#!/usr/bin/env Rscript

wants <- c('tidyverse', 'edgeR', 'limma', 'Glimma', 'ggplot2', 'piano', 'pheatmap', 'RColorBrewer', 'viridis', 'reshape2')
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
count_matrix_df <- count_matrix_df[!duplicated(count_matrix_df$Gene),]
count_matrix_df <- count_matrix_df[!is.na(count_matrix_df$Gene),]
count_matrix_df$Gene <- toupper(count_matrix_df$Gene)

start_column = grep("\\Glucose\\b", colnames(count_matrix_df))
end_column = grep("\\Fructose\\b", colnames(count_matrix_df))

start_column_unc = grep("\\Glucose.1\\b", colnames(count_matrix_df))
end_column_unc = grep("\\Fructose.1\\b", colnames(count_matrix_df))


# extract uncertainty values from data
unc_matrix = data.matrix(count_matrix_df[,start_column_unc[1]:end_column_unc[1]])



####################################################################################
# Data computed from mean and std values  
####################################################################################

# read in file provided on the command line
args = commandArgs(trailingOnly=TRUE)
print(args)
data = args[1]

results_path = args[2]


## Calculate the Counts Per Million measure
# myCPM <- cpm(countdata)


# Subset the rows of countdata without nan values
counts.keep <- count_matrix_df[,start_column[1]:end_column[1]]
# set rownames to uniprot id
rownames(counts.keep) <- count_matrix_df$Gene
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

# p-values and fch as directions
results_pValues = data.frame(matrix(nrow=nrow(counts.keep), ncol=0))
rownames(results_pValues) = rownames(counts.keep)
results_directions = data.frame(matrix(nrow=nrow(counts.keep), ncol=0))
rownames(results_directions) = rownames(counts.keep)

for(i in 1:dim(contr)[2]){
  # fit using generalized linear model
  fit <- glmFit(dgeObj, design)
  lrt.BvsL <- glmLRT(fit, contrast = contr[,i])
  
  # results for overrepresentation analysis 
  results[colnames(contr)[i]] = decideTests(lrt.BvsL, p.value=0.01)
  
  # results for Gene set enrichment analysis
  # results need to be ordered as defined in the results df above
  results_pValues[colnames(contr)[i]] = topTags(lrt.BvsL, n = dim(dgeObj$counts)[1])$table[order(rownames(topTags(lrt.BvsL, n = dim(dgeObj$counts)[1])$table)),]$FDR
  results_directions[colnames(contr)[i]] = topTags(lrt.BvsL, n = dim(dgeObj$counts)[1])$table[order(rownames(topTags(lrt.BvsL, n = dim(dgeObj$counts)[1])$table)),]$logFC
  
  }

# differentially expressed genes
up_down_regulated_table <- as.data.frame(apply(results, MARGIN = 2, table))
up_down_regulated_table["direction"] <- c("downregulated", "not significant", "upregulated")
up_down_regulated_table_molten <- melt(up_down_regulated_table, id.vars="direction")

ggplot(data=up_down_regulated_table_molten, aes(x=variable, y=value, fill=direction)) +
  geom_bar(stat="identity") +
  ggtitle("N genes Up- and Downregulated at alpha = 0.01") +
  theme(axis.text.x = element_text(angle = 90))



####################################################################################
# piano GO analysis
####################################################################################
# load gene set, http://www.go2msig.org/cgi-bin/prebuilt.cgi?taxid=511145 , gene identifier: entrez
# goSets <- loadGSC(file = paste0(RAW_EXTERNAL, "go_terms/Escherichia_coli_str_K12_substr_MG1655_GSEA_GO_sets_all_ids_April_2015.gmt"))

# load gene set, http://ge-lab.org/gskb/ , gene identifier: gene names
goSets <- loadGSC(file = paste0(RAW_EXTERNAL, "go_terms/Escherichia_coli_k12_intestinal-bacteria_gmt.gmt"))


gsa_results <- list()
# run gsa
for(i in colnames(results_pValues)){
  gsa_results[[i]] <- GSAsummaryTable(runGSA(results_pValues[i], directions = results_directions[i] ,geneSetStat = "fisher",  gsc=goSets, verbose=TRUE))
  
}



# present results
# upregulated 
gsa_results_up = matrix(unlist(lapply(gsa_results, function(x) return(x["p (mix.dir.up)"]))), ncol = length(gsa_results))

pheatmap(mat = gsa_results_up, cluster_rows = F, color = inferno(10), labels_col = names(gsa_results))

# downregulated 
gsa_results_down = matrix(unlist(lapply(gsa_results, function(x) return(x["p (mix.dir.dn)"]))), ncol = length(gsa_results))

pheatmap(mat = gsa_results_down, cluster_rows = F, color = inferno(10), labels_col = names(gsa_results))

# both directions
gsa_results_both = matrix(unlist(lapply(gsa_results, function(x) return(x["p (non-dir.)"]))), ncol = length(gsa_results))

pheatmap(mat = gsa_results_both, cluster_rows = F, color = inferno(10), labels_col = names(gsa_results))



# BPs
# find BPs
bps <- unlist(lapply(strsplit(gsa_results$`LB-Glucose`$Name, "_"), function(x) return (if(x[2]=="BP") T else F)))
gsa_bps <- lapply(gsa_results, function(x) x[bps,])

# plot top BPs
as.data.frame(lapply(gsa_bps ,function (x) x$Name[order(x$`p (non-dir.)`)][1:10]))

# MFs
# find MFs
mfs <- unlist(lapply(strsplit(gsa_results$`LB-Glucose`$Name, "_"), function(x) return (if(x[2]=="MF") T else F)))
gsa_mfs <- lapply(gsa_results, function(x) x[mfs,])

# plot top BPs
as.data.frame(lapply(gsa_mfs ,function (x) x$Name[order(x$`p (non-dir.)`)][1:10]))









####################################################################################
# save plots and files
####################################################################################


# Differential expression plots
# differentially expressed genes
up_down_regulated_table <- as.data.frame(apply(results, MARGIN = 2, table))
up_down_regulated_table["direction"] <- c("downregulated", "not significant", "upregulated")
up_down_regulated_table_molten <- melt(up_down_regulated_table, id.vars="direction")

####################################################################################0.05
# results = as.data.frame(matrix(0, ncol = 0, nrow = dim(dgeObj$counts)[1]))
results = data.frame(matrix(nrow=nrow(counts.keep), ncol=0))
rownames(results) = rownames(counts.keep)

for(i in 1:dim(contr)[2]){
  # fit using generalized linear model
  fit <- glmFit(dgeObj, design)
  lrt.BvsL <- glmLRT(fit, contrast = contr[,i])
  
  # results for overrepresentation analysis 
  results[colnames(contr)[i]] = decideTests(lrt.BvsL, p.value=0.05)
  
}
# differentially expressed genes
up_down_regulated_table <- as.data.frame(apply(results, MARGIN = 2, table))
up_down_regulated_table["direction"] <- c("downregulated", "not significant", "upregulated")
up_down_regulated_table_molten <- melt(up_down_regulated_table, id.vars="direction")

# cutoff 0.01
png(file=paste0(FIGURES,"go_analysis/differential-expression_basic-proteomics-dataset_n-up-downregulated05.png"), width=900, height=600)
ggplot(data=up_down_regulated_table_molten, aes(x=variable, y=value, fill=direction)) +
  geom_bar(stat="identity") +
  ggtitle("N genes Up- and Downregulated at alpha = 0.05") +
  theme(axis.text.x = element_text(angle = 90))
dev.off()
####################
################## 0.01
# results = as.data.frame(matrix(0, ncol = 0, nrow = dim(dgeObj$counts)[1]))
results = data.frame(matrix(nrow=nrow(counts.keep), ncol=0))
rownames(results) = rownames(counts.keep)

for(i in 1:dim(contr)[2]){
  # fit using generalized linear model
  fit <- glmFit(dgeObj, design)
  lrt.BvsL <- glmLRT(fit, contrast = contr[,i])
  
  # results for overrepresentation analysis 
  results[colnames(contr)[i]] = decideTests(lrt.BvsL, p.value=0.01)
  
}
# differentially expressed genes
up_down_regulated_table <- as.data.frame(apply(results, MARGIN = 2, table))
up_down_regulated_table["direction"] <- c("downregulated", "not significant", "upregulated")
up_down_regulated_table_molten <- melt(up_down_regulated_table, id.vars="direction")

# cutoff 0.01
png(file=paste0(FIGURES,"go_analysis/differential-expression_basic-proteomics-dataset_n-up-downregulated01.png"), width=900, height=600)
ggplot(data=up_down_regulated_table_molten, aes(x=variable, y=value, fill=direction)) +
  geom_bar(stat="identity") +
  ggtitle("N genes Up- and Downregulated at alpha = 0.01") +
  theme(axis.text.x = element_text(angle = 90))
dev.off()
####################
################## 0.001
# results = as.data.frame(matrix(0, ncol = 0, nrow = dim(dgeObj$counts)[1]))
results = data.frame(matrix(nrow=nrow(counts.keep), ncol=0))
rownames(results) = rownames(counts.keep)

for(i in 1:dim(contr)[2]){
  # fit using generalized linear model
  fit <- glmFit(dgeObj, design)
  lrt.BvsL <- glmLRT(fit, contrast = contr[,i])
  
  # results for overrepresentation analysis 
  results[colnames(contr)[i]] = decideTests(lrt.BvsL, p.value=0.001)
  
}
# differentially expressed genes
up_down_regulated_table <- as.data.frame(apply(results, MARGIN = 2, table))
up_down_regulated_table["direction"] <- c("downregulated", "not significant", "upregulated")
up_down_regulated_table_molten <- melt(up_down_regulated_table, id.vars="direction")

# cutoff 0.01
png(file=paste0(FIGURES,"go_analysis/differential-expression_basic-proteomics-dataset_n-up-downregulated001.png"), width=900, height=600)
ggplot(data=up_down_regulated_table_molten, aes(x=variable, y=value, fill=direction)) +
  geom_bar(stat="identity") +
  ggtitle("N genes Up- and Downregulated at alpha = 0.001") +
  theme(axis.text.x = element_text(angle = 90))
dev.off()
####################################################################################


# Go terms plots
# upregulated
png(file=paste0(FIGURES,"go_analysis/go-p-values-fisher-combined_basic-proteomics-dataset_upregulated.png"), width=450, height=900)
pheatmap(mat = gsa_results_up, cluster_rows = F, color = inferno(10), labels_col = names(gsa_results), main = "P-values of GO terms for upregulated genes only")
dev.off()

# Downregulated
png(file=paste0(FIGURES,"go_analysis/go-p-values-fisher-combined_basic-proteomics-dataset_downregulated.png"), width=450, height=900)
pheatmap(mat = gsa_results_down, cluster_rows = F, color = inferno(10), labels_col = names(gsa_results), main = "P-values of GO terms for downregulated genes only")
dev.off()

# both directions
png(file=paste0(FIGURES,"go_analysis/go-p-values-fisher-combined_basic-proteomics-dataset_up-and-downregulated.png"), width=450, height=900)
pheatmap(mat = gsa_results_both, cluster_rows = F, color = inferno(10), labels_col = names(gsa_results), main = "P-values of GO terms for both up- and downregulated genes")
dev.off()





