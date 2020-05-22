

make_dge_and_contrasts = function(count_data_frame, contrast_column) {
  
  
  #### Imports ####
  if (!requireNamespace("pkg", quietly = TRUE)) {
    stop("Package \"pkg\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  treatments = unique(str_replace(colnames(count_data_frame), "_[0-9]$", ""))
  group <- rep(1:length(treatments), each = 3)
  dge_object <-
    DGEList(
      data.matrix(count_data_frame),
      genes = rownames(count_data_frame),
      group = group
    )
  
  # Perform TMM normalisation
  dge_object <- calcNormFactors(dge_object, method = "TMM")
  
  # defining the design matrix
  factor <- factor(group)
  treatments_factors <-
    unlist(lapply(treatments, function(x)
      rep(x, 3)))
  design <- model.matrix( ~ 0 + group, data = dgeObj$samples)
  
  # all against the defined column
  combs <- sapply(treatments[1:length(treatments)], function(x) {
    paste0(x, paste0("-", contrast_column))
  }, USE.NAMES = FALSE)
  combs <- combs[-c(2)]
  
  #contr = makeContrasts(contrasts = combs, levels = factor(treatments))
  contr <- makeContrasts(contrasts = combs, levels = factor(treatments))
  
  dge_object <-
    estimateCommonDisp(dge_object, design = design, robust = TRUE)
  dge_object <- estimateGLMTrendedDisp(dge_object, design = design)
  dge_object <- estimateTagwiseDisp(dge_object, design = design)
  
  
  plotBCV(dgeObj)
  plotMDS(dgeObj)
  
  return(list(dge_object = dge_object, contrasts = contr))
}




volcano_plot_pipeline = function(dge_object, contrasts){
  
  results_fc_volc = data.frame(logFC = c(), neg_log_pval = c())
  
  for (i in 1:dim(contrasts)[2]) {
    # fit using generalized linear model
    fit <- glmFit(dge_object)
    lrt.BvsL <- glmLRT(fit, contrast = contrasts[, i])
    topResults = topTags(lrt.BvsL, n = dim(dge_object$counts)[1])$table
    
    # add results to 
    results_fc_volc = rbind(results_fc_volc,
                            data.frame(
                              logFC = topResults$logFC ,
                              neg_log_pval = -log10(topResults$FDR)
                            ))
  }
  # volcano plot of negative log p values vs log FC
  results_fc_volc = results_fc_volc[results_fc_volc$neg_log_pval != Inf , ]
  
  # plot volcanoplot 
  plot(
    results_fc_volc,
    pch = 20,
    col = sapply(results_fc_volc$logFC, function(x)
      ifelse(abs(x) > 2, "red", "blue")),
    main = "Volcano plot of differential expression, blue: abs logFC > 2 "
  )
}


run_diff_exp <- function(dge_object, design_matrix, contrasts) {
  # Compute an object containing
  
  results = data.frame(matrix(nrow = nrow(dge_object), ncol = 0))
  rownames(results) = toupper(rownames(dge_object))
  
  results_pValues = data.frame(matrix(nrow = nrow(dge_object), ncol = 0))
  rownames(results_pValues) = toupper(rownames(dge_object))
  results_directions = data.frame(matrix(nrow = nrow(dge_object), ncol =
                                           0))
  rownames(results_directions) = toupper(rownames(dge_object))
  for (i in 1:dim(contrasts)[2]) {
    # fit using generalized linear model
    fit <- glmFit(dge_object, design_matrix)
    lrt.BvsL <- glmLRT(fit, contrast = contrasts[, i])
    
    # results for overrepresentation analysis
    results[colnames(contrasts)[i]] = decideTests(lrt.BvsL, p.value = 0.05)
    
    # results for Gene set enrichment analysis
    # results need to be ordered as defined in the results df above
    topResults = topTags(lrt.BvsL, n = dim(dge_object$counts)[1])$table
    
    results_pValues[colnames(contrasts)[i]] = topResults[order(rownames(topResults)), ]$FDR
    results_directions[colnames(contrasts)[i]] = topResults[order(rownames(topResults)), ]$logFC
  }
  return(list("pValues" = results_pValues, "directions" = results_directions))
}