library(limma)

da_analysis <- function(data_matrix, condition, batch=NULL) {
    # data matrix is expected with features as rows and samples as columns, with rownames as feature names and colnames as sample names
    # run da analysis with limma, assuming exactly two conditions in 
    # class_label, optionally a batch_vector can be given for batch correction

    # Create design matrix for Limma
    conditions <- factor(condition)
    # if batch is given, add it as a fixed effect to the model
    if (!is.null(batch)) {
        batch <- factor(batch)
        design <- model.matrix(~conditions + batch)

    } else {
        design <- model.matrix(~conditions)
    }
    colnames(design) <- make.names(colnames(design))
    # Fit the linear model
    fit <- lmFit(data_matrix, design)
    efit <- eBayes(fit)
    # return the top table of results for the condition coefficient (the second column of the design matrix)
    results <- topTable(efit, coef=2, number=Inf)
  return(results)
}