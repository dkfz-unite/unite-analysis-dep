library(limma)

da_analysis <- function(data_matrix, condition, ref_category = NULL, batch=NULL) {
    # data matrix is expected with features as rows and samples as columns, with rownames as feature names and colnames as sample names
    # run da analysis with limma, assuming exactly two conditions in 
    # class_label, optionally a batch_vector can be given for batch correction

    # convert to factor
    conditions <- factor(condition)
    # check expected two levels
    if (nlevels(conditions) != 2) {
        stop("da_analysis requires exactly two conditions, got: ", nlevels(conditions))
    }
    # if reference category is set explicitly set it here
    if (!is.null(ref_category)) {
        conditions <- relevel(conditions, ref=as.character(ref_category))
    }
    # if batch is given, add it as a fixed effect to the model
    if (!is.null(batch)) {
        batch <- factor(batch)
        design <- model.matrix(~conditions + batch)

    } else {
        design <- model.matrix(~conditions)
    }
    colnames(design) <- make.names(sub("^conditions", "condition", colnames(design)))
    # Fit the linear model
    fit <- lmFit(data_matrix, design)
    efit <- eBayes(fit)
    # return the top table of results for the condition coefficient (the second column of the design matrix)
    results <- topTable(efit, coef=2, number=Inf)
    # write the name of the contrast of the end of the output table
    ref = levels(conditions)[1]
    test = levels(conditions)[2]

    results$contrast <- paste0(test, " - ", ref)
    return(results)
}