#' Rank CpG probes
#'
#' Apply a linear model to Illumina's 450k or EPIC methylation data to get the
#' t-value of each CpG probe
#'
#' @param methData A data frame or a matrix containing Illumina's CpG probes in
#' rows and samples in columns. A SummarizedExperiment object can be used too
#' @param pheno A data frame or a matrix containing samples in rows and
#' covariates in columns. If NULL (default), pheno is extracted from the
#' SummarizedExperiment object
#' @param explanatory The column name or position from pheno used to perform the
#'  comparison between groups (default = first column)
#' @param covariates A list or character vector with column names from pheno
#' used as data covariates in the linear model
#' @param refGroup The group name or position from explanatory variable used to
#' perform the comparison (default = first group)
#' @param continuous A list or character vector with columns names from pheno
#' which should be treated as continuous variables (default = none)
#' @param typeInput Type of input methylation data. "beta" for Beta-values and
#' "M" for M-values
#' @param typeAnalysis "M" to use M-values to rank the CpG probes (default).
#' "beta" to use Beta-values instead
#'
#' @return A named vector containing the t-values from the linear model for each
#'  CpG probe
#'
#' @author Jordi Martorell Marugán, \email{jordi.martorell@@genyo.es}
#'
#' @references Smyth, G. K. (2005). \emph{Limma: linear models for microarray
#' data}. Bioinformatics and Computational Biology Solutions using R and
#' Bioconductor, 397-420.
#'
#' @seealso \code{\link{mCSEATest}}
#'
#' @examples
#' data(mcseadata)
#' myRank <- rankProbes(betaTest, phenoTest, refGroup = "Control")
#' head(myRank)
#' @export

rankProbes <- function(methData, pheno = NULL, explanatory = 1,
                    covariates = c(), refGroup = 1, continuous = NULL,
                    typeInput = "beta", typeAnalysis = "M")
    {

    # Check input objects
    if (!any(class(methData) == "data.frame" | class(methData) == "matrix" |
            class(methData) == "SummarizedExperiment" |
            class(methData) == "RangedSummarizedExperiment")){
        stop("methData must be a data frame, a matrix or a SummarizedExperiment
            object")
    }

    if (!any(class(pheno) == "data.frame" | class(pheno) == "matrix" |
            is.null(pheno))){
        stop("pheno must be a data frame, a matrix or NULL")
    }

    if (!any(class(explanatory) != "character" |
            !is.numeric(explanatory))){
        stop("explanatory must be a character or numeric object")
    }

    if (!any(class(covariates) == "character" | class(covariates) == "list" |
            is.null(covariates))){
        stop("covariates must be a character vector, a list or NULL")
    }

    if (!any(class(refGroup) != "character" |
            class(refGroup) != "numeric")){
        stop("refGroup must be a character or numeric object")
    }

    if (!any(class(continuous) == "character" | class(continuous) == "list" |
            is.null(continuous))){
        stop("continuous must be a character vector, a list or NULL")
    }

    if (class(typeInput) != "character"){
        stop("typeInput must be a character object")
    }

    if (class(typeAnalysis) != "character"){
        stop("typeAnalysis must be a character object")
    }

    # Get data from SummarizedExperiment objects

    if (class(methData) == "SummarizedExperiment" |
        class(methData) == "RangedSummarizedExperiment" ){
        if (is.null(pheno)){
            pheno <- SummarizedExperiment::colData(methData)
        }
        methData <- SummarizedExperiment::assay(methData)
    }
    else {
        if (is.null(pheno)) {
            stop("If methData is not a SummarizedExperiment, you must provide
                pheno parameter")
        }
    }

    if (is.null(continuous)){
        continuous <- c()
        categorical <- colnames(pheno)
    }
    else {
        if (class(continuous) != "character") {
            continuous <-colnames(pheno)[continuous]
        }
        categorical <- setdiff(colnames(pheno), continuous)
    }

    # Ensure all categorial variables are factors, continuous variables are
    # numeric and there are no lists
    for (column in colnames(pheno)) {
        if (column %in% categorical) {
            if (typeof(pheno[,column]) == "list") {
                message(paste(column, "variable skipped due to it is a list"))
                pheno[, -which(names(pheno) == column)]
            }
            else {
                pheno[,column] <- factor(pheno[,column])
            }
        }
        else {
            pheno[,column] <- as.numeric(as.character(pheno[,column]))
        }
    }

    typeInput <- match.arg(typeInput, c("beta", "M"))
    typeAnalysis <- match.arg(typeAnalysis, c("M", "beta"))


    if (class(explanatory) == "numeric") {
        explanatory <- colnames(pheno)[explanatory]
    }

    if (is.numeric(covariates)) {
        covariates <- colnames(pheno)[covariates]
    }

    if (class(refGroup) == "numeric") {
        refGroup <- levels(pheno[,explanatory])[refGroup]
    }

    if (length(intersect(explanatory, covariates)) > 0) {
        stop("You specified some variable(s) as both explanatory and covariate")
    }

    pheno <- data.frame(pheno[,c(explanatory, covariates)])
    colnames(pheno) <- c(explanatory, covariates)

    # Prepare methylation data for limma
    if (typeInput == "beta") {
        if (any(min(methData, na.rm=TRUE) < 0 | max(methData, na.rm=TRUE) > 1)){
            warning("Introduced beta-values are not between 0 and 1. Are you
                    sure these are not M-values?")
        }

        if (typeAnalysis == "beta") {
            dataLimma <- methData
        }
        else {
            message("Transforming beta-values to M-values")
            dataLimma <- log2(methData) - log2(1 - methData)
        }
    }
    else {
        if (min(methData, na.rm=TRUE) >= 0 && max(methData, na.rm=TRUE) <= 1) {
            warning("Introduced M-values are between 0 and 1. Are you sure these
                    are not beta-values?")
        }

        if (typeAnalysis == "beta") {
            message("Transforming M-values to beta-values")
            dataLimma <- 2^(methData)/(1 + 2^(methData))
        }
        else {
            dataLimma <- methData
        }
    }

    # Perform linear model
    message("Calculating linear model...")
    message(paste("\tExplanatory variable:", explanatory))

    if (is.factor(pheno[,explanatory])){
        message(paste("\tReference group:", refGroup))
        pheno[,explanatory] <- relevel(pheno[,explanatory], ref=refGroup)
    }

    if (is.null(covariates)){
        message("\tCovariates: None")
        message(paste("\tCategorical variables:",
                    paste(categorical, collapse=" ")))
        if (length(continuous) > 0) {
            message(paste("\tContinuous variables:",
                            paste(continuous, collapse=" ")))
        }
        else {
            message("\tContinuous variables: None")
        }

        model <- model.matrix(~get(explanatory), data=pheno)
        }
    else {
        message(paste("\tCovariates:", paste(covariates, collapse=" ")))
        message(paste("\tCategorical variables:",
                    paste(categorical, collapse=" ")))
        message(paste("\tContinuous variables:",
                    paste(continuous, collapse=" ")))
        model <- model.matrix(~., data=pheno)
    }

    linearModel <- limma::eBayes(limma::lmFit(dataLimma, model))

    tValues <- linearModel[["t"]][,2]
    return(tValues)
}
