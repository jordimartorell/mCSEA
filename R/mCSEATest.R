#' mCSEA core analysis
#'
#' Perform a methylated CpG sites enrichment analysis in predefined genomic
#' regions
#'
#' @param rank A named numeric vector with the ranking statistic of each CpG
#' site
#' @param methData A data frame or a matrix containing Illumina's CpG probes in
#' rows and samples in columns. A SummarizedExperiment object can be used too
#' @param pheno A data frame or a matrix containing samples in rows and
#' covariates in columns. If NULL (default), pheno is extracted from the
#' SummarizedExperiment object
#' @param column The column name or number from pheno used to split the samples
#' into groups (first column is used by default)
#' @param regionsTypes A character or character vector indicating the predefined
#'  regions to be analyzed. NULL to skip this step and use customAnnotation.
#' @param customAnnotation An optional list with the CpGs associated to each
#' feature (default = NULL)
#' @param minCpGs Minimum number of CpGs associated to a region. Regions below
#' this threshold are not tested
#' @param nproc Number of processors to use in GSEA step (default = 1)
#' @param nperm Number of permutations to do in GSEA step (default = 10000)
#' @param platform Platform used to get the methylation data (450k or EPIC)
#'
#' @return A list with the results of each of the analyzed regions. For each
#' region type, a data frame with the results and a list with the probes
#' associated to each region are generated. In addition, this list also contains
#' the input methData, pheno and platform objects
#'
#' @author Jordi Martorell Marugán, \email{jordi.martorell@@genyo.es}
#'
#' @seealso \code{\link{rankProbes}}, \code{\link{mCSEAPlot}},
#' \code{\link{mCSEAPlotGSEA}}
#'
#' @references Subramanian, A. et al (2005). \emph{Gene set enrichment analysis:
#'  A knowledge-based approach for interpreting genome-wide expression profiles}
#'  . PNAS 102, 15545-15550.
#'
#' @examples
#' \dontrun{
#' library(mCSEAdata)
#' data(mcseadata)
#' myRank <- rankProbes(betaTest, phenoTest, refGroup = "Control")
#' set.seed(123)
#' myResults <- mCSEATest(myRank, betaTest, phenoTest,
#' regionsTypes = "promoters", platform = "EPIC")
#' }
#' data(precomputedmCSEA)
#' head(myResults[["promoters"]])
#' head(myResults[["promoters_association"]])
#' @export


mCSEATest <- function(rank, methData, pheno = NULL, column = 1,
                        regionsTypes = c("promoters", "genes", "CGI"),
                        customAnnotation = NULL, minCpGs = 5, nproc = 1,
                        nperm = 10000, platform = "450k")
    {

    output <- list()

    # Check input objects
    if (!any(class(rank) == "numeric" | class(rank) == "integer")){
        stop("rank must be a numeric vector")
    }

    if (!typeof(rank) == "double"){
        stop("rank must be a named vector")
    }

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

    if (!identical(colnames(methData), rownames(pheno))) {
        if (setdiff(colnames(methData),  rownames(pheno)) == 0 &&
            setdiff( rownames(pheno), colnames(methData)) == 0) {
            pheno <- pheno[colnames(methData),]
        }

        else {
            stop("Sample labels of methData and pheno must be the same")
        }
    }


    if (!any(class(column) != "character" |
             !is.numeric(column))){
        stop("column must be a character or numeric object")
    }

    if (class(regionsTypes) != "character" & !is.null(regionsTypes)){
        stop("regionsTypes must be a character")
    }

    if (is.null(regionsTypes) & is.null(customAnnotation)){
        stop("Either regionsTypes or customAnnotations must be specified")
    }

    if (!any(class(customAnnotation) == "list" | is.null(customAnnotation))){
        stop("customAnnotation must be a list or NULL")
    }

    if (class(platform) != "character"){
        stop("platform must be a character object")
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


    platform <- match.arg(platform, c("450k", "EPIC"))

    if (platform == "450k") {
        assocPromoters <- mCSEAdata::assocPromoters450k
        assocGenes <- mCSEAdata::assocGenes450k
        assocCGI <- mCSEAdata::assocCGI450k

    }
    else {
        assocPromoters <- mCSEAdata::assocPromotersEPIC
        assocGenes <- mCSEAdata::assocGenesEPIC
        assocCGI <- mCSEAdata::assocCGIEPIC
    }

    if (!is.null(regionsTypes)){
        regionsTypes <- match.arg(regionsTypes,
                                choices=c("promoters","genes","CGI"),
                                several.ok=TRUE)

        for (region in regionsTypes) {
            if (region == "promoters") {

                resGSEA <- .performGSEA(region, rank, platform, assocPromoters,
                                        minCpGs, nperm, nproc)

                output[["promoters"]] <- resGSEA[[1]]
                output[["promoters_association"]] <- resGSEA[[2]]

            }

            else if (region == "genes") {

                resGSEA <- .performGSEA(region, rank, platform, assocGenes,
                                        minCpGs, nperm, nproc)

                output[["genes"]] <- resGSEA[[1]]
                output[["genes_association"]] <- resGSEA[[2]]

            }
            else if (region == "CGI") {

                resGSEA <- .performGSEA(region, rank, platform, assocCGI,
                                        minCpGs, nperm, nproc)

                output[["CGI"]] <- resGSEA[[1]]
                output[["CGI_association"]] <- resGSEA[[2]]

            }
        }

    }

    if (!is.null(customAnnotation)) {

        resGSEA <- .performGSEA(region, rank, platform, customAnnotation,
                                minCpGs, nperm, nproc)

        output[["custom"]] <- resGSEA[[1]]
        output[["custom_association"]] <- resGSEA[[2]]
    }

    output[["methData"]] <- methData
    output[["pheno"]] <- data.frame(Group = factor(pheno[,column]),
                                    row.names = rownames(pheno))
    output[["platform"]] <- platform
    return(output)
}


.performGSEA <- function(region, rank, platform, assoc, minCpGs, nperm, nproc) {

    if (region != "custom"){
        message(paste("Associating CpG sites to", region))

        if (platform == "450k" & length(rank) > 242500 |
            platform == "EPIC" & length(rank) > 433000){
            dataDiff <- setdiff(unlist(assoc), names(rank))
            genes <- lapply(assoc,
                            function(x) {x[!x %in% dataDiff]})
        }

        else {
            genes <- lapply(assoc,
                            function(x) {x[x %in% names(rank)]})
        }
    }
    else {
        genes <- assoc
    }

    message(paste("Analysing", region))

    fgseaRes <- fgsea::fgsea(genes, rank, minSize=minCpGs,
                            nperm=nperm, nproc=nproc)
    fgseaDataFrame <- as.data.frame(fgseaRes)
    rownames(fgseaDataFrame) <- fgseaDataFrame[,1]
    fgseaDataFrame <- fgseaDataFrame[,-1]
    message(paste(sum(fgseaDataFrame[["padj"]] < 0.05),
                "DMRs found (padj < 0.05)"))
    fgseaSorted <- fgseaDataFrame[order(fgseaDataFrame[["NES"]]),]
    fgseaSorted[,7] <- sapply(fgseaSorted[,7],
                            function(x) {paste(unlist(x),
                                                collapse=", ")})

    return(list(fgseaSorted, genes))
}
