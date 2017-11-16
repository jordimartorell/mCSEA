#' mCSEA core analysis
#'
#' Perform a methylated CpG sites enrichment analysis in predefined genomic
#' regions
#'
#' @param rank A named numeric vector with the ranking statistic of each CpG
#' site
#' @param regionsTypes A character or character vector indicating the predefined
#'  regions to be analyzed. NULL to skip this step
#' @param customAnnotation An optional list with the CpGs associated to each
#' feature (default = NULL)
#' @param minCpGs Minimum number of CpGs associated to a region. Regions below
#' this threshold are not tested
#' @param nproc Number of processors to use in GSEA step (default = 1)
#' @param nperm Number of permutations to do in GSEA step (default = 100000)
#' @param platform Platform used to get the methylation data (450k or EPIC)
#'
#' @return A list with the results of each of the analyzed regions. For each
#' region type, a data frame with the results and a list with the probes
#' associated to each region are generated
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
#' myResults <- mCSEATest(myRank, regionsTypes = "promoters",
#' platform = "EPIC")
#' }
#' data(precomputedmCSEA)
#' head(myResults$promoters)
#' head(myResults$promoters_association)
#' @export


mCSEATest <- function(rank, regionsTypes = c("promoters", "genes", "CGI"),
                    customAnnotation = NULL, minCpGs = 5, nproc = 1,
                    nperm = 100000, platform = "450k")
    {

    output <- list()

    # Check input objects
    if (!any(class(rank) == "numeric" | class(rank) == "integer")){
        stop("rank must be a numeric vector")
    }

    if (!typeof(rank) == "double"){
        stop("rank must be a named vector")
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
        stop("typeInput must be a character object")
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

                message("Associating CpG sites to promoters")

                if (platform == "450k" & length(rank) > 242500 |
                    platform == "EPIC" & length(rank) > 433000){
                    dataDiff <- setdiff(unlist(assocPromoters), names(rank))
                    genes <- lapply(assocPromoters,
                                    function(x) {x[!x %in% dataDiff]})
                }

                else {
                    genes <- lapply(assocPromoters,
                                    function(x) {x[x %in% names(rank)]})
                }


                message("Analysing promoters")

                set.seed(123)
                fgseaRes <- fgsea::fgsea(genes, rank, minSize=minCpGs,
                                        nperm=nperm, nproc=nproc)
                fgseaDataFrame <- as.data.frame(fgseaRes)
                rownames(fgseaDataFrame) <- fgseaDataFrame[,1]
                fgseaDataFrame <- fgseaDataFrame[,-1]
                message(paste(sum(fgseaDataFrame$padj < 0.05),
                            "differentially methylated promoters found",
                            "(padj < 0.05)"))
                fgseaSorted <- fgseaDataFrame[order(fgseaDataFrame$NES),]
                fgseaSorted[,7] <- sapply(fgseaSorted[,7],
                                        function(x) {paste(unlist(x),
                                                        collapse=", ")})

                output$promoters <- fgseaSorted
                output$promoters_association <- genes

            }

            else if (region == "genes") {
                message("Associating CpG sites to gene bodies")

                if (platform == "450k" & length(rank) > 242500 |
                    platform == "EPIC" & length(rank) > 433000){
                    dataDiff <- setdiff(unlist(assocGenes), names(rank))
                    genes <- lapply(assocGenes,
                                    function(x) {x[!x %in% dataDiff]})
                }

                else {
                    genes <- lapply(assocGenes,
                                    function(x) {x[x %in% names(rank)]})
                }


                message("Analysing gene bodies")

                set.seed(123)
                fgseaRes <- fgsea::fgsea(genes, rank, minSize=minCpGs,
                                        nperm=nperm, nproc=nproc)
                fgseaDataFrame <- as.data.frame(fgseaRes)
                rownames(fgseaDataFrame) <- fgseaDataFrame[,1]
                fgseaDataFrame <- fgseaDataFrame[,-1]
                message(paste(sum(fgseaDataFrame$padj < 0.05), "differentially
                            methylated gene bodies found (padj < 0.05)"))
                fgseaSorted <- fgseaDataFrame[order(fgseaDataFrame$NES),]
                fgseaSorted[,7] <- sapply(fgseaSorted[,7],
                                        function(x) {paste(unlist(x),
                                                        collapse=", ")})

                output$genes <- fgseaSorted
                output$genes_association <- genes

            }
            else if (region == "CGI") {
                message("Associating CpG sites to CpG Islands")

                if (platform == "450k" & length(rank) > 242500 |
                    platform == "EPIC" & length(rank) > 433000){
                    dataDiff <- setdiff(unlist(assocCGI), names(rank))
                    genes <- lapply(assocCGI, function(x) {x[!x %in% dataDiff]})
                }

                else {
                    genes <- lapply(assocCGI,
                                    function(x) {x[x %in% names(rank)]})
                }


                message("Analysing CpG Islands")

                set.seed(123)
                fgseaRes <- fgsea::fgsea(genes, rank, minSize=minCpGs,
                                        nperm=nperm, nproc=nproc)
                fgseaDataFrame <- as.data.frame(fgseaRes)
                rownames(fgseaDataFrame) <- fgseaDataFrame[,1]
                fgseaDataFrame <- fgseaDataFrame[,-1]
                message(paste(sum(fgseaDataFrame$padj < 0.05),
                            "differentially methylated CpG Islands found",
                            "(padj < 0.05)"))
                fgseaSorted <- fgseaDataFrame[order(fgseaDataFrame$NES),]
                fgseaSorted[,7] <- sapply(fgseaSorted[,7],
                                        function(x){paste(unlist(x),
                                                        collapse=", ")})

                output$CGI <- fgseaSorted
                output$CGI_association <- genes

            }
        }

    }

    if (!is.null(customAnnotation)) {

        genes <- customAnnotation

        message("Analysing custom regions")

        set.seed(123)
        fgseaRes <- fgsea::fgsea(genes, rank, minSize=minCpGs, nperm=nperm,
                                nproc=nproc)
        fgseaDataFrame <- as.data.frame(fgseaRes)
        rownames(fgseaDataFrame) <- fgseaDataFrame[,1]
        fgseaDataFrame <- fgseaDataFrame[,-1]
        message(paste(sum(fgseaDataFrame$padj < 0.05),
                    "differentially methylated custom regions found",
                    "(padj < 0.05)"))
        fgseaSorted <- fgseaDataFrame[order(fgseaDataFrame$NES),]
        fgseaSorted[,7] <- sapply(fgseaSorted[,7],
                                function(x) {paste(unlist(x),
                                                collapse=", ")})

        output$custom <- fgseaSorted
        output$custom_association <- genes
    }
    return(output)
}





