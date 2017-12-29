test_mCSEAIntegrate <- function() {
    library(leukemiasEset)
    data(leukemiasEset)
    ALLexpr <- Biobase::exprs(leukemiasEset
                            [,leukemiasEset[["LeukemiaType"]] == "ALL"])[,1:10]
    NoLexpr <- Biobase::exprs(leukemiasEset
                            [,leukemiasEset[["LeukemiaType"]] == "NoL"])[,1:10]
    exprTest <- cbind(ALLexpr, NoLexpr)
    colnames(exprTest) <- 1:20
    resultsInt <- mCSEAIntegrate(myResults, exprTest, "promoters", "ENSEMBL", "GATA2",
                                    makePlot = FALSE)
    checkEqualsNumeric(resultsInt[1,4],
                    -0.8908771, tolerance=1.0e-6)
}
