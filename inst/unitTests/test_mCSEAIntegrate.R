test_mCSEAIntegrate <- function() {
    library(leukemiasEset)
    data(exprTest)
    resultsInt <- mCSEAIntegrate(myResults, exprTest, "promoters", "ENSEMBL", "GATA2",
                                    makePlot = FALSE)
    checkEqualsNumeric(resultsInt[1,4],
                    -0.8908771, tolerance=1.0e-6)
}
