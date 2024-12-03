test_mCSEATest <- function() {
    set.seed(123)
    checkEqualsNumeric(mCSEATest(myRank[1:100], betaTest, phenoTest, 
                                    regionsTypes = "p", platform = "E")[["promoters"]]["PKN3","ES"],
                    0.6470448, tolerance=1.0e-6)
}
