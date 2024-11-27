test_rankProbes <- function() {
    checkEqualsNumeric(rankProbes(betaTest, phenoTest)[[1]],
                    2.262585, tolerance=1.0e-6)
}
