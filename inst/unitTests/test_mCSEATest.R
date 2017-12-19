test_mCSEATest <- function() {
    set.seed(123)
    checkEqualsNumeric(mCSEATest(myRank[1:100], "p", platform = "E")[[1]]["PKN3","ES"],
                    -0.6900556, tolerance=1.0e-6)
}
