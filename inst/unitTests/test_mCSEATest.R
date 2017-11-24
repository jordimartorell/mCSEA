test_mCSEATest <- function() {
    checkEqualsNumeric(mCSEATest(myRank[1:100], "p", platform = "E")[[1]]["PKN3","ES"],
                    -0.6900556, tolerance=1.0e-6)
}
