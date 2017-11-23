test_mCSEATest <- function() {
    checkEqualsNumeric(mCSEATest(myRank[1:100], "p", platform = "E")[[1]][1,1],
                    0.08897591, tolerance=1.0e-6)
}
