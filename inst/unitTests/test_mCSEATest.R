test_mCSEATest <- function() {
    checkEqualsNumeric(mCSEATest(myRank[1:100], "p", platform = "E")[[1]]["PKN3","pval"],
                    0.08897591, tolerance=1.0e-6)
}
