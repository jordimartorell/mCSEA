\name{mCSEA-package}
\alias{mCSEA-package}
\alias{mCSEA}
\docType{package}
\title{
    Methylated CpGs Set Enrichment Analysis
}
\description{
    Identification of diferentially methylated regions (DMRs) in predefined 
    regions (promoters, CpG islands...) from the human genome using
    Illumina's 450K or EPIC microarray data. Provides methods to rank CpG 
    probes based on linear models and includes plotting functions.
}
\author{
Jordi Martorell Marugán

Maintainer: Jordi Martorell Marugán<jordi.martorell@genyo.es>
}

\examples{
\dontrun{
library(mCSEA)
data(mcseadata)
myRank <- rankProbes(betaTest, phenoTest, refGroup = "Control")
myResults <- mCSEATest(myRank, regionsTypes = "promoters", platform = "EPIC")
}
data(precomputedmCSEA)
head(myResults$promoters)
}

