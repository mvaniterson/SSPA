\name{Nutrigenomics}
\alias{Nutrigenomics}
\docType{data}
\title{Test statistics from a Nutrigenomics gene expression profiling experiment}
\description{There are five sets of test statistics each represents a different compound and exposure time. Test statistics were obtained by using an empirical Bayes linear model.   
}
\usage{data(Nutrigenomics)}
\format{
  A data frame with 16539 test statistics for five experiments.
  \describe{
    First row indicates the effective sample size of the experiment.
    Column names refer to the compound and exposure time (see details).
  }
}
\details{
In this experiment the outcome of specific PPAR-alpha activation on murine small intestinal gene expression was examined 
using Affymetrix GeneChip Mouse 430 2.0 arrays. PPAR-alpha was activated by several PPAR-alpha-agonists that differed in activating potency. 
In this paper the data of three agonists were used, namely Wy14,643, fenofibrate and trilinolenin (C18:3). The first two compounds belong to the 
fibrate class of drugs that are widely prescribed to treat dyslipidemia, whereas trilinolenin is an agonist frequently found in the human diet. 
For intestinal PPAR-alpha, Wy14,643 is the most potent agonist followed by C18:3 and fenofibrate. Since time of exposure also affects the effect size, intestines were collected 6 hrs (all three agonists) or 5  days (Wy14,643 and fenofibrate only) after exposure.
}
\source{
 van Iterson, M. 't Hoen, P.A.C. Pedotti, P. Hooiveld, G.J.E.J. den Dunnen, J.T. van Ommen, G.J.B. Boer, J.M. Menezes, R.X., Relative power and  sample size analysis on gene expression profiling data, BMC Genomics, (2009).  
}
%\references{}
\examples{
data(Nutrigenomics)
str(Nutrigenomics)
}
\keyword{datasets}
