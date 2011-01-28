\name{main}
\alias{main}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
main fits a brownian motion model to the given data
}
\description{
main fits a brownian motion model to the given data
}
\usage{
main(x.tre, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.tre}{
A rooted phylogenetic network
}
  \item{Y}{
Phenotypic data for extant species.
}
}
\details{
##  ~~ If necessary, more details than the description above ~~
}
\value{
main returns the maximum likelihood estimators for the mean, vairation and hybrid parameter.
}
\references{
A new comparative analysis for investigating rate
of trait evolution using phylogenetic networks. manuscript. 
}
\author{
Tony Jhwueng<djhwueng@nimbios.org>
}
\note{
##  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
## ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x.tre, Y) 
{
    res <- newick2phylog(x.tre)
    data <- getntn(res)
    branchlength <- getbrnlen(res)
    tipnames <- sort(names(res$droot[which(res$droot == 1)]))
    nleaves <- length(tipnames)
    n <- nleaves
    one <- matrix(1, nrow = n)
    SS <- optimize(NegLogLike, c(-100, 100))
    mle <- MLE(SS$minimum)
    mu_hat <- mle[2]
    sigmasq_hat <- mle[3]
    print(mle)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

