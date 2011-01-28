\name{ord_fcn}
\alias{ord_fcn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ord_fcn function gets the acenstor-descendants relationship.}
\description{
The function read the output from eNewick2phylog and then convert the acenstor-descendants relationship in an array format.
}
\usage{
ord_fcn()
}
%- maybe also 'usage' for other objects documented here.
\details{
##  ~~ If necessary, more details than the description above ~~
}
\value{
##  ~Describe the value returned
##  If it is a LIST, use
##  \item{comp1 }{Description of 'comp1'}
##  \item{comp2 }{Description of 'comp2'}
## ...
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
function () 
{
    bmtp <- matrix(rev(res$nde), ncol = 1)
    rvpt <- rev((res$parts))
    rept <- array(0, c(length(rvpt), 2))
    for (i in 1:length(rvpt)) {
        rept[i, ] = unlist(rvpt[i])
    }
    cmb <- cbind(bmtp, rept)
    brnlen <- res$droot[(length(tipnames) + 1):length(res$droot)]
    root <- matrix(cmb[1, ], nrow = 1)
    cmb <- cmb[-1, ]
    brnlen <- brnlen[1:(length(brnlen) - 1)]
    new_ord <- order(brnlen, decreasing = TRUE)
    cmb <- cmb[new_ord, ]
    cmb <- rbind(root, cmb)
    return(cmb)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

