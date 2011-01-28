\name{pair_fcn}
\alias{pair_fcn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
pair_fcn 
}
\description{
pair_fcn is an assist function that return pair for tmp sequences.}
\usage{
pair_fcn(tmp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tmp}{
%%     ~~Describe \code{tmp} here~~
}
}
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
function (tmp) 
{
    numl = comb(length(tmp), 2)
    count = 0
    posit <- array(0, c(numl))
    for (i in 1:length(tmp)) {
        for (j in 1:length(tmp)) {
            if (i < j) {
                count = count + 1
                posit[count] = paste(c(tmp[i]), c(tmp[j]), sep = ",")
            }
        }
    }
    return(posit)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

