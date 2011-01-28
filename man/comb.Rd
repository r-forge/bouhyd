\name{comb}
\alias{comb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
comb simply returns the combination 
}
\description{
comb simply returns the combination which is a way of selecting several things out of a larger group, where order does not matter.
}
\usage{
comb(n, r)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
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
of trait evolution using phylogenetic networks. manuscript submmitted in Evolution. 
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
function (n, r) 
{
    return(factorial(n)/(factorial(n - r) * factorial(r)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

