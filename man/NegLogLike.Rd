\name{NegLogLike}
\alias{NegLogLike}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
NegLogLike is the negative log likelihood function for the model.
}
\description{
The function represents the negative logarithm likelihood function. 
}
\usage{
NegLogLike(bt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bt}{
%%     ~~Describe \code{bt} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
A new comparative analysis for investigating rate
of trait evolution using phylogenetic networks. manuscript. 
}
\author{
Tony Jhwueng<djhwueng@nimbios.org>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (bt) 
{
    W <- cov_mtx(bt)
    solveW <- solve(W)
    mu <- c((t(one) \%*\% solveW \%*\% Y)/(t(one) \%*\% solveW \%*\% 
        one))
    Y <- Y - mu
    sigma_sq <- c(abs((t(Y) \%*\% solveW \%*\% Y)/n))
    neglogML <- n/2 * log(2 * pi) + n/2 * log(sigma_sq) + 1/(2 * 
        sigma_sq) * t(Y) \%*\% solveW \%*\% Y + 1/2 * log(abs(det(W)))
    return(neglogML)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
