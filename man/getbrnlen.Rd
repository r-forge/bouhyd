\name{getbrnlen}
\alias{getbrnlen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getbrnlen gets branch lengths from the network.}
\description{
getbrnlen gets branch lengths from the network for the purpose of biulding the covaraince matrix among species.}
}
\usage{
getbrnlen(res)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
%%     ~~Describe \code{res} here~~
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
function (res) 
{
    ndlen <- res$droot[!(res$droot == 1)]
    nam <- names(ndlen)
    ck1 <- array(0, c(length(nam)))
    count <- 0
    for (ele in c(nam)) {
        count <- count + 1
        len <- length(unlist(strsplit(ele, "_")))
        if (len == 2) {
            ck1[count] <- 1
        }
    }
    ndlen <- ndlen[!ck1]
    ndlen <- sort(ndlen)
    ck2 <- array(0, c(length(ndlen)))
    for (i in 1:(length(ndlen) - 1)) {
        if (abs(ndlen[i] - ndlen[i + 1]) < 10^(-5)) {
            ck2[i] = 1
        }
    }
    ndlen <- ndlen[!ck2]
    brnlen <- array(0, c(length(ndlen)))
    tmplen <- ndlen
    for (i in 1:(length(brnlen) - 1)) {
        brnlen[i] <- tmplen[i + 1] - tmplen[i]
    }
    brnlen[length(brnlen)] <- 1 - tmplen[(length(tmplen))]
    return(brnlen)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
