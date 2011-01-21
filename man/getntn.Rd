\name{getntn}
\alias{getntn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getntn is an assisting function for function cov_mtx
}
\description{
This function gets rid of unnecessarily "_" symbol. 
}
\usage{
getntn(res)
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
    size <- length(res$parts)
    relarr <- array(0, c(size, 3))
    rvpt <- (res$parts)
    rept <- array(0, c(length(rvpt), 2))
    for (i in 1:length(rvpt)) {
        rept[i, ] = unlist(rvpt[i])
    }
    for (i in 1:size) {
        relarr[i, 1] <- names(res$parts)[i]
    }
    relarr[, 2:3] <- rept
    temp <- matrix(0, row <- size)
    for (j in 2:3) {
        for (i in 1:size) {
            stmp <- unlist(strsplit(relarr[, j][i], "_"))
            temp[i] <- stmp[1]
        }
        relarr[, j] <- temp
    }
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
    new_ord <- order(ndlen)
    relarr <- relarr[new_ord, ]
    return(relarr)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
