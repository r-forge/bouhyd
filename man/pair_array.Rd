\name{pair_array}
\alias{pair_array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove  "" from pair-fcn }
\description{
}
\usage{
pair_array(tmp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tmp}{
%%     ~~Describe \code{tmp} here~~
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
function (tmp) 
{
    pair <- pair_fcn(tmp)
    p_arr <- matrix(c(0), nrow = length(pair), ncol = 2)
    for (i in 1:length(pair)) {
        p_arr[i, 1] = unlist(strsplit(pair[i], ","))[1]
        p_arr[i, 2] = unlist(strsplit(pair[i], ","))[2]
    }
    return(p_arr)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
