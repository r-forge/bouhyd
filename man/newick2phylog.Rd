\name{newick2phylog}
\alias{newick2phylog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
newick2phylog read phylogenetic network in eNewwick format 
}
\description{
newick2phylog read phylogenetic network in eNewwick format }
\usage{
newick2phylog(x.tre, add.tools = FALSE, call = match.call())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.tre}{
x.tre is an network data in eNewick format
}
  \item{add.tools}{
%%     ~~Describe \code{add.tools} here~~
}
  \item{call}{
%%     ~~Describe \code{call} here~~
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
function (x.tre, add.tools = FALSE, call = match.call()) 
{
    lvs_num <- function(x.tree) {
        x.tree <- unlist(strsplit(x.tre, ""))
        lvs <- 0
        for (i in 1:length(x.tree)) {
            if (x.tree[i] == ":") {
                lvs = lvs + 1
            }
        }
        return((lvs + 1)/2)
    }
    complete <- function(x.tre) {
        if (length(x.tre) > 1) {
            w <- ""
            for (i in 1:length(x.tre)) {
                w <- paste(w, x.tre[i], sep = "")
            }
            x.tre <- w
        }
        ndroite <- nchar(gsub("[^)]", "", x.tre))
        ngauche <- nchar(gsub("[^(]", "", x.tre))
        if (ndroite != ngauche) {
            stop(paste(ngauche, "( versus", ndroite, ")"))
        }
        if (regexpr(";", x.tre) == -1) {
            stop("';' not found")
        }
        i <- 0
        kint <- 0
        kext <- 0
        arret <- FALSE
        lvs <- lvs_num(x.tree)
        if (regexpr("\\[", x.tre) != -1) {
            x.tre <- gsub("\\[[^\\[]*\\]", "", x.tre, ext = FALSE)
        }
        x.tre <- gsub(" ", "", x.tre, ext = FALSE)
        while (!arret) {
            i <- i + 1
            if (substr(x.tre, i, i) == ";") {
                arret <- TRUE
            }
            if (substr(x.tre, i, i + 1) == "(,") {
                kext <- kext + 1
                add <- paste("Ext", kext, sep = "")
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
            else if (substr(x.tre, i, i + 1) == ",,") {
                kext <- kext + 1
                add <- paste("Ext", kext, sep = "")
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
            else if (substr(x.tre, i, i + 1) == ",)") {
                kext <- kext + 1
                add <- paste("Ext", kext, sep = "")
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
            else if (substr(x.tre, i, i + 1) == "(:") {
                kext <- kext + 1
                add <- paste("Ext", kext, sep = "")
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
            else if (substr(x.tre, i, i + 1) == ",:") {
                kext <- kext + 1
                add <- paste("Ext", kext, sep = "")
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
            else if (substr(x.tre, i, i + 1) == "),") {
                kint <- kint + 1
                add <- kint + lvs
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
            else if (substr(x.tre, i, i + 1) == "))") {
                kint <- kint + 1
                add <- kint + lvs
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
            else if (substr(x.tre, i, i + 1) == "):") {
                kint <- kint + 1
                add <- kint + lvs
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
            else if (substr(x.tre, i, i + 1) == ");") {
                add <- "Root"
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
        }
        lab.points <- strsplit(x.tre, "[(),;]")[[1]]
        lab.points <- lab.points[lab.points != ""]
        no.long <- (regexpr(":", lab.points) == -1)
        if (all(no.long)) {
            lab.points <- paste(lab.points, ":", c(rep("1", length(no.long) - 
                1), "0.0"), sep = "")
        }
        else if (no.long[length(no.long)]) {
            lab.points[length(lab.points)] <- paste(lab.points[length(lab.points)], 
                ":0.0", sep = "")
        }
        else if (any(no.long)) {
            stop("Non convenient data leaves or nodes with and without length")
        }
        w <- strsplit(x.tre, "[(),;]")[[1]]
        w <- w[w != ""]
        leurre <- make.names(w, unique = TRUE)
        leurre <- gsub("[.]", "_", leurre, ext = FALSE)
        for (i in 1:length(w)) {
            old <- paste(w[i])
            x.tre <- sub(old, leurre[i], x.tre, ext = FALSE)
        }
        w <- strsplit(lab.points, ":")
        label <- function(x) {
            lab <- x[1]
            lab <- gsub("[.]", "_", lab, ext = FALSE)
            return(lab)
        }
        longueur <- function(x) {
            long <- x[2]
            return(long)
        }
        labels <- unlist(lapply(w, label))
        longueurs <- unlist(lapply(w, longueur))
        labels <- make.names(labels, TRUE)
        labels <- gsub("[.]", "_", labels, ext = FALSE)
        w <- labels
        for (i in 1:length(w)) {
            new <- w[i]
            x.tre <- sub(leurre[i], new, x.tre, ext = FALSE)
        }
        cat <- rep("", length(w))
        for (i in 1:length(w)) {
            new <- w[i]
            if (regexpr(paste(")", new, sep = ""), x.tre, ext = FALSE) != 
                -1) 
                cat[i] <- "int"
            else if (regexpr(paste(",", new, sep = ""), x.tre, 
                ext = FALSE) != -1) 
                cat[i] <- "ext"
            else if (regexpr(paste("(", new, sep = ""), x.tre, 
                ext = FALSE) != -1) 
                cat[i] <- "ext"
            else cat[i] <- "unknown"
        }
        return(list(tre = x.tre, noms = labels, poi = as.numeric(longueurs), 
            cat = cat))
    }
    res <- complete(x.tre)
    poi <- res$poi
    nam <- res$noms
    names(poi) <- nam
    cat <- res$cat
    res <- list(tre = res$tre)
    res$leaves <- poi[cat == "ext"]
    names(res$leaves) <- nam[cat == "ext"]
    res$nodes <- poi[cat == "int"]
    names(res$nodes) <- nam[cat == "int"]
    listclass <- list()
    dnext <- c(names(res$leaves), names(res$nodes))
    res$nde <- c(names(res$nodes))
    res$nme <- c(names(res$leaves))
    listpath <- as.list(dnext)
    names(listpath) <- dnext
    x.tre <- res$tre
    print(x.tre)
    while (regexpr("[(]", x.tre) != -1) {
        a <- regexpr("([^()]*)", x.tre, ext = FALSE)
        n1 <- a[1] + 1
        n2 <- n1 - 3 + attr(a, "match.length")
        chasans <- substring(x.tre, n1, n2)
        chaavec <- paste("(", chasans, ")", sep = "")
        nam <- unlist(strsplit(chasans, ","))
        w1 <- strsplit(x.tre, chaavec, ext = FALSE)[[1]][2]
        parent <- unlist(strsplit(w1, "[,);]", ext = FALSE))[1]
        listclass[[parent]] <- nam
        x.tre <- gsub(chaavec, "", x.tre, ext = FALSE)
        w2 <- which(unlist(lapply(listpath, function(x) any(x[1] == 
            nam))))
        for (i in w2) {
            listpath[[i]] <- c(parent, listpath[[i]])
        }
    }
    res$parts <- listclass
    res$paths <- listpath
    dnext <- c(res$leaves, res$nodes)
    names(dnext) <- c(names(res$leaves), names(res$nodes))
    res$droot <- unlist(lapply(res$paths, function(x) sum(dnext[x])))
    res$call <- call
    class(res) <- "phylog"
    if (!add.tools) 
        return(res)
    return(newick2phylog.addtools(res))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
