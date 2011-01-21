\name{cov_mtx}
\alias{cov_mtx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cov_mtx returns the variance-covariance matrix to describe the relationship among species.
}
\description{
cov_mtx returns the variance-covariance matrix to describe the relationship among species.
}
\usage{
cov_mtx(bt)
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
    ins_fcn <- function(ist, sqc) {
        ist <- as.numeric(unlist(strsplit(ist, "X"))[2])
        arr <- array(0, c(length(otmp)))
        for (i in 1:length(arr)) {
            arr[i] <- as.numeric(unlist(strsplit(sqc[i], "X"))[2])
        }
        insp <- which(arr == (ist - 1)) + 1
        return(insp)
    }
    var_fcn <- function() {
        for (i in 1:length(otmp)) {
            newi <- which(rownames(mtx) \%in\% otmp[i])
            oldi <- which(rownames(omtx) \%in\% otmp[i])
            mtx[newi, newi] <- omtx[oldi, oldi]
        }
        prn1 <- tmp[which(tmp \%in\% ins) - 1]
        prn2 <- tmp[which(tmp \%in\% ins) + 1]
        prn1 <- which(rownames(omtx) \%in\% prn1)
        prn2 <- which(rownames(omtx) \%in\% prn2)
        vhii <- bt^2 * (omtx[prn1, prn1] + omtx[prn2, prn2] + 
            2 * omtx[prn1, prn2])
        hii <- which(!(tmp \%in\% otmp))
        mtx[hii, hii] <- vhii
        return(mtx)
    }
    fillspcmtx <- function() {
        elm <- function() {
            ck <- c(tmp[nsi], tmp[nsj])
            for (i in 1:dim(pn_arr)[1]) {
                if (sum(pn_arr[i, ] == ck) == 2) {
                  break
                }
            }
            return(i)
        }
        pn_arr <- pair_array(tmp)
        po_arr <- pair_array(otmp)
        nsi <- which(!(tmp \%in\% otmp))[1]
        nsj <- which(!(tmp \%in\% otmp))[2]
        osii <- which(!(otmp \%in\% tmp))
        mtx[nsi, nsj] <- omtx[osii, osii]
        pn_arr <- pn_arr[-elm(), ]
        while (length(pn_arr[, 1]) > 0) {
            newi <- which(rownames(mtx) \%in\% pn_arr[1, 1])
            newj <- which(rownames(mtx) \%in\% pn_arr[1, 2])
            if (tmp[nsi] \%in\% pn_arr[1, ]) {
                otg <- which(!(pn_arr[1, ] \%in\% tmp[nsi]))
                oldi <- which(rownames(omtx) \%in\% otmp[osii])
                oldj <- which(rownames(omtx) \%in\% pn_arr[1, otg])
            }
            if (tmp[nsj] \%in\% pn_arr[1, ]) {
                otg <- which(!(pn_arr[1, ] \%in\% tmp[nsj]))
                oldi <- which(rownames(omtx) \%in\% otmp[osii])
                oldj <- which(rownames(omtx) \%in\% pn_arr[1, otg])
            }
            if (!(tmp[nsi] \%in\% pn_arr[1, ]) && !(tmp[nsj] \%in\% 
                pn_arr[1, ])) {
                oldi <- which(rownames(omtx) \%in\% pn_arr[1, 1])
                oldj <- which(rownames(omtx) \%in\% pn_arr[1, 2])
            }
            mtx[newi, newj] <- omtx[oldi, oldj]
            pn_arr <- pn_arr[-1, ]
            if (length(pn_arr) == 2) {
                pn_arr <- matrix(pn_arr, nrow = 1)
            }
        }
        mtx <- mtx + t(mtx)
        mtx[nsi, nsi] <- omtx[osii, osii] + branchlength[length(tmp) - 
            1]
        mtx[nsj, nsj] <- omtx[osii, osii] + branchlength[length(tmp) - 
            1]
        dianew <- which(tmp \%in\% otmp)
        diaold <- which(otmp \%in\% tmp)
        for (i in 1:length(dianew)) {
            mtx[dianew[i], dianew[i]] <- omtx[diaold[i], diaold[i]] + 
                branchlength[length(tmp) - 1]
        }
        return(mtx)
    }
    fillhydmtx <- function() {
        pn_arr <- pair_array(tmp)
        while (length(pn_arr[, 1]) > 0) {
            newi <- which(rownames(mtx) \%in\% pn_arr[1, 1])
            newj <- which(rownames(mtx) \%in\% pn_arr[1, 2])
            if (ins \%in\% pn_arr[1, ]) {
                otg <- pn_arr[1, which(!(pn_arr[1, ] \%in\% ins))]
                otgj <- which(rownames(omtx) \%in\% otg)
                prn1 <- tmp[which(tmp \%in\% ins) - 1]
                prn2 <- tmp[which(tmp \%in\% ins) + 1]
                prn1 <- which(rownames(omtx) \%in\% prn1)
                prn2 <- which(rownames(omtx) \%in\% prn2)
                mtx[newi, newj] <- bt * (omtx[prn1, otgj] + omtx[prn2, 
                  otgj])
            }
            else {
                oldi <- which(rownames(omtx) \%in\% pn_arr[1, 1])
                oldj <- which(rownames(omtx) \%in\% pn_arr[1, 2])
                mtx[newi, newj] <- omtx[oldi, oldj]
            }
            pn_arr <- pn_arr[-1, ]
            if (length(pn_arr) == 2) {
                pn_arr <- matrix(pn_arr, nrow = 1)
            }
        }
        return(mtx)
    }
    ckins <- FALSE
    rept <- data[, 2:3]
    bmtp <- matrix((data)[, 1], ncol = 1)
    loop <- 2
    tmp = array(0, c(loop))
    if (loop == 2) {
        tmp = rept[1, ]
        otmp <- tmp
        mtx <- diag(branchlength[1], c(length(tmp)))
        rownames(mtx) <- c(tmp)
        colnames(mtx) <- c(tmp)
        omtx <- mtx
    }
    while (loop < length(bmtp)) {
        loop <- loop + 1
        tmp = array(0, c(length(otmp) + 1))
        mtx <- matrix(0, nrow = length(tmp), ncol = length(tmp))
        q = loop - 1
        op <- which(otmp == bmtp[q])
        if (length(op) != 0) {
            tmp[op:(op + 1)] = rept[q, ]
            if (op == 1) {
                tmp[(op + 2):length(tmp)] = otmp[(op + 1):(length(tmp) - 
                  1)]
            }
            if ((op + 1) == length(tmp)) {
                tmp[1:(op - 1)] = otmp[1:(op - 1)]
            }
            if (op != 1 && (op + 1) != length(tmp)) {
                tmp[(op + 2):length(tmp)] = otmp[(op + 1):(length(tmp) - 
                  1)]
                tmp[1:(op - 1)] = otmp[1:(op - 1)]
            }
            rownames(mtx) <- c(tmp)
            colnames(mtx) <- c(tmp)
            mtx <- fillspcmtx()
            otmp <- tmp
            omtx <- mtx
        }
        else {
            ins <- (bmtp[q])
            insp <- ins_fcn(ins, otmp)
            tmp[insp] <- ins
            tmp[(insp + 1):length(tmp)] = otmp[insp:(length(tmp) - 
                1)]
            tmp[1:(insp - 1)] = otmp[1:(insp - 1)]
            rownames(mtx) <- c(tmp)
            colnames(mtx) <- c(tmp)
            diamtx <- var_fcn()
            mtx <- fillhydmtx()
            mtx <- mtx + t(mtx) + diamtx
            otmp <- tmp
            omtx <- mtx
            ckins <- TRUE
        }
        if (ckins) {
            tmp <- array(0, c(length(tmp)))
            tmp[which(otmp == ins)] <- rept[loop - 1, 1]
            tmp[which(otmp == bmtp[loop])] = rept[loop, which(rept[loop, 
                ] != ins)]
            tmp[which(otmp == bmtp[loop + 1])] = rept[loop + 
                1, which(rept[loop + 1, ] != ins)]
            tx1 <- which(otmp == ins)
            tx2 <- which(otmp == bmtp[loop])
            tx3 <- which(otmp == bmtp[loop + 1])
            for (i in 1:length(tmp)) {
                if (i != tx1 && i != tx2 && i != tx3) {
                  tmp[i] = otmp[i]
                }
            }
            otmp <- tmp
            rownames(mtx) <- c(tmp)
            colnames(mtx) <- c(tmp)
            mtx <- mtx + diag(branchlength[length(tmp) - 1], 
                c(length(tmp)))
            omtx <- mtx
            ckins <- FALSE
            loop <- loop + 2
        }
    }
    if (sum(tipnames \%in\% tmp) != nleaves) {
        tmp <- tipnames
        mtx <- matrix(0, nrow = length(tmp), ncol = length(tmp))
        rownames(mtx) <- c(tmp)
        colnames(mtx) <- c(tmp)
        mtx <- fillspcmtx()
    }
    return(mtx)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
