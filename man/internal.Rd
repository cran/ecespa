\name{ecespa-internal}
\alias{check}
\alias{ginv}
\alias{mNNinfo}
\alias{mNNinfo2}
\alias{mNNtest}
\alias{NNid}
\title{ Internal ecespa functions. }
\description{Internal ecespa functions.}
\usage{
check(x, v, l1, l2)
ginv(m)
mNNinfo2(n, R, Q)
mNNinfo(xy, label, nnid = NULL, splancs = TRUE)
mNNtest(info, obsN = NULL)
NNid(xy, splancs = TRUE)
}
\details{
These fucntions are called by \code{\link{dixon2002}} and are usually not to be called by the user.
}
\author{ Philip M. Dixon \url{http://www.public.iastate.edu/~pdixon/}
}
\keyword{ internal }
