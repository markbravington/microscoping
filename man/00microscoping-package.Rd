\name{microscoping-package}
\alias{microscoping}
\alias{microscoping-package}
\docType{package}
\title{Rough feasiibility calculations for CKMR on typical fish}
\description{
If you don't know what CKMR is, stop reading now.

The main function \code{\link{ckmr_laugh_test}} (qv) predicts how many kin-pairs of various types might be expected for Your species of fish, given Your proposed CKMR sampling scheme and Your guess as to what the current stock size etc might be. The clue is in the name \code{\link{ckmr_laugh_test}} (qv). It is \emph{deliberately} crude because you are not supposed to really trust the results, and certainly not to try fine-tuning them by microtweaking the sampling scheme, because the model assumptions are bound to be somewhat wrong.

I wrote \code{\link{ckmr_laugh_test}} for private use by me & a few chosen colleagues, but it seems to have escaped into the wider world, and predictably to have been used more seriously than I ever intended. Be aware that AFAIK \code{\link{ckmr_laugh_test}} has \emph{never been simulation-tested} and might have substantial errors (I hope not). I don't feel it's my job to test it thoroughly, since I never meant it for general use! Please \bold{feel free} to do that simulation-testing yourself, e.g. using Shane Baylis' \code{fishSim} package or Eric Anderson's \code{ckmrSim} package (and let me know the result); it would be quite helpful..

There is a sex-specific version \code{\link{clt2}}, which I reluctantly added for species with unignorably strong sexual dimorphism (e.g. swordfish) where \code{\link{ckmr_laugh_test}} is really inadequate and cannot be tricked into giving even a semi-reliable answer. \code{\link{clt2}} is harder to use, and even less tested.
}
\keyword{misc}
