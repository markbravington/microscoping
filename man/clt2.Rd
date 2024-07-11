\name{clt2}
\alias{clt2}
\title{Sex-aware close-kin scoping}
\description{Like \code{\link{ckmr_laugh_test}} (qv) but everything is split by sex, for better handling of strongly dimorphic species. This two-sex version \code{clt2} entails considerably more tedium when setting up the input data, and runs about eight times slower, so stick to the single-sex version \code{\link{ckmr_laugh_test}} if you can bear to. It's all very rough anyway. And I haven't really tested it. Feel free to do so by simulation (not sarcastic!).
}
\usage{
clt2( n_sa, w_sa, pmat_sa, nsamp_sya= NULL, samp_p_sa,
  nsamp_y, MAX_JUVE_AGE, MIN_CATCH_CURVE_AGE= 0, beqy_POP= FALSE,
  etriv= 0.01, extraplus= FALSE, want_gory_details= FALSE,
  .BLOCKSIZE= 1e4, .REPORTSEC= 5)
}
\arguments{
\item{ n_sa,  w_sa, pmat_sa, samp_p_sa}{all as per \code{ckmr_laugh_test}, except now sex-specific; so they should each be a two-row matrix or \code{\link{offarray}}, one row for each sex. Rownames, if present, must be "F" and "M" in either order, but must be consistent. With \code{offarray} arguments, I \emph{think} that AGES \emph{must} start at 0, but I'm not sure this is checked...}
\item{ nsamp_sya}{Sample sizes by sex, year, age (3D \code{offarray}, or possibly normal \code{array}). A value of NULL, the default, means it will be constructed from \code{samp_p_sa[S,A]*nsamp_y[Y]} in the obvious notation.}
\item{ nsamp_y, MAX_JUVE_AGE, MIN_CATCH_CURVE_AGE, beqy_POP, etriv, extraplus, want_gory_details}{as per \code{ckmr_laugh_test}}
\item{ .BLOCKSIZE, .REPORTSEC}{passed to \code{\link{autoloop}}. If the main loop is slow, you can set these to get a progress update every \code{.REPORTSEC} seconds. I wouldn't recommend messing with \code{.BLOCKSIZE} but it exists for a reason--- and so does the \code{autoloop} documentation ;).}
}
\value{
Mostly as per \code{\link{ckmr_laugh_test}}. Results are not disaggregated by sex (though the computations internally are sex-specific) because the focus is on total kin-pairs.
Attributes \code{call}, \code{params}, and \code{nsamp_sya} are attached. The last two won't display in full if \code{mvbutils} has been explicitly loaded, because they are class \code{nullprint}. That is entirely deliberate. If you want to see them, use \code{unclass} but that will zap their offarray-ness too; safer is eg
%#
\item{  res <- clt2(...)}{}
\item{  nsampo <- res@nsamp_sya # requires library( atease); already imported}{}
\item{  oldClass( nsampo) <- oldClass( nsampo)[-1]}{}
\item{  # can now print nsampo}{}
}
\examples{
# This really SHOULD have an example
# but it doesn't yet
}
\keyword{misc}
