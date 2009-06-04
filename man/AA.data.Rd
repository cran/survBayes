\name{AA.data}
\alias{AA.data}
\docType{data}
\title{Shrinkage of aneurisms}
\description{
  Data on the shrinkage of
aneurisms associated with cerebral arteriovenous malformations
(cAVM) after treatment. The time to a shrinkage of the aneurism to
below 50\% 
of the baseline volume was of interest. There is one random inspection 
time (current status). The data is given in interval notation.
Several patients had multiple aneurisms. 
}
\usage{data(AA.data)}
\format{
  A data frame with 149 observations on the following 6 variables.
  \describe{
    \item{z}{censoring variable}
    \item{mo}{the degree of cAMV occlusion by embolization (dichotomized at 50\%
    )}
    \item{lok}{the location of the aneurism, whether at the midline arteries
or at other afferent cerebral arteries}
    \item{gr}{The single aneurisms are not independent because aneurisms within
a patient may shrink in the same way. Multiple aneurisms were
observed per patient. This clustering of aneurisms is indicated by
this grouping variable.}
    \item{t.left}{time of the begin of the interval}
    \item{t.right}{time of the end of the interval}
  }
}
\source{
H.~J. Meisel, U.~Mansmann, H.~Alvarez, G.~Rodesch, M.~Brock, and
P.~Lasjaunias.
Cerebral arteriovenous malformations and associated aneurysms:
  Analysis of 305 cases from a series of 662 patients.
\emph{Neurosurgery}, 46:793--802, 2000.
}
\examples{
data(AA.data)
}
\keyword{datasets}
