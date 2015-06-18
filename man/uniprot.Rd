\name{uniprot}
\alias{uniprot}
\alias{data}
\docType{data}
\title{
Database from DMarker
}
\description{
Database from DMarker \cr
}
\usage{data(uniprot)}
\format{
  The format is:
 num [1:43462, 1:5] ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:43462] "1" "2" "3" "4" ...
  ..$ : chr [1:5] "Protein" "Gene" "blood" "urine" ...
}
\examples{
data(uniprot);
}
\author{
Yu Shang (JLU & UGA) \email{yushang@uga.edu} \cr 
Qiong Yu (JLU & UGA) \email{yuqiong@uga.edu} \cr 
Wei Du (JLU & UGA) \email{weidu@uga.edu} \cr 
Ying Xu (JLU & UGA) \email{xyn@bmb.uga.edu} \cr \cr 
Maintainer: Yu Shang (JLU & UGA) \email{yushang@uga.edu} \cr 
}
\references{
\code{[1]} Yan Wang, et al. (2009) \emph{DMarker: A Bio-Marker Inference System for Human Diseases based
on Microarray Gene Expression Data} 2009 \cr
\code{[2]} Juan Cui, et al. (2008) \emph{Computational prediction of human proteins that can be secreted into the bloodstream} BIOINFORMATICS, Vol.24 no. 20 2008 pages 2370-2375 \cr
\code{[3]} Jiaxin Wang, et al. (2013) \emph{Computational Prediction of Human Salivary Proteins from Blood Circulation and Application to Diagnostic Biomarker Identification} PLoS ONE, DOI: 10.1371/journal.pone.0080211, 2013 \cr
\code{[4]} S Hong, et al. (2011) \emph{A Computational Method for Prediction of Excretory Proteins and Application to Identification of cancer markers in urine and application to gastric cancer} PLoS ONE,6(2):e16875, 2011 \cr
\code{[5]} \url{http://bioinfosrv1.bmb.uga.edu/DMarker/} \cr
}
\keyword{ datasets }
\keyword{ matrix }
\keyword{ uniprot }
\seealso{
\code{\link{DMarker-package}} \code{\link{DMarker-method}} \cr
}
