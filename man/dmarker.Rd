\name{DMarker}
\alias{dmarker}
\alias{DMarker-method}
\alias{DMarker}
\title{
DMarker main function
}
\description{
\code{DMarker} is a package to predict whether the protein or gene can be detected in blood and urine. \cr \cr 
The fuction \code{dmarker} is the main function of the package. It predicts the biomarker in a data vector, and helps to know whether the protein or gene can be detected in blood and urine. \cr \cr
results = dmarker(data, \ldots); \cr \cr
}
\arguments{
\item{data}{ Data vector. }
\item{pre}{ -pre should be ( "gene", "protein" ), is a data type variable, telling DMarker gene or protein you are using. \cr default: "GENE" }
\item{file}{ -file should be a filename to output the results if want to save that, or use "FALSE" to ignore. \cr default: "FALSE" }
\item{type}{ -type should be ( "both", "blood", "urine" ), is a predict type variable, telling DMarker where the data can be detected in blood, urine or both. \cr default: "both" }
\item{Seq}{ -Seq should be ( "T", "F" ). If Seq is True, the DMarker would predic the sequence whether can be detect in blood or urine, and the data vector must be the sequence vector. \cr default: "F" }
}
\value{
The R function, \code{dmarker} returns an object of \code{list}: \cr 
\tabular{ll}{
BLOOD \tab An object of \code{list}, having \code{GENE} and \code{PROTEIN}. \cr 
URINE \tab An object of \code{list}, having \code{GENE} and \code{PROTEIN}. \cr 
GENE \tab A result vector, containing the gene names which can be detected in \code{BLOOD} or \code{URINE}. \cr 
PROTEIN \tab A result vector, containing the protein names which can be detected in \code{BLOOD} or \code{URINE}. \cr 
}
The R function, \code{dmarker} will write the results in a file named by -file, if -file is not FALSE by default. The file will be in the \code{Working Directory} . \cr
}
\details{
results = dmarker(data, pre = "Gene", file = "FALSE", type = "both");
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
\keyword{ DMarker }
\keyword{ blood }
\keyword{ unine }
\seealso{
\code{\link{DMarker-package}} \code{\link{DMarker-method}} \code{\link{uniprot}} \cr
}
\examples{

data(data);
# load the example vector

results = dmarker(data);
# run dmarker with default parameters
# results is a variable of list with,
#   $BLOOD,
#   $BLOOD$GENE, which genes can be detected in BLOOD,
#   $BLOOD$PROTEIN, which proteins can be detected in BLOOD,
#   $URINE,
#   $URINE$GENE, which genes can be detected in BLOOD,
#   $URINE$PROTEIN, which proteins can be detected in URINE,
# or results = dmarker(data, pre = "Gene", file = "FALSE", type = "both");

results = dmarker(data, pre = "Gene", file = "DMarkerResults.txt", type = "both", Seq = "F");
# Saving the results both in R variable results and the file named DMarkerResults.txt by -file.
}
