\name{HumanTF}
\alias{HumanTF}
\alias{HumanTF_entrezgene}
\docType{data}
\title{
List of Human Transcription Factors.
}
\description{
A list of 2020 Human Transcription Factors in EntrezGene or Official Gene Symbol (HUGO) from the FANTOM consortium (Ravasi et al, 2010).
}
\usage{data(HumanTF)}
\format{
  Two different objects :
  \describe{
    \item{\code{HumanTF_entrezgene}}{A character vector of EntrezGene}
    \item{\code{HumanTF}}{A character vector of Official (HUGO) gene symbols}
  }
}
\references{
Ravasi T, Suzuki H, Cannistraci CV, Katayama S, Bajic VB, Tan K, Akalin A, Schmeier S, Kanamori-Katayama M, Bertin N,
Carninci , Daub 
et al (2010) An Atlas of Combinatorial Transcriptional Regulation in Mouse and Man. Cell 140: 744-752}

\examples{
data(HumanTF)
head(HumanTF)
data(HumanTF_entrezgene)
head(HumanTF_entrezgene)
}
\keyword{datasets}
