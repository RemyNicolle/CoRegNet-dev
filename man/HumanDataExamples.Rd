\name{Human Data Examples}
\alias{CIT_BLCA_EXP}
\alias{CIT_BLCA_CNV}
\alias{CIT_BLCA_Subgroup}
\alias{ENCODE_sub}
\alias{CHEA_sub}
\alias{STRING_sub}
\alias{HIPPIE_sub}
\docType{data}
\title{
Human Transcription Factor data and Bladder cancer dataset.
}
\description{
Internal datasets for the use of the package containing :
A list of 2020 Human Transcription Factors in EntrezGene or Official Gene Symbol (HUGO) from the FANTOM consortium (Ravasi et al, 2010),
A gene expression and copy number dataset of bladder cancer,
Datasets of human TF protein interactions and regulatory interactions.
}

\usage{data(HumanTF)
data(HumanTF_entrezgene)
data(CIT_BLCA_EXP)
data(CIT_BLCA_CNV)
data(CIT_BLCA_Subgroup)
data(ENCODE_sub)
data(CHEA_sub)
data(STRING_sub)
data(HIPPIE_sub)
}

\format{
  \describe{
    \item{\code{HumanTF}}{A character vector of Official (HUGO) gene symbols}
    \item{\code{HumanTF_entrezgene}}{A character vector of EntrezGene}
    \item{\code{CIT_BLCA_EXP}}{A matrix of bladder cancer gene expression from the CIT program. Only contains the 1000 genes with highest standard deviation. rownames : genes, colnames : samples. (from Rebouissou et al. 2014)}
    \item{\code{CIT_BLCA_CNV}}{A matrix of bladder cancer DNA copy number aberrations from the CIT program. Only contains the TF and the tumor samples of \code{CIT_BLCA_EXP}. rownames : TF, colnames : samples.}
  \item{\code{CIT_BLCA_Subgroup}}{A data.frame with sample name in first column and classification from the TCGA in the second column.}
    \item{\code{ENCODE_sub}}{A data.frame of regulatory interactions from the ENCODE project.}
    \item{\code{CHEA_sub}}{A data.frame of regulatory interactions from the ChEA2 database. (http://amp.pharm.mssm.edu/ChEA2)}
    \item{\code{STRING_sub}}{A data.frame of Protein interaction between TF from the STRING database. (http://string.embl.de)}
    \item{\code{HIPPIE_sub}}{A data.frame of Protein interaction between TF from the HIPPIE database. (http://cbdm.mdc-berlin.de/tools/hippie)}
       
  }
}

\references{
Rebouissou S, Bernard-Pierrot I, de Reynies A, Lepage M-L, Krucker C, Chapeaublanc E, Herault A, Kamoun A, Caillault A, Letouze E, Elarouci N, Neuzillet Y, Denoux Y, Molinie V, Vordos D, Laplanche A, Maille P, Soyeux P, Ofualuka K, Reyal F, et al (2014)
EGFR as a potential therapeutic target for a subset of muscle-invasive bladder cancers presenting a basal-like phenotype. 
Science Translational Medicine 6: 244ra91-244ra91

Ravasi T, Suzuki H, Cannistraci CV, Katayama S, Bajic VB, Tan K, Akalin A, Schmeier S, Kanamori-Katayama M, Bertin N,
Carninci , Daub 
et al (2010) An Atlas of Combinatorial Transcriptional Regulation in Mouse and Man. Cell 140: 744-752

Kou Y, Chen EY, Clark NR, Duan Q, Tan CM and Ma ayan A (2013) ChEA2: Gene-Set Libraries from ChIP-X Experiments to Decode the Transcription Regulome. In ... pp 416-430. Berlin, Heidelberg: Springer Berlin Heidelberg

Franceschini A, Szklarczyk D, Frankild S, Kuhn M, Simonovic M, Roth A, Lin J, Minguez P, Bork P, Mering von C and Jensen LJ (2012) STRING v9.1: protein-protein interaction networks, with increased coverage and integration. Nucleic Acids Research 41: D808-D815

Schaefer MH, Fontaine J-F, Vinayagam A, Porras P, Wanker EE and Andrade-Navarro MA (2012) HIPPIE: Integrating Protein Interaction Networks with Experiment Based Quality Scores. PLoS ONE 7: e31826

}
\examples{
data(HumanTF)
head(HumanTF)
data(HumanTF_entrezgene)
head(HumanTF_entrezgene)
}
\keyword{datasets}
