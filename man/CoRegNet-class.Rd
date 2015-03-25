\name{coregnet-class}
\alias{coregnet-class}
\docType{class}


\title{Class coregnet --- Specifying the structure of the network used throughout the package.}


\description{
  \code{coregnet} is an S4 type class which specifies the structure containing the necessary data of a Co-Regulatory network. \code{coregnet} objects are constructed by the \code{hLICORN} function which infers a regulatory network from gene expression data. A \code{coregnet} object can also be initiliazed by a network constructed by any other reverse engeneering method as specified in the manual of the \code{coregnet} function.
}


\section{Initialization}{
  An object of type CoRegNet can be either built from a gene expression data set using the \code{hLICORN} function, using a data.frame containing a network specifying the activators and inhibitors of each gene or using an adjacency matrix (see the \link{coregnet} function). 
}


\section{Slots}{
  \describe{
    \item{\code{GRN}:}{A data.frame with three columns, Target Gene, co-activators and co-repressors. Co-regulators are seperated by spaces. }
    \item{\code{adjacencyList}:}{ A description of the network through a adjacency matrix.
      \describe{
        \item{bytf}{A named list of lists. The first entry point, the names of the list, are the Transcription Factors of the network. For each TF a list with two entries, act and rep, contains the set of activated and repressed genes respectively.        
        }
        \item{bygene}{
        A named list of lists. The first entry point, the names of the list, are the target genes of the network. For each gene a list with two entries, act and rep, contains the set of activators and repressors of the genes respectively.
        }
      }
    }
    \item{\code{coRegulators}:}{A data.frame specifying the set of inferred co-regulators with several measures and statistics for each pairs.}
    \item{\code{evidences}:}{A list containing all the added regulatory or co-regulatory evidences if any.}
    \item{\code{evidenceDescription}:}{A list containing all the added regulatory or co-regulatory evidences if any. }
    \item{\code{inferenceParameters}:}{A list of parameters of hLICORN. }
  }
}


\section{Methods}{
  Generic methods to be used with this class :

  \describe{
    \item{print}{ Prints the number of genes, transcription factors and regulatory interactions.}
    \item{summary}{ Same as print}
    \item{length}{ outputs the number of genes, transcription factors and regulatory interactions in the form of a numeric vector}
    \item{dim}{ same as length}
    \item{targets}{Returns the set of target genes or the targets of a given list of regulators}
    \item{regulators}{Returns the set of regulators of the entwork or the set of regulators of a given list of genes.}
    \item{activators}{Same as regulators but returns only activators of a gene}
    \item{repressors}{Same as activtors for repressor regulators}
    \item{as.list}{Gets the network in the form of an adjacency list}
    \item{as.data.frame}{Get the network in the form of an edge list}
    \item{addEvidences}{Enriches the network with external regulation evidences}
    \item{addCooperativeEvidences}{Enriches the network with external evidence of cooperative interaction between regulators}
    \item{refine}{Refine the network based on a score which can take into account external regulation or co-regulation data.}
    \item{regulatorInfluence}{Predicts the influence of the regulators in the network. Returns a predicted regulatory activity sample by sample}
    \item{fitCoregnet}{ Returns the fitness of the regulatory network given an expression dataset}
    \item{display}{Starts a shiny application to display the co-regulation graph }    
  }   
}

\author{
Remy Nicolle <remy.c.nicolle AT gmail.com>
}

\seealso{
  \code{\link{hLICORN}}, \link{coregnet}
}

