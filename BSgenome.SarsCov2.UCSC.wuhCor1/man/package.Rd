\name{BSgenome.SarsCov2.UCSC.wuhCor1}
\docType{package}

\alias{BSgenome.SarsCov2.UCSC.wuhCor1-package}
\alias{BSgenome.SarsCov2.UCSC.wuhCor1}
\alias{SarsCov2}

\title{Full genome sequence for SARS-Cov-2 (UCSC version wuhCor1)}

\description{
  Contains sequence of Wuhan Coronavirus.
}

\details{
  
}

\note{
  This BSgenome data package was made from the following source data files:
  \preformatted{
wuhCor1.2bit from http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/bigZips
  }

  See \code{?\link[BSgenome]{BSgenomeForge}} and the BSgenomeForge
  vignette (\code{vignette("BSgenomeForge")}) in the \pkg{BSgenome}
  software package for how to make a BSgenome data package.
}

\author{Aditya Bhagwat}

\seealso{
  \itemize{
    \item \link[BSgenome]{BSgenome} objects and the
          \code{\link[BSgenome]{available.genomes}} function
          in the \pkg{BSgenome} software package.
    \item \link[Biostrings]{DNAString} objects in the \pkg{Biostrings}
          package.
    \item The BSgenomeForge vignette (\code{vignette("BSgenomeForge")})
          in the \pkg{BSgenome} software package for how to make a BSgenome
          data package.
  }
}

\examples{
BSgenome.SarsCov2.UCSC.wuhCor1
genome <- BSgenome.SarsCov2.UCSC.wuhCor1
head(seqlengths(genome))


## ---------------------------------------------------------------------
## Genome-wide motif searching
## ---------------------------------------------------------------------
## See the GenomeSearching vignette in the BSgenome software
## package for some examples of genome-wide motif searching using
## Biostrings and the BSgenome data packages:
if (interactive())
    vignette("GenomeSearching", package="BSgenome")
}

\keyword{package}
\keyword{data}
