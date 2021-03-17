\name{DiscoverMotif}
\alias{DiscoverMotif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DNA Methylation Motif Discovery
}
\description{
DNA Methylation Motif finding algorithm
}
\usage{
DiscoverMotif(Species, Growing.Mode = "Central", Clust.Metric = "Cosine", InputDataFileName = "ffipsc1911", Num.CPU = 10, Prone.Threshold = 0.95, Resistant.Threshold = 0.15, Min.MotifLength = 6, Max.MotifLength = 15, SeqRepetition.Frequency = 5, CG.Displacement = 0, Coord.Displacement = 1, DrawLogo = TRUE, AssignGenes = TRUE, FullInputDataFile, PathReferenceSequence, PathOutput = NA, ThrFDR = NA)
}
\arguments{
  \item{Species}{
  Species Parameter
}
  \item{Growing.Mode}{
  Growing Mode
}
  \item{Clust.Metric}{
  Clust Metric
}
  \item{InputDataFileName}{
  InputDataFileName Parameter
}
  \item{Num.CPU}{
  Num CPU
}
  \item{Prone.Threshold}{
  Prone Threshold
}
  \item{Resistant.Threshold}{
  Resistant Threshold
}
  \item{Min.MotifLength}{
  Min MotifLength
}
  \item{Max.MotifLength}{
  Max MotifLength
}
  \item{SeqRepetition.Frequency}{
  SeqRepetition Frequency
}
  \item{CG.Displacement}{
  CG Displacement
}
  \item{Coord.Displacement}{
  Coord Displacement
}
  \item{DrawLogo}{
  Draw Logo
}
  \item{AssignGenes}{
  Assign Genes
}
  \item{FullInputDataFile}{
  FullInputDataFile Parameter
}
  \item{PathReferenceSequence}{
  PathReferenceSequence Parameter
}
  \item{PathOutput}{
  PathOutput Parameter
}
 \item{ThrFDR}{
  FDR Threshold
}
}
\details{
CG centered Motifs
}
\value{
  \item{ListFDR}{Returned Variable}
}
\author{
Inigo Reiriz and Marcos Arauzo-Bravo
}
\note{
Package for finding methylation motifs
}
