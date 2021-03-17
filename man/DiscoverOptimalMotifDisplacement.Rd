\name{DiscoverOptimalMotifDisplacement}
\alias{DiscoverOptimalMotifDisplacement}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DNA Methylation Motif Discovery
}
\description{
DNA Methylation Motif finding algorithm
}
\usage{
DiscoverOptimalMotifDisplacement(Species, Growing.Mode = "Central", Clust.Metric = "Cosine", InputDataFileName = "ffipsc1911", Num.CPU = 10, Prone.Threshold = 0.95, Resistant.Threshold = 0.15, Min.MotifLength = 6, Max.MotifLength = 15, SeqRepetition.Frequency = 5, CG.Displacement = 0, Coord.Displacement = 1, FullInputDataFile, PathReferenceSequence, PathOutput = NA, VecDisplacements = seq(1,100,1))
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
  \item{FullInputDataFile}{
  FullInputDataFile Parameter
}
  \item{PathReferenceSequence}{
  PathReferenceSequence Parameter
}
  \item{PathOutput}{
  PathOutput Parameter
}
  \item{VecDisplacements}{
  Vector Displacements
}
}
\details{
CG centered Motifs
}
\value{
Motifs
}
\author{
Inigo Reiriz and Marcos Arauzo-Bravo
}
\note{
Package for finding methylation motifs
}
