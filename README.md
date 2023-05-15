# RMetSeq

RMetSeq package for processing targeted MRE-seq The evaluation is based on a comparison of the number of cut and uncut reads in the sites of methyl sensitive restrictases.

## Installation

You can install RMetSeq from GitHub using devtools:

```R
# Install devtools if not already installed
install.packages("devtools")
library(devtools)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("bamsignals")
BiocManager::install("GenomicRanges")
BiocManager::install("Biostrings")


# Install RMetSeq
install_github("alekseizarubin/RMetSeq")
```
## Usage
```R
# Load RMetSeq and other necessary libraries
library("RMetSeq")
library(GenomicAlignments)
library(Rsamtools)
library(bamsignals)
library(msgbsR)

#Selection of coordinates of restriction sites
genome<- readDNAStringSet(filepath = "path/genome.fa",format = "fasta")
Coordinates_CCGG<-vmatchPattern("CCGG", genome,
                        max.mismatch=0, min.mismatch=0,
                        with.indels=FALSE, fixed=TRUE,
                        algorithm="auto")

Coordinates_CCGG<-GRanges(Coordinates_CCGG,strand = "*")
start(Coordinates_CCGG)<-start(Coordinates_CCGG)+1
end(Coordinates_CCGG)<-end(Coordinates_CCGG)-1


bamFiles<-dir("path_bam/",pattern = ".bam",full.names = T)

#Exclusion of restriction sites with insufficient coverage
GR<-GRCovFiltration(bamFiles, Coordinates_CCGG,  tread = 20,  mincov = 20,  minSample = 1,  duplicate = T)

#Estimation of the level of methylation
RMetSeq(bamFiles, GR, tread = 1, duplicate = T)

# Or in methylKit format
metKitExport(bamFiles, GR, mincov = 30, save_dir = "./", duplicate = T)
```
