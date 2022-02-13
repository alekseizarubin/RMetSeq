#' Exclusion of restriction sites with insufficient coverage
#'
#' @param bamFiles Vector with bam file names
#' @param GR GRanges object with coordinates of restriction sites
#' @param tread The number of running threads at the same time, doesn't work in windows.
#' @param mincov Minimum coverage.
#' @param minSample Minimum number of samples for which the restriction site has the minimum amount of coverage.
#' @param duplicate Remove duplicates
#' @return Filtered Granges object
#' @import Biostrings
#' @import GenomicRanges
#' @import bamsignals
#' @export

GRCovFiltration<- function(bamFiles,GR,tread = 1,mincov = 30,minSample = 1, duplicate = T){

  if (length(bamFiles) < 1){
    stop("Missing bam files")
  }

  Methylation<-function(bamFile,GR,duplicate){

    if (duplicate){
      Coverage<-do.call(rbind,bamCoverage(bampath = bamFile,GR,paired.end = "extend",filteredFlag = 1024)@signals)
      Profile<-do.call(rbind,bamProfile(bampath = bamFile,GR,ss =T,filteredFlag = 1024)@signals)
    } else {

      Coverage<-do.call(rbind,bamCoverage(bampath = bamFile,GR,paired.end = "extend")@signals)
      Profile<-do.call(rbind,bamProfile(bampath = bamFile,GR,ss =T)@signals)

    }
    Coverage[,1]<-Coverage[,1]+Profile[seq(1,nrow(Profile),2),2]
    Coverage[,2]<-Coverage[,2]+Profile[seq(2,nrow(Profile),2),1]
    Coverage<-(Coverage[,1]+Coverage[,2])/2

    return(Coverage)

  }


  if (is.numeric(tread)) {

    tread <- ceiling(tread)

  } else {

    stop("The threads must be a numeric value")

  }


  if(tread == 1)
  {

    Coverage_List<-lapply(X = bamFiles, FUN = function(x){Methylation(x,GR,duplicate)})
    Coverage_List<-do.call(cbind,Coverage_List)
  }

  if (tread > 1){

    Coverage_List <-mclapply(X = bamFiles, FUN = function(x){Methylation(x,GR,duplicate)}, mc.cores = tread)
    Coverage_List<-do.call(cbind,Coverage_List)
  }

  GR<-GR[rowSums(Coverage_List >= mincov) >= minSample]



  return(GR)
}
