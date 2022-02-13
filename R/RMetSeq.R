#' Estimation of the level of methylation
#'
#' @param bamFiles Vector with bam file names
#' @param GR GRanges object with coordinates of restriction sites
#' @param tread The number of running threads at the same time, doesn't work in windows.
#' @param duplicate Remove duplicates
#' @return Data frame with the level of methylation and coverage by restriction sites
#' @import Biostrings
#' @import GenomicRanges
#' @import bamsignals
#' @import parallel
#' @export



RMetSeq<-function(bamFiles,GR,tread = 1, duplicate = T){

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
    Profile<-Profile[seq(1,nrow(Profile),2),]+Profile[seq(2,nrow(Profile),2),]
    Profile<-Profile[,1]+Profile[,2]
    Metil<-Profile/Coverage
    Metil[Metil>1]<-1
    Metil[Metil<0]<-0
    Metil<-(-1*Metil)+1
    Metil2<-cbind(round(Coverage,0),Metil)
    colnames(Metil2)<-c("Coverage","M_level")
    row.names(Metil2)<-paste0(seqnames(GR),".",start(GR))

    return(data.frame(Metil2,stringsAsFactors = F))

  }


  if (is.numeric(tread)) {

    tread <- ceiling(tread)

  } else {

    stop("The threads must be a numeric value")

  }


  if(tread == 1)
  {

    Methyl_List<-lapply(X = bamFiles, FUN = function(x){Methylation(x,GR,duplicate)})
    Methyl_List<-do.call(cbind,Methyl_List)
  }

  if (tread > 1){

    Methyl_List <-mclapply(X = bamFiles, FUN = function(x){Methylation(x,GR,duplicate)}, mc.cores = tread)
    Methyl_List<-do.call(cbind,Methyl_List)
  }

  if (length(names(bamFiles))<length(bamFiles)){
    warning("Sample names are missing, they will be assigned automatically")
    names(bamFiles)<- 1:length(bamFiles)
  }
  colnames(Methyl_List)[seq(1,ncol(Methyl_List),2)]<-paste0(names(bamFiles),".","Coverage")
  colnames(Methyl_List)[seq(2,ncol(Methyl_List),2)]<-paste0(names(bamFiles),".","M_level")

  return(Methyl_List)
}









