#' Export results in a format compatible with methylKit
#'
#' @param bamFiles Vector with bam file names
#' @param GR GRanges object with coordinates of restriction sites
#' @param mincov Minimum coverage.
#' @param save_dir Directory for saving results.
#' @param duplicate Remove duplicates
#' @import Biostrings
#' @import GenomicRanges
#' @import bamsignals
#' @export



metKitExport<-function(bamFiles,GR,mincov = 30,save_dir = "./", duplicate = T){



  if (length(bamFiles) < 1){
    stop("Missing bam files")
  }

  if (length(names(bamFiles))<length(bamFiles)){
    warning("Sample names are missing, they will be assigned automatically")
    names(bamFiles)<- 1:length(bamFiles)
  }

  for (i in 1:length(bamFiles)){
    bamFile<-bamFiles[i]

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
    keep<-(Coverage >= mincov)
    Coverage<-Coverage[keep]
    Profile<-Profile[seq(1,nrow(Profile),2),]+Profile[seq(2,nrow(Profile),2),]
    Profile<-Profile[keep,1]+Profile[keep,2]
    Metil<-Profile/Coverage
    Metil[Metil>1]<-1
    Metil[Metil<0]<-0
    Metil<-(-1*Metil)+1
    Metil2<-cbind(round(Coverage,0),Metil)
    colnames(Metil2)<-c("Coverage","M_level")
    row.names(Metil2)<-paste0(seqnames(GR)[keep],".",start(GR)[keep])

          mKit<-data.frame(chrBase = row.names(Metil2),
                       chr = seqnames(GR)[keep],
                       base = start(GR)[keep],
                       strand = "F",
                       coverage = Metil2[,1],
                       freqC = round(Metil2[,2] * 100,2),stringsAsFactors = F)
      mKit$freqT<-round(100 - mKit$freqC,2)

      write.table(mKit,file = paste0(save_dir,names(bamFile),"_CpG.txt"),quote = F,sep = "\t",row.names = F)

}






}
