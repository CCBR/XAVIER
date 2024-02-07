library(dplyr)
### VerifyBAM Output
args <- commandArgs(trailingOnly = TRUE)
user.input.1 <- args[1]
# user.input.2=args[2]
user.input.3 <- args[2]

# SomaliaOutput
somaliaDistance <- read.table(user.input.1, sep = "\t", header = F)
# somaliaDistance<-read.table("~/relatedness.pairs.tsv",sep = "\t",header = F)
predictedPairs <- list()
for (sample in sort(unique(c(somaliaDistance$V1, somaliaDistance$V2))))
{
  samplePairs <- somaliaDistance %>% dplyr::filter(V1 %in% sample | V2 %in% sample)
  maxRelatedess <- which(samplePairs$V3 == max(samplePairs$V3, na.rm = T))
  maxHomCon <- which(samplePairs$V6 == max(samplePairs$V6, na.rm = T))
  m <- intersect(maxRelatedess, maxHomCon)
  if (length(m) > 0) {
    for (i in m)
    {
      if (sample == unlist(samplePairs[i, ]$V2)) {
        t <- unlist(samplePairs[i, c(2, 1, 3, 6)])
        names(t) <- NULL
        predictedPairs[[sample]] <- (t)
      } else {
        predictedPairs[[sample]] <- (unlist(samplePairs[i, c(1, 2, 3, 6)]))
      }
    }
  } else {
    print(paste0("No consensous between the Relatedness and Homology for sample:", sample))
    maxVal <- c(maxRelatedess, maxHomCon)
    for (i in maxVal)
    {
      if (sample == unlist(samplePairs[i, ]$V2)) {
        t <- unlist(samplePairs[i, c(2, 1, 3, 6)])
        names(t) <- NULL
        predictedPairs[[sample]] <- (t)
      } else {
        predictedPairs[[sample]] <- (unlist(samplePairs[i, c(1, 2, 3, 6)]))
      }
    }
  }
}
finalPredPairs <- data.frame(do.call("rbind", predictedPairs))
colnames(finalPredPairs) <- c("Sample1", "Sample2", "Som:relatedness", "Som:hom_concordance")

# VerifyBAMID
# verifyBAMID<-read.table(user.input.2,sep = "\t",header = T)
# verifyBAMID<-read.table("~/IntendedSamplesPCs.cor.pairs.tsv",sep = "\t",header = T)
# predictedPairsVerifyBAMID<-list()
# repSam<-c()
# for(sample in sort(unique(c(verifyBAMID$Sample1,verifyBAMID$Sample2))))
# {
#   samplePairs<-verifyBAMID %>% dplyr::filter(Sample1 %in% sample | Sample2 %in% sample)
#   print(sample)
#   maxCorrelation<-which(samplePairs$VerfiyBAMId.Correlation == max(samplePairs$VerfiyBAMId.Correlation,na.rm = T))
#   for(i in maxCorrelation)
#   {
#     if(sample == unlist(samplePairs[i,]$Sample2))
#     {
#       t<-unlist(samplePairs[i,c(2,1,3)])
#       names(t)<-NULL
#       predictedPairsVerifyBAMID[[sample]]<-t
#     }else
#     {
#       predictedPairsVerifyBAMID[[sample]]<-cbind(samplePairs[i,])
#     }
#   }
# }
# finalpredictedPairsVerifyBAMID<-do.call("rbind",predictedPairsVerifyBAMID)

## Combine the output from both the tools
# mergedDF<-merge(x=finalPredPairs,y=finalpredictedPairsVerifyBAMID,by = "Sample1",all = TRUE)
# write.table(mergedDF[,c(1:4,6)],file = user.input.3,sep = "\t",quote = FALSE,row.names = FALSE)
write.table(finalPredPairs, file = user.input.3, sep = "\t", quote = FALSE, row.names = FALSE)
