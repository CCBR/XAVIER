require(dplyr)

args <- commandArgs(trailingOnly = TRUE)
user.input.1 <- args[1]
user.input.2 <- args[2]

somaliaDistance <- read.table(user.input.1, sep = "\t", header = T, comment.char = "")
sampleSexChrMD <- somaliaDistance[, c("sample_id", "X_depth_mean", "Y_depth_mean")]
sampleSexChrMD$Scale_X_depth_mean <- unlist(scale(sampleSexChrMD$X_depth_mean)[, 1])
sampleSexChrMD$Scale_Y_depth_mean <- unlist(scale(sampleSexChrMD$Y_depth_mean)[, 1])
sampleSexChrMD$Gender <- unlist(lapply(sampleSexChrMD$Y_depth_mean, FUN = function(x) {
  if (x > 0) {
    return("Male")
  } else {
    return("Female")
  }
}))
write.table(sampleSexChrMD[, c("sample_id", "Gender")], user.input.2, quote = F, sep = "\t", row.names = F)
