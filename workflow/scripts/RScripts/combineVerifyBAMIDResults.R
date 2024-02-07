library(dplyr)
### VerifyBAM Output
args <- commandArgs(trailingOnly = TRUE)
user.input.1 <- args[1]
user.input.2 <- args[2]
sample_files <- scan(user.input.1, character())
# mywd<-"/Users/jaina13/myPART/"
# sample_files <- list.files(path = mywd,pattern = ".Ancestry",recursive = F, full.names = T)
fileNames <- unlist(lapply(sample_files, function(currfile) {
  unlist(stringr::str_split(fs::path_file(currfile), pattern = ".Ancestry"))[1]
}))
names(sample_files) <- fileNames
verifyBAMIdResults <- lapply(sample_files, function(currfile) {
  r <- read.table(currfile, header = T, sep = "\t")
  t <- r[, 3, drop = FALSE]
  colnames(t) <- names(currfile)
  return(t)
})
v <- do.call("cbind", verifyBAMIdResults)
corRes <- cor(v)
finalSamplePairsCorrelation <- cor(v) %>%
  as.data.frame() %>%
  dplyr::mutate(var1 = rownames(.)) %>%
  tidyr::gather(var2, value, -var1) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::group_by(value) %>%
  dplyr::filter(var1 != var2) # %>% dplyr::filter(row_number()==1)
colnames(finalSamplePairsCorrelation) <- c("Sample1", "Sample2", "VerfiyBAMId:Correlation")
write.table(finalSamplePairsCorrelation, file = user.input.2, sep = "\t", quote = F, row.names = F)
