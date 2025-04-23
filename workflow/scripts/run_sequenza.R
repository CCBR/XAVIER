args <- commandArgs(trailingOnly = TRUE)

library(sequenza)

if (length(args) == 0) {
  stop("Must provide a seqz file")
} else {
  seqz_file <- args[1]
  if (!file.exists(seqz_file)) {
    stop(paste0("Can't find this SEQZ output file: ", seqz_file))
  }
}

if (length(args) > 1) {
  out_dir <- args[2]
} else {
  out_dir <- dirname(seqz_file)
}

if (length(args) > 2) {
  sampleid <- args[3]
} else {
  sampleid <- gsub(".seqz.gz", "", basename(seqz_file))
}

if (length(args) > 3) {
  n_cores <- as.numeric(args[4])
} else {
  n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
}

if (is.na(n_cores)) {
  n_cores <- 1
}


print(paste0("Using ", n_cores, " cores..."))
date()
print("Extracting seqz data...")
seqzdata <- sequenza.extract(seqz_file, min.reads = 30, min.reads.normal = 10, parallel = n_cores)

date()
print("Fitting model...")
CP.example <- sequenza.fit(seqzdata, mc.cores = n_cores)

## Sequenza.extract seems to fail if too few mutations
num_mutations <- unlist(lapply(seqzdata$mutations, nrow))
chrom_list <- names(num_mutations)[num_mutations > 3]
## But it might actually be segments, idk?
# num_segments <- unlist(lapply(seqzdata$segments, nrow))
# chrom_list <- names(num_mutations)[num_segments > 1]

not_included <- setdiff(names(num_mutations), chrom_list)
print("Printing results...")
if (length(not_included) > 0) {
  print("Excluding these chromosomes because of too few mutations...")
  print(not_included)
}
sequenza.results(sequenza.extract = seqzdata, cp.table = CP.example, sample.id = sampleid, out.dir = out_dir)

date()
print("Done")
