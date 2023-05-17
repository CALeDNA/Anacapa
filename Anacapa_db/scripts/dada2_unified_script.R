#!/usr/bin/env Rscript
# Interpret the input variables -----

##### command arguments for running script ADD Soon!
## the dada2 commands come directly from https://benjjneb.github.io/dada2/tutorial.html

args = commandArgs(trailingOnly=TRUE)

barC = args[1]  #barcode target
odirpath = args[2]  #path to the fastq files
paired_or_not = args[3] # type of reads- should be "paired", "forward", or "reverse

# confirm that the user has specified paired_or_not properly
if (!(paired_or_not %in% c("paired", "forward", "reverse"))) {
  cat("Please specify sequence type as 'paired', 'forward', or 'reverse'")
  quit()
}

## path to output

if (paired_or_not == "paired") {
  path = paste(odirpath,"/", barC, "/", barC, "_sort_by_read_type/paired",  sep='')
  mergedoutpath=paste(odirpath,"/", barC, "/", barC, "dada2_out", sep='')
  unmergedoutpath=paste(odirpath,"/", barC, "/", barC, "dada2_out", sep='')
} else if(paired_or_not == "forward") {
  path = paste(odirpath,"/", barC, "/", barC, "_sort_by_read_type/unpaired_F", sep='')
  outpath=paste(odirpath,"/", barC, "/", barC, "dada2_out", sep='')
} else {
  path = paste(odirpath,"/", barC, "/", barC, "_sort_by_read_type/unpaired_R", sep='')
  outpath=paste(odirpath,"/", barC, "/", barC, "dada2_out", sep='')

}

# Confirm that the user has Write access to the path
if (file.access(path, mode = 2) != 0) {
  stop("Please make sure that you have write access to the supplied path.")
}

library("dada2")

# Set up paths to files ----------

list.files(path)

if(paired_or_not == "paired") {
  fnFs <- sort(list.files(path, pattern="_Paired_1_pairs_R1.fastq"))
  fnRs <- sort(list.files(path, pattern="_Paired_2_pairs_R2.fastq"))
  all_sample_names <- sapply(strsplit(fnFs, "_Paired"), `[`, 1)
  fnFs <- file.path(path, fnFs)
  fnRs <- file.path(path, fnRs)

} else {
  fnFs <- sort(list.files(path, pattern="_Paired_\\d_singles.fastq"))
  all_sample_names <- sapply(strsplit(fnFs, "_Paired_\\d_singles.fastq"), `[`, 1)
  fnFs <- file.path(path, fnFs)
}

# Make the path to which filtered sequences should be outputted ---------
filt_path <- file.path(path, "filtered")

if(paired_or_not == "paired") {
  filtered_seqs_name <- file.path(filt_path, paste0(all_sample_names, "_F_filt.fastq.gz"))
  filtered_seqs_name_R <- file.path(filt_path, paste0(all_sample_names, "_R_filt.fastq.gz"))
} else if (paired_or_not == "forward") {
  filtered_seqs_name <- file.path(filt_path, paste0(all_sample_names, "_F_filt.fastq.gz"))
} else {
  filtered_seqs_name <- file.path(filt_path, paste0(all_sample_names, "_R_filt.fastq.gz"))
}

# Run the filtering step ------
if(paired_or_not == "paired") {
  filtered_seqs <- filterAndTrim(fnFs, filtered_seqs_name, fnRs, filtered_seqs_name_R, minLen = 10,
                                 maxN=0, maxEE=c(2,2), truncQ=0, rm.phix=TRUE,
                                 compress=F,matchIDs=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
} else {
  filtered_seqs <- filterAndTrim(fnFs, filtered_seqs_name, minLen = 10,
                                 maxN=0, maxEE=c(2), truncQ=0, rm.phix=TRUE,
                                 compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

}
head(filtered_seqs)
