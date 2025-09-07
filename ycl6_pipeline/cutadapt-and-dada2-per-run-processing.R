#

library("dada2")
library("ShortRead")
library("Biostrings")
library("ggplot2")
options(width=190)

# Get script location when running with Rscript
args <- commandArgs(trailingOnly = FALSE)
script.path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
script.dir <- dirname(script.path)

# Change fastq, trimmed, and filt paths accordingly
fastq <- file.path(script.dir, "..", "data", "raw", "20250810", "all_fastq")
fastq <- normalizePath(fastq)  
trimmed <- file.path(script.dir, "..", "data", "processed", "20250810", "trimmed")
filt    <- file.path(script.dir, "..", "data", "processed", "20250810", "filt")

# Create trim/filt directories if they don’t exist
if (!dir.exists(trimmed)) dir.create(trimmed, recursive = TRUE)
if (!dir.exists(filt)) dir.create(filt, recursive = TRUE)

fns = sort(list.files(fastq, full.names = TRUE))
fnFs = fns[grep("1.fq.gz", fns)]
fnRs = fns[grep("2.fq.gz", fns)]


fnFs.cut = file.path(trimmed, basename(fnFs))
fnRs.cut = file.path(trimmed, basename(fnRs))
log.cut = gsub("1.fq.gz", ".log", fnFs.cut)
sample.names = gsub("1.fq.gz", "", basename(fnFs.cut))

# Define the primer set used to perform PCR
FWD = "CCTACGGGNGGCWGCAG"       # 341F
REV = "GACTACHVGGGTATCTAATCC"   # 805R

# Get reverse complement DNA sequences
FWD.RC = dada2::rc(FWD)
REV.RC = dada2::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags = paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags = paste("-G", REV, "-A", FWD.RC)

# Run cutadapt to remove primers
# Note to change the PATH to cutadapt accordingly
cutadapt = "/home/patwuch/miniforge3/envs/16S/bin/cutadapt"

for(i in seq_along(fnFs)) {
	print(paste0("[", i ,"/", length(sample.names), "] ", sample.names[i]))

	system2(cutadapt,
		stdout = log.cut[i], stderr = log.cut[i],	# log file
		args = c(R1.flags, R2.flags,
			"-n 2",							# -n 2 required to remove FWD and REV from reads
			"--match-read-wildcards",		# enable IUPAC nucleotide codes (wildcard characters)
			"--length 300",					# Truncate reads to 300 bp
			"-m 150",						# discard reads shorter than LEN (avoid length zero sequences)
			"--overlap 10",					# min overlap between read and adapter for an adapter to be found
			"-j 0",							# auto-detection of CPU cores, only available on Python 3
			"-o", fnFs.cut[i], "-p", fnRs.cut[i],	# trimmed files
			fnFs[i], fnRs[i])				# input files
	)
}

# Check list of trimmed files
head(list.files(trimmed),15)

fns = sort(list.files(trimmed, full.names = TRUE))
fnFs = fns[grep("1.fq.gz", fns)]
fnRs = fns[grep("2.fq.gz", fns)]
sample.names = gsub(".1.fq.gz", "", basename(fnFs))

# Plot quality profile of fastq files
ii = 1:length(sample.names)
pdf("plotQualityProfile.pdf", width = 8, height = 8, pointsize = 12)
for(i in ii) {
	message(paste0("[", i ,"/", length(sample.names), "] ", sample.names[i]))
	print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd"))
	print(plotQualityProfile(fnRs[i]) + ggtitle("Rev"))
}
dev.off()

# Set paths to the dada2-filterd files
filtFs = file.path(filt, basename(fnFs))
filtRs = file.path(filt, basename(fnRs))

# Perform filtering and trimming
# Review "plotQualityProfile.pdf" to select the best paramters for 'truncLen'
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
	# Need to keep paramters consistent between runs of the same study
	truncLen = c(260,200), minLen = 200, maxN = 0, truncQ = 2, maxEE = c(2,5),
	rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = TRUE)

out = as.data.frame(out)
rownames(out) = sample.names
head(out, 10)

# Dereplication and learn the error rates (Default: nbases = 1e8)
# derepFastq() has been intergrated into learnErrors()
errF = learnErrors(filtFs, multithread = TRUE)
errR = learnErrors(filtRs, multithread = TRUE)

pdf("plotErrors.pdf", width = 10, height = 10, pointsize = 12)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()

# Sample Inference
# By default, the `dada` function processes each sample independently (pool = FALSE)
# Use `pool = TRUE` or `pool = pseudo` (recommended) if samples are from an extremely diverse community (e.g. soil)
dadaFs = dada(filtFs, err = errF, pool = FALSE, multithread = TRUE)
dadaRs = dada(filtRs, err = errR, pool = FALSE, multithread = TRUE)

# Merge paired reads (Default: minOverlap = 12; maxMismatch = 0)
mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Construct sequence table
seqtab = makeSequenceTable(mergers)

# View the length frequency distribution
table(nchar(getSequences(seqtab)))

# Save sequence table
saveRDS(seqtab, "seqtab.rds") # or as an example, use seqtab[c(1:5),] to save data for a subset of the first 5 samples

# Save current workspace
save.image(file = "image1.RData")
