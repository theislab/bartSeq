# author: richa batra 
# modified original code of steffen saas
# last modified: 20 Dec 2017
suppressPackageStartupMessages({
	library(Biostrings)
	library(ShortRead)
	library(parallel)
	library(stringi)
	library(gplots)
	library(reshape)
	library(stringr)
	library(foreach)
	library(doMC)
	library(DESeq2)
	library(plyr)
	library(reshape2)
	library(glue)
})

shell_transformer <- function(code, envir) shQuote(glue::evaluate(code, envir))

# quality control of the reads with FASTQC
get_qc <- function(lib_loc, res_loc) {
	qc_loc <- glue("{res_loc}/qc")
	dir.create(path = qc_loc, recursive = TRUE, showWarnings = FALSE)
	
	reads <- list.files(path = lib_loc, pattern = ".gz", full.names = TRUE)
	for (read in reads) {
		out_name <- basename(sub('\\.fastq\\.gz$', '_fastqc.html', read))
		if (file.exists(file.path(qc_loc, out_name))) {
			message(glue("QC files exist. Skipping creation of {out_name} and *.zip"))
			next
		}
		system2("fastqc", c(read, "-o", qc_loc)) 
	}
}

# merging the overlapping paired reads using FLASH
merge_paired_ends <- function(lib_loc, res_loc) {
	pdata_loc <- glue("{res_loc}/processed_data")
	dir.create(path = pdata_loc, recursive = TRUE, showWarnings = FALSE)
	
	read1s <- list.files(path = lib_loc, pattern = "R1.+gz", full.names = TRUE)
	for (reads_1.fq in read1s) {
		reads_2.fq <- sub("R1", "R2", reads_1.fq)
		out_prefix <- sub("_R1.*", "", basename(reads_1.fq))
		if (file.exists(reads_1.fq) && file.exists(reads_2.fq)) {
			if (file.exists(glue("{pdata_loc}/{out_prefix}.hist"))) {
				message(glue("Paired files exist. Skipping creation of {pdata_loc}/{out_prefix}.*"))
				next
			}
			system2("flash", paste(reads_1.fq, reads_2.fq, "-M", "150", "-o", out_prefix, "-d", pdata_loc))
		}
	}
}

# formatting the reference amplicons to be used as blast reference database
generate_inserts <- function(res_loc, run_blast) {
	am_loc <- glue("{res_loc}/amplicons")
	sample_desc <- read.delim(glue("{am_loc}/amplicons.txt"))
	fasta <- glue_data(sample_desc, "> {Gene}\n{toupper(Amplicon)}")
	write(fasta, sep = "\n", file = glue("{am_loc}/inserts.txt"))
	inserts <- readDNAStringSet(glue("{am_loc}/inserts.txt"))
	# blast to match each insert to human genes
	dir.create(path = glue("{res_loc}/blast"), recursive = TRUE)
	if (run_blast) {
		writeXStringSet(inserts, glue("{res_loc}/blast/inserts.fa"), width = 10000)
		blast_out <- try(system2(
			"makeblastdb", c(
				"-dbtype", "nucl",                          # nucleotides
				"-in", glue("{res_loc}/blast/inserts.fa"),  # Input file
				"-out", glue("{res_loc}/blast/inserts")     # index file basename
			), stdout = TRUE))
	}
}

# barcodes to be searched forward and reverse complementary
get_barcode_combinations <- function(res_loc) {
	bc_loc <- glue("{res_loc}/barcodes")
	bcs <- read.delim(glue("{bc_loc}/barcodes.txt"), header = TRUE)
	names(bcs) <- c("ids", "seqs")
	write.table(data.frame(id = bcs$ids, seq = as.character(bcs$seqs)), file = glue("{bc_loc}/bctable_bol.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
	write.table(data.frame(id = bcs$ids, seq = as.character(reverseComplement(DNAStringSet(bcs$seqs)))), file = glue("{bc_loc}/bctable_eol.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# blast the reads and in parallel perform barcode identification
blast_reads <- function(res_loc, bc_splitter, bc_mismatch, run_blast, queue) {
	dir.create(path = glue("{res_loc}/blast"), recursive = TRUE)
	# read in fastq
	files <- list.files(glue("{res_loc}/processed_data"), pattern = "extendedFrags", full.names = TRUE)
	# split basename with "."
	m <- regexpr(".+/", files)
	# names the files with the representation
	names(files) <- sapply(strsplit(substring(files, attr(m, "match.length") + 1), ".", fixed = TRUE), function(x) x[1])
	# foreach file run blast
	foreach(id = names(files)) %dopar% {
		cat(id, "\n")
		if (run_blast) {
			path_reads <- glue("{res_loc}/blast/{id}_reads_sh.fa")
			path_inserts <- glue("{res_loc}/blast/inserts")
			path_res_inserts <- glue("{res_loc}/blast/{id}_res_inserts.txt")
			script <- glue(
				"cat {files[id]}",
					"| awk '{{",
						"if (NR%4==1) print \">\"((NR-1)/4)+1;",
						"if (NR%4==2) print",
					"}}' > {path_reads}\n",
				"blastn",  # https://www.ncbi.nlm.nih.gov/books/NBK279684/
					"-reward 1 -gapopen -1 -gapextend -1",
					"-outfmt 6",                  # alignment view: tabular
					"-db {path_inserts}",         # Database file
					"-query {path_reads}",        # Query file
					"-out {path_res_inserts}\n",  # Report output file
				.sep = " ",
				.transformer = shell_transformer)
			cat(script)
			
			blast_out <- try(
				if (queue) system2("qsub", c("-o", ".", "-e", ".", "-q", "long@@core", "-hard", "-l", "job_mem=40G"), input = script, stdout = TRUE)
				else system2("bash", c("-c", script), stdout = TRUE))
		}
		#barcode matching from the read
		bc_loc <- glue("{res_loc}/barcodes")
		dir.create(path = glue("{bc_loc}/bc_identification"), recursive = TRUE, showWarnings = FALSE)
		cat("running", id)
		bc_ident_out <- try(system2(
			bc_splitter, c(
				"--partial", "3",
				"--eolbcfile", glue("{bc_loc}/bctable_eol.txt"),
				"--bolbcfile", glue("{bc_loc}/bctable_bol.txt"),
				"--mismatches ", bc_mismatch,
				"--o ", glue("{bc_loc}/bc_identification/{id}.txt")),
			input = files[id], stdout = TRUE))
	}
	files
}

# formate the output into count matrix
tabulate_mapping_results <- function(files, lib_loc, res_loc, eval_cut) {
	results <- lapply(names(files), function(id) {
		file <- files[id]
		bc_ident <- read.delim(glue("{res_loc}/barcodes/bc_identification/{id}.txt"), stringsAsFactors = FALSE, header = FALSE, na.strings = "unmatched")
		colnames(bc_ident) <- c("bol", "eol")
		insertmatch <- read.table(glue("{res_loc}/blast/{id}_res_inserts.txt"), stringsAsFactors = FALSE)
		colnames(insertmatch) <- c(
			"query", "subject", "percentage.id", "alignment.length", "mismatches", "gap.openings",
			"query.start", "query.end", "subject.start", "subject.end", "E.value", "bit.score")
		# filter the low quality matches
		insertmatch <- insertmatch[insertmatch$E.value < eval_cut, ]
		insertmatch <- insertmatch[order(insertmatch$E.value), ]
		# if there are more than one match for a query
		insertmatch <- insertmatch[which(!duplicated(insertmatch$query)), ]
		insertmatch <- insertmatch[order(insertmatch$query), ]
		# match to the barcodes given the query numbers are corresponding to the order of the barcodes
		insertmatch$bol <- bc_ident$bol[insertmatch$query] # query hold the oder or reads in for blast and barcode identification  
		insertmatch$eol <- bc_ident$eol[insertmatch$query] # and therefore can be used to map the two tables
		res <- data.frame(
			amplicon = insertmatch$subject,
			left.bc = pmin(insertmatch$eol,insertmatch$bol),
			right.bc = pmax(insertmatch$eol,insertmatch$bol)) # left barcodes minimun in eol and right barcodes max in eol
		# blast file has duplicates of reads and we 
		rownames(res) <- insertmatch$query 
		res <- res[grepl("L[[:digit:]]+", res$left.bc) & grepl("R[[:digit:]]+", res$right.bc), ] # removing the nas
		# list(
		#   fractions = c(
		#   	frac.compl = nrow(res) / nrow(bc_ident),
		#   	frac.insert = nrow(insertmatch) / nrow(bc_ident),
		#   	frac.bol = sum(!is.na(bc_ident$bol)) / nrow(bc_ident),
		#   	frac.eol = sum(!is.na(bc_ident$eol)) / nrow(bc_ident)
		#   ),
		#   result.counts = table(res),
		#   result = res)
	})
	names(results) <- names(files)
	results
}
