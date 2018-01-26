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

here <- function() {
	args <- commandArgs(trailingOnly = FALSE)
	needle <- "--file="
	match <- grep(needle, args)
	file <-
		if (length(match) > 0) # Rscript
			sub(needle, "", args[match])
		else # 'source'd via R console
			sys.frames()[[1]]$ofile
	normalizePath(dirname(file))
}

parse_results <- function(results, res_loc, run_id) {
	count.dir <- glue("{res_loc}/counts")
	dir.create(path = count.dir, recursive = TRUE, showWarnings = FALSE)
	
	for (r in names(results)) {
		res <- results[[r]]
		res$cell.id <- paste(res$left.bc, res$right.bc, sep = "-")
		m <- unlist(strsplit(r, "_"))
		m2 <- tolower(m[grep("lib*\\d+", m, ignore.case = TRUE)])
		lib_info <- read.delim(glue("{res_loc}/sample_desc{m2}_info.txt"))
		lib_info$cell.id <- glue_data(lib_info, "{left}-{right}")
		
		res2 <- res[which(res$cell.id %in% lib_info$cell.id), ] # filter out the combinations that do not exist in the library
		res3 <- t(table(res2[, c(1, 4)]))
		res3 <- data.frame(cbind(row.names(res3), res3))
		names(res3)[1] <- "cell.id"
		
		res4 <- merge(lib_info, res3, by = "cell.id", all.x = TRUE)
		
		write.table(res4, file = glue("{count.dir}/{r}_count_matrix_{run_id}.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
		
		res2 <- res[which(!(res$cell.id %in% lib_info$cell.id)), ] # filter out the combinations that do not exist in the library
		res3 <- t(table(res2[, c(1, 4)]))
		res3 <- data.frame(cbind(row.names(res3), res3))
		names(res3)[1] <- "cell.id"
		write.table(res3, file = glue("{count.dir}/{r}_count_matrix_extrabc_{run_id}.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
	}
}

# call the functions for processing
call_mapping <- function(lib_loc, res_loc, bc_splitter, bc_mismatch, run_blast, eval_cut, run_id, queue) {
	get_qc(lib_loc, res_loc) 
	merge_paired_ends(lib_loc, res_loc)
	generate_inserts(res_loc, run_blast)
	get_barcode_combinations(res_loc)
	files <- blast_reads(res_loc, bc_splitter, bc_mismatch, run_blast, queue)
	files <- list.files(glue("{res_loc}/processed_data"), pattern = "extendedFrags", full.names = TRUE)
	# split basename with "."
	m <- regexpr(".+/", files)
	# names the files with the representation
	names(files) <- sapply(strsplit(substring(files, attr(m, "match.length") + 1), ".", fixed = TRUE), function(x) x[1])
	results <- tabulate_mapping_results(files, lib_loc, res_loc, eval_cut) 
	names(results) <- names(files)
	save(results, file = glue("{res_loc}/allres_{run_id}.RData"))
	load(glue("{res_loc}/alllibs_res_{run_id}.RData"))
	names(results) <- sub("-", "", sub("NGS\\d+-", "", names(results)))
	parse_results(results, res_loc, run_id)
}

HERE <- here()

options(java.parameters = "-Xmx8000m")
registerDoMC(cores = 4)

### expected
# res_loc/barcodes/barcodes.txt
# res_loc/amplicons/amplicons.txt
# res_loc/sample_desc/library info
# res_loc/raw_data
res_loc <- "data/ngs21"

# other scripts
source(glue("{HERE}/src/fastq2count.R"))

call_mapping(
	lib_loc = glue("{res_loc}/rawdata"),
	res_loc = res_loc,
	bc_splitter = "{HERE}/src/barcode_splitter.pl",
	bc_mismatch = 0,
	run_blast = TRUE,
	eval_cut = 1e-35,
	run_id = "e35_bcmm0",
	queue = FALSE)
