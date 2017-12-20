suppressMessages({library(Biostrings)
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
})
options(java.parameters = "-Xmx8000m")
registerDoMC(cores=4)

# parse results
parseRes <- function (results, res.loc, run.id) {
  count.dir <- paste0(res.loc, "/counts/")
  if(!file.exists(count.dir)) {
    dir.create(path = count.dir, recursive = T)
  }
  for (i in c(1:length(results))) {
    res= results[[i]]
    res$cell.id <- paste(res$left.bc, res$right.bc, sep="-")
    m <- unlist(strsplit(names(results)[i], "_"))
    lib.info <- read.delim(paste0(res.loc, "/sample_desc/", tolower(m[grep("lib*\\d+", m, ignore.case = T)]), "_info.txt"))
    lib.info$cell.id <- paste0(lib.info$left, "-", lib.info$right)
    
    res2 <- res[(which(res$cell.id%in%lib.info$cell.id==T)), ] # filter out the combinations that do not exist in the library
    res3 <- t(table(res2[, c(1, 4)]))
    res3 <- data.frame(cbind(row.names(res3), res3))
    names(res3)[1] <- "cell.id"
    
    res4 <- merge(lib.info, res3, by="cell.id", all.x=T)
    
    write.table(res4,file=paste0(count.dir,"/" ,names(results)[i], "_count_matrix_", run.id, ".txt"), sep="\t", col.names=T, row.names=F, quote=F)
    
    res2 <- res[(which(res$cell.id%in%lib.info$cell.id==F)), ] # filter out the combinations that do not exist in the library
    res3 <- t(table(res2[, c(1, 4)]))
    res3 <- data.frame(cbind(row.names(res3), res3))
    names(res3)[1] <- "cell.id"
    write.table(res3,file=paste0(count.dir,"/" ,names(results)[i], "_count_matrix_extrabc_", run.id, ".txt"), sep="\t", col.names=T, row.names=F, quote=F)
  }
}
# call the functions for processing
callMapping <- function(lib.loc, res.loc, bc.splitter, bc.missmatch, run.blast, eval.cut, run.id){
 getQC(lib.loc, res.loc) 
 mergePairedEnds(lib.loc, res.loc)
 generateInserts (lib.loc, res.loc, run.blast)
 getBarcodeCombinations (lib.loc, res.loc)
 files <- blastReads (lib.loc, res.loc, bc.splitter, bc.missmatch, run.blast)
 files = list.files(paste0(res.loc, "processed_data"),pattern = "extendedFrags",full.names = T)
 # split basename with "."
 m=regexpr(".+/",files)
 # names the files with the representation
 names(files) = sapply(strsplit(substring(files,attr(m,"match.length")+1),".",fixed=T),function(x)x[1])
 results <- tabulateMappingResults (files, lib.loc, res.loc, eval.cut) 
 names(results) <- names(files)
 save(results, file=paste0(res.loc, "allres_", run.id, ".RData"))
 load(paste0(res.loc, "alllibs_res_", run.id, ".RData"))
 names(results) <- sub("-", "", sub("NGS\\d+-", "", names(results)))
 parseRes(results, res.loc, run.id)
}

### expected
# res.loc/barcodes/barcodes.txt
# res.loc/amplicons/amplicons.txt
# res.loc/sample_desc/library info
# raw data
ngs.id <- "ngs21"

# flags
run.blast = T
queue = T
run.id="e35_bcmm0"

# other scripts
source("/src/fastq2count.R")
bc.splitter <- "/src/barcode_splitter.pl"
# parameters 
bc.missmatch <- 0
eval.cut <- 1e-35
# data locations
lib.loc <- paste0("/", ngs.id, "/rawdata/")
res.loc <- paste0("/", ngs.id, "/")
if(!file.exists(res.loc)) {
  dir.create(path = res.loc, recursive = T)
}
callMapping (lib.loc, res.loc, bc.splitter, bc.missmatch, run.blast, eval.cut, run.id)
