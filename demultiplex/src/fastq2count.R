# author: richa batra 
# modified original code of steffen saas
# last modified: 20 Dec 2017
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

# quality control of the reads with FASTQC
getQC <- function (lib.loc, res.loc) {
 qc.loc <- paste0(res.loc, "/qc/")
 if(!file.exists(qc.loc)) {
  dir.create(path = qc.loc, recursive = T)
 }
 reads <- list.files(path=lib.loc, pattern=".gz", full.names=T)
 lapply(1:length(reads), FUN=function(i){
 i.read <- reads[i]
 sys.out <- paste0("/path/to/FASTQC/dir/fastqc ", i.read, " -o ", qc.loc)
 system(sys.out) 
 })
}
# merging the overlapping paired reads using FLASH
mergePairedEnds<- function (lib.loc, res.loc) {
 pdata.loc <- paste0(res.loc, "/processed_data/")
 if(!file.exists(pdata.loc)) {
   dir.create(path = pdata.loc, recursive = T)
 }
 read1s <- list.files(path=lib.loc, pattern="R1.+gz", full.names=T)
 lapply(1:length(read1s), FUN=function(i){
  reads_1.fq <- read1s[i]
  reads_2.fq <- sub("R1", "R2", read1s[i])
  if(file.exists(reads_1.fq) & file.exists(reads_2.fq)) {
   sys.out <- paste0("/path/to/FLASH/dir/flash ", reads_1.fq, " ", reads_2.fq, " -M 150 -o", sub("_R1.*", "", basename(reads_1.fq)), " -d", pdata.loc)
   system(sys.out)
  }
 })
}
# formatting the reference amplicons to be used as blast reference database
generateInserts <- function (lib.loc, res.loc, run.blast) {
 am.loc <- paste0(res.loc, "/amplicons/")
 sample.desc <- read.delim(paste0(am.loc, "amplicons.txt"))
 fasta <- unlist(lapply(1:nrow(sample.desc), FUN=function(i) {
    c(paste0(">",sample.desc$Gene[i]),
      toupper(sample.desc$Amplicon[i]))}))
 write(fasta, sep="\n", file=paste0(am.loc, "inserts.txt"))
 inserts = readDNAStringSet(paste0(am.loc, "inserts.txt"))
 # blast to match each insert to human genes
 dir.create(path = paste0(res.loc, "blast"), recursive = T)
 if(run.blast)
 {
   writeXStringSet(inserts, paste0(res.loc,"blast/inserts.fa"), width=10000)
   cmd = paste0("/path/to/FORMATDB/dir/formatdb -p F -i ", paste0(res.loc, "blast/inserts.fa"), " -n ", paste0(res.loc, "blast/inserts"))
   blast.out=try(system(cmd, intern = TRUE))
 }
}
# barcodes to be searched forward and reverse complementary
getBarcodeCombinations <- function (lib.loc, res.loc) {

 bc.loc <- paste0(res.loc, "/barcodes/")
 bcs <- read.delim(paste0(bc.loc, "barcodes.txt"), header=T)
 names(bcs) <- c("ids", "seqs")
 write.table(data.frame(id=bcs$ids,seq=as.character(bcs$seqs)), file=paste0(bc.loc, "bctable_bol.txt"),sep="\t",quote=F,row.names=F,col.names=F)
 write.table(data.frame(id=bcs$ids,seq=as.character(reverseComplement(DNAStringSet(bcs$seqs)))), file=paste0(bc.loc, "bctable_eol.txt"), sep="\t",quote=F,row.names=F,col.names=F)
}
# blast the reads and in parallel perform barcode identification
blastReads <- function(lib.loc, res.loc, bc.splitter, bc.missmatch, run.blast) {
  # read in fastq
  files = list.files(paste0(res.loc, "processed_data"),pattern = "extendedFrags",full.names = T)
  # split basename with "."
  m=regexpr(".+/",files)
  # names the files with the representation
  names(files) = sapply(strsplit(substring(files,attr(m,"match.length")+1),".",fixed=T),function(x)x[1])
  # foreach file run blast
  foreach (id = names(files)) %dopar%
  {
    cat(id,"\n")
    if(run.blast)
    {
      reads.path = paste0(res.loc, "blast/", id, "_reads_sh.fa")
      if (queue) cmd = paste0("echo \"cat ", files[id]," | awk '{ if (NR%4==1) print \\\">\\\"((NR-1)/4)+1 ; if (NR%4==2) print }' > ",reads.path,"; /path/to/BLAST/dir/blastall -p blastn -r 1 -G -1 -E -1 -m 8 -d ", paste0(res.loc,"blast/inserts")," -i ",reads.path," -o ", paste0(res.loc,"blast/",id,"_res_inserts.txt"),"\" | /path/to/QSUB/dir/qsub -o \".\" -e \".\" -q long@@core -hard -l job_mem=40G")
      cat(cmd)
      blast.out=try(system(cmd, intern = TRUE))
    }
    #barcode matching from the read
    bc.loc <- paste0(res.loc, "/barcodes/")
    if(!file.exists(paste0(bc.loc, "/bc_identification/"))) {
      dir.create(path = paste0(bc.loc, "/bc_identification/"), recursive = T)
    }
    cat("running",id)
    cmd = paste0("cat ", files[id], " | ", bc.splitter, " --partial 3 --eolbcfile ", bc.loc, "/bctable_eol.txt --bolbcfile ", bc.loc, "/bctable_bol.txt --mismatches ", bc.missmatch, " --o ", bc.loc, "/bc_identification/", id, ".txt")
    bcident.out=try(system(cmd, intern = TRUE))
  }
  return(files)
}
# formate the output into count matrix
tabulateMappingResults <- function (files, lib.loc, res.loc, eval.cut) {
  results <- lapply (names(files),function(id)
  {
    file <- files[id]
    bc.ident <- read.delim(paste0(res.loc, "/barcodes/bc_identification/", id, ".txt"), stringsAsFactors = F, header=F, na.strings = "unmatched")
    colnames(bc.ident) <- c("bol", "eol")
    insertmatch <- read.table(paste0(res.loc, "/blast/", id, "_res_inserts.txt"),stringsAsFactors=F)
    colnames(insertmatch) <- c("query","subject","percentage.id","alignment.length","mismatches","gap.openings","query.start","query.end","subject.start","subject.end","E.value","bit.score")
    # filter the low quality matches
    insertmatch <- insertmatch[insertmatch$E.value<eval.cut,]
    insertmatch <- insertmatch[order(insertmatch$E.value),]
    # if there are more than one match for a query
    insertmatch <- insertmatch[which(!duplicated(insertmatch$query)),]
    insertmatch <- insertmatch[order(insertmatch$query),]
    # match to the barcodes given the query numbers are corresponding to the order of the barcodes
    insertmatch$bol <- bc.ident$bol[insertmatch$query] # query hold the oder or reads in for blast and barcode identification  
    insertmatch$eol <- bc.ident$eol[insertmatch$query] # and therefore can be used to map the two tables
    res <- data.frame(amplicon=insertmatch$subject,left.bc=pmin(insertmatch$eol,insertmatch$bol),right.bc=pmax(insertmatch$eol,insertmatch$bol)) # left barcodes minimun in eol and right barcodes max in eol
    # blast file has duplicates of reads and we 
    rownames(res) <- insertmatch$query 
    res <- res[grepl("L[[:digit:]]+",res$left.bc) & grepl("R[[:digit:]]+",res$right.bc),] # removing the nas
   # list(fractions=c(frac.compl=nrow(res)/nrow(bc.ident),frac.insert=nrow(insertmatch) / nrow(bc.ident),frac.bol=sum(!is.na(bc.ident$bol))/nrow(bc.ident),frac.eol=sum(!is.na(bc.ident$eol))/nrow(bc.ident)),result.counts=table(res),result=res)
  })
  names(results) <- names(files)
  return(results)
}