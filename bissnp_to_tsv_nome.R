# supply bedGraph filename
stopifnot(length(commandArgs(TRUE))==2)
samp <- commandArgs(TRUE)[1]
filename <- paste0(samp, ".bissnp.", samp, ".meth.bed")
stopifnot(file.exists(filename))

library(GenomicRanges)
library(data.table)

# Find the appropriate BSgenome package
build.name <- commandArgs(TRUE)[2]
packages <- installed.packages()[,"Package"]
bs.build <- grep(paste0("^BSgenome.*UCSC.", build.name), packages)
stopifnot(length(bs.build)==1)
library(packages[bs.build], character.only=TRUE)
# get object of the build
build <- get(sub("BSgenome.", "", gsub(".UCSC.*", "", packages[bs.build])))

overlapSums <- function(x, y, z) {
    stopifnot(class(x)=="GRanges")
    stopifnot(class(y)=="GRanges")
    stopifnot(class(z)=="numeric" | class(z)=="integer")
    ov <- as.matrix(findOverlaps(y, x))
    ovSums <- rep(as.numeric(0), length(y))
    ovSums[unique(ov[,1])] <- viewSums(Views(z[ov[,2]], ranges(Rle(ov[,1]))))
    ovSums
}

#GCH = GCs methylated by M.CviPI, not GCG as confounded by CpG methylation
GCH <- c(resize(vmatchPattern("GCA", build), 1, fix="center"), 
         resize(vmatchPattern("GCC", build), 1, fix="center"), 
         resize(vmatchPattern("GCT", build), 1, fix="center"))
# Interleave + and - strand
GCH <- GCH[order(as.factor(seqnames(GCH)), start(GCH), as.factor(strand(GCH)))]
GCH$cov <- GCH$C <- 0L

#WCG = not HCG cos CCG is methylated by M.CviPI
WCG <- c(resize(vmatchPattern("ACG", build), 1, fix="center"), 
         resize(vmatchPattern("TCG", build), 1, fix="center"))
# Interleave + and - strand
WCG <- WCG[order(as.factor(seqnames(WCG)), start(WCG), as.factor(strand(WCG)))]
WCG$cov <- WCG$C <- 0L

# Read bissnp table
tab <- fread(filename)
tab <- GRanges(tab[[1]], IRanges(tab[[2]]+1, width=1), C=tab[[5]], T=tab[[6]], cov=tab[[5]]+tab[[6]], context=tab[[7]])

# Match up GCHs/WCGs and bissnp table
GCH.ov <- as.matrix(findOverlaps(GCH, tab))
WCG.ov <- as.matrix(findOverlaps(WCG, tab))

# Sanity checks
stopifnot(all(!duplicated(GCH.ov[,1])))
stopifnot(all(!duplicated(WCG.ov[,1])))
stopifnot(all(!duplicated(c(GCH.ov[,2], WCG.ov[,2]))))

# Assign bissnp covered data to GCH/WCGs
GCH$C[GCH.ov[,1]] <- tab$C[GCH.ov[,2]]
GCH$cov[GCH.ov[,1]] <- tab$cov[GCH.ov[,2]]
WCG$C[WCG.ov[,1]] <- tab$C[WCG.ov[,2]]
WCG$cov[WCG.ov[,1]] <- tab$cov[WCG.ov[,2]]

# Sanity check
stopifnot(all(GCH$cov>=GCH$C))
stopifnot(all(WCG$cov>=WCG$C))

# Export stranded GCH and WCG calls
tmp <- data.frame("chr"=as.character(seqnames(GCH)), "position"=start(GCH), "strand"=as.character(strand(GCH)), GCH$C, GCH$cov, stringsAsFactors=FALSE)
names(tmp)[4:5] <- paste(samp, c("C", "cov"), sep=".")
write.table(tmp, file=gzfile(paste0(samp, ".GCH.strand.tsv.gz")), quote=FALSE, sep="\t", row.names=FALSE)

tmp <- data.frame("chr"=as.character(seqnames(WCG)), "position"=start(WCG), "strand"=as.character(strand(WCG)), WCG$C, WCG$cov, stringsAsFactors=FALSE)
names(tmp)[4:5] <- paste(samp, c("C", "cov"), sep=".")
write.table(tmp, file=gzfile(paste0(samp, ".CpGs.strand.tsv.gz")), quote=FALSE, sep="\t", row.names=FALSE)

# Pool strands for GCH - -ve strand is first so shift +ve strand -1bp
GCH.pool <- GCH
values(GCH.pool) <- NULL
GCH.pool[strand(GCH.pool)=="+"] <- shift(GCH.pool[strand(GCH.pool)=="+"], -1 )
strand(GCH.pool) <- "*"
GCH.pool <- reduce(GCH.pool)
GCH.pool$C <- overlapSums(GCH, resize(GCH.pool, 2, fix="start"), GCH$C)
GCH.pool$cov <- overlapSums(GCH, resize(GCH.pool, 2, fix="start"), GCH$cov)

# Pool strands for WCG - +ve strand is first so shift -ve strand -1bp
WCG.pool <- WCG
values(WCG.pool) <- NULL
WCG.pool[strand(WCG.pool)=="-"] <- shift(WCG.pool[strand(WCG.pool)=="-"], -1)
strand(WCG.pool) <- "*"
WCG.pool <- reduce(WCG.pool)
WCG.pool$C <- overlapSums(WCG, resize(WCG.pool, 2, fix="start"), WCG$C)
WCG.pool$cov <- overlapSums(WCG, resize(WCG.pool, 2, fix="start"), WCG$cov)

# Sanity check
stopifnot(all(GCH.pool$cov>=GCH.pool$C))
stopifnot(all(WCG.pool$cov>=WCG.pool$C))

# Export pooled data
GCH <- GCH.pool
tmp <- data.frame("chr"=as.character(seqnames(GCH)), "position"=start(GCH), GCH$C, GCH$cov, stringsAsFactors=FALSE)
names(tmp)[3:4] <- paste(samp, c("C", "cov"), sep=".")
write.table(tmp, file=gzfile(paste0(samp, ".GCH.tsv.gz")), quote=FALSE, sep="\t", row.names=FALSE)

WCG <- WCG.pool
tmp <- data.frame("chr"=as.character(seqnames(WCG)), "position"=start(WCG), WCG$C, WCG$cov, stringsAsFactors=FALSE)
names(tmp)[3:4] <- paste(samp, c("C", "cov"), sep=".")
write.table(tmp, file=gzfile(paste0(samp, ".CpGs.tsv.gz")), quote=FALSE, sep="\t", row.names=FALSE)
