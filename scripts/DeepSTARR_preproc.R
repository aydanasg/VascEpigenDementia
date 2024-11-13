library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)

#Loading the peak summits from a bed file and creating a 250bp region around them

summits <- import(summitbed)
start(summits) <- start(summits) - 125
end(summits) <- end(summits) + 125

#Loading the bam files
bam <- import(bam)
#Counting overlaps between the reads and the 250bp regions
overlaps <- countOverlaps(summits, bam)

seq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38,summits)
seq <- as.list(as.character(seq))
chr <- as.vector(seqnames(summits))
df <- cbind(chr, unlist(seq), overlaps)
df <- data.frame(df)
colnames(df) <- c("chr", "seq", "score")

#Save to csv file to use in model training
write.csv(df, outputPath, row.names = FALSE)