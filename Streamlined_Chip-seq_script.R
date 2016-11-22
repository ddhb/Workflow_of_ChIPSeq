# script 01 , Downloading sample data from Sequence Read Archive
library(SRAdb)
if(!file.exists('SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile()
sqlfile <- getSRAdbFile()
file.info('SRAmetadb.sqlite')
sra_con <- dbConnect(SQLite(),sqlfile)
getSRAfile( c("SRR227391","SRR227441"), sra_con, fileType = 'fastq' )

# script 02 , Fastq files trimming
library("ShortRead")
Trimming <- function(reads, minQuality, firstBase=1, minLength=5){
  qualMat <- as(FastqQuality(quality(quality(reads))),'matrix')
  qualList <- split(qualMat,row(qualMat))
  ends <- as.integer(lapply(qualList,
                          function(x){which(x < minQuality)[1]-1}))
  
  starts <- as.integer(lapply(ends,function(x){min(x+1,firstBase)}))
  newQ <- ShortReadQ(sread=subseq(sread(reads),start=starts,end=ends),
                   quality=new(Class=class(quality(reads)),
                               quality=subseq(quality(quality(reads)),
                                              start=starts,end=ends)),
                   id=id(reads))
  
  lengthCutoff <- srFilter(function(x) {
    width(x)>=minLength
  },name="length cutoff")
  newQ[lengthCutoff(newQ)]
} 

reads <- readFastq("SRR227391.fastq.gz") 
trimmedReads <- Trimming(reads=reads, minQuality=5, firstBase=4, minLength=3) 
writeFastq(trimmedReads, file="SRR227391_trimmed.fastq.gz")

reads <- readFastq("SRR227441.fastq.gz") 
trimmedReads <- Trimming(reads=reads, minQuality=5, firstBase=4, minLength=3) 
writeFastq(trimmedReads, file="SRR227441_trimmed.fastq.gz")

# script 03 , Alignment to reference genome 
library(QuasR)
library(BSgenome)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)
library(Gviz)

genomeName <- "BSgenome.Hsapiens.UCSC.hg19"
Aligned_control <- qAlign("SRR227391_trimmed.fastq.gz", genome=genomeName)
files_control <- Sys.glob("SRR227391*.bam")
file.rename(from=files_control[[1]], "SRR227391_trimmed.bam")

Aligned_H3K4me3 <- qAlign("SRR227441_trimmed.fastq.gz", genome=genomeName)
files_H3K4 <- Sys.glob("SRR227441*.bam")
file.rename(from=files_H3K4[[1]], "SRR227441_trimmed.bam")

# script 04 , Qualityt check of aligned files
qQCReport(Aligned_control, "./QCreport_control.pdf")
qQCReport(Aligned_H3K4me3, "./QCreport_H3K4me3.pdf")

# script 05 , Constructing bin-level files from aligned read file
library("mosaics")
contructBins(infile= "SRR227391_trimmed.bam", 
              fileFormat="bam", outfileLoc="./", bychr=FALSE
              useChrfile=FALSE, chrfile=NULL, excludeChr=NULL,
              PET=FALSE, fragLen=200, binsize=200, capping=0)
contructBins(infile= "SRR227441_trimmed.bam", 
              fileFormat="bam", outfileLoc="./", bychr=FALSE
              useChrfile=FALSE, chrfile=NULL, excludeChr=NULL,
              PET=FALSE, fragLen=200, binsize=200, capping=0)

# script 06 , Reading bin-level data into the R enrivonment
fileName <- c("SRR227391_trimmed.bam_fragL200_bin200.txt",
               "SRR227441_trimmed.bam_fragL200_bin200.txt")
binTFBS <- readBins(type=c("input","chip"), fileName=fileName)

# script 07 , Fitting the MOSAiCS model
fitTFBS <- mosaicsFit (binTFBS, analysisType="IO", bgEst="rMOM")

# script 08 , Identifying peaks based on the fitted model
peakTFBS <- mosaicsPeak (fitTFBS, signalModel="2S", FDR=0.05, maxgap=200, minsize=50, thres=10)

# script 09 , Export the result
export (peakTFBS, type="bed", filename="sharp_peakTFBS.bed")

# script 10 , Reading bed file into ChIP-seeker
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)
peak <- readPeakFiles("sharp_peakTFBS.bed")

# script 11 , ChIP peaks coverage plot
covplot(peak, weightCol="V5")

# script 12 , Heatmap of ChIP binding to TSS regions
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")

# script 13 , Peak annotation
peakAnno <- annotatePeak ("sharp_peakTFBS.bed", tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

# script 14 , Functional enrichment analysis
library(ReactomePA)
gene <- seq2gene(peak, tssRegion=c(-1000,1000), flankDistance=3000, TxDb=txdb)
pathway <- enrichPathway(gene)
dotplot(pathway)





