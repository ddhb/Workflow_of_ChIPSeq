##################################################################
## Step 0 | Installment and set the working directory
####################################################################
working_dir="K:/LAB4"
setwd(working_dir)
library("BiocParallel")
library("parallel")
library("ShortRead")
library("RSQLite")
library("QuasR")
library("BSgenome")
library("Rsamtools")
library("rtracklayer")
library("GenomicFeatures")
library("Hmisc")
library("Gviz")
library("XML")
library("mosaics")
library("ChIPseeker")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("clusterProfiler")
library("ReactomePA")

download_filename_Control = "AGS_Input.fastq.gz"
download_filename_Case = "AGS_H3K4me3.fastq.gz"

####################################################################
## Step 1 | Downloading sample data from EBI (2.8GB and 1.5GB respectively) ## 30 minutes
####################################################################
sample_control = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR227/SRR227391/SRR227391.fastq.gz"
sample_chip = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR227/SRR227441/SRR227441.fastq.gz"

download.file(sample_control, destfile=download_filename_Control, method="internal")
download.file(sample_chip, destfile=download_filename_Case, method="internal")
####################################################################
## Step 2 | Align to reference genome # 1.2 hours using 7 cores
##        | Recommend to start at this point
####################################################################
## Step 2-1. Preprocessing
filename_Control = strsplit(download_filename_Control, ".fastq.gz")[[1]]
filename_Case = strsplit(download_filename_Case, ".fastq.gz")[[1]]

## Step 2-2. Making a matrix file to use qAlign function. Please refer the "sample_*.txt" file in GitHub
Matrix_file_control <- matrix(c("FileName",download_filename_Control,"SampleName","Sample1"),nrow=2,ncol=2)
write.table(Matrix_file_control, file="sample_control.txt", sep="\t",row.names=FALSE, col.names=FALSE)
controlFile <- "sample_control.txt"
genomeFile <- "BSgenome.Hsapiens.UCSC.hg19"

Matrix_file_case <- matrix(c("FileName",download_filename_Case,"SampleName","Sample1"),nrow=2,ncol=2)
write.table(Matrix_file_case, file="sample_Case.txt", sep="\t",row.names=FALSE, col.names=FALSE)
caseFile <- "sample_Case.txt"
genomeFile <- "BSgenome.Hsapiens.UCSC.hg19"

## Step 2-3. Assigning core and Aligning to reference genome
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
Aligned_control<-qAlign(sampleFile=controlFile, genome=genomeFile,
                        clObj=cl, alignmentsDir = working_dir)
Aligned_case<-qAlign(sampleFile=caseFile, genome=genomeFile,
                     clObj=cl, alignmentsDir = working_dir)
stopCluster(cl)
file.remove(controlFile)
file.remove(caseFile)

####################################################################
## Step 3 | Quality check ## 1.2 hours using 7 cores
####################################################################
## Step 3-1. The QC plotting of case and control files
qQCReport(Aligned_control, pdfFilename = paste(filename_Control,"_QCReport.pdf",sep=""))
qQCReport(Aligned_case, pdfFilename = paste(filename_Case,"_QCReport.pdf",sep=""))

## Step 3-2. Renaming aligned output files
file.rename(list.files(pattern=paste("^",filename_Control,"(.*).bam$",sep="")), to=paste(filename_Control,".bam", sep=""))
file.rename(list.files(pattern=paste("^",filename_Control,"(.*).bam.bai$",sep="")), to=paste(filename_Control,".bam.bai", sep=""))
file.rename(list.files(pattern=paste("^",filename_Control,"(.*).bam.txt$",sep="")), to=paste(filename_Control,".bam.txt", sep=""))
file.rename(list.files(pattern=paste("^",filename_Case,"(.*).bam$",sep="")), to=paste(filename_Case,".bam", sep=""))
file.rename(list.files(pattern=paste("^",filename_Case,"(.*).bam.bai$",sep="")), to=paste(filename_Case,".bam.bai", sep=""))
file.rename(list.files(pattern=paste("^",filename_Case,"(.*).bam.txt$",sep="")), to=paste(filename_Case,".bam.txt", sep=""))

####################################################################
## Step 4 | Peak calling ## 13 minutes
####################################################################
## 4-1. Preprocessing
filename_Control_bam = paste(filename_Control,".bam",sep="")
filename_Case_bam = paste(filename_Case,".bam",sep="")

## 4-2. Setting the fragment length and binsize and contructing bin-level files from bam files
fragment_length = 200
binSize = 200
constructBins(infile=filename_Control_bam, fileFormat="bam", outfileLoc="./",byChr=FALSE, useChrfile=FALSE,
              chrfile=NULL, excludeChr=NULL, PET=FALSE, fragLen=fragment_length,binSize=binSize,capping=0)
constructBins(infile=filename_Case_bam, fileFormat="bam", outfileLoc="./",byChr=FALSE, useChrfile=FALSE,
              chrfile=NULL, excludeChr=NULL, PET=FALSE, fragLen=fragment_length,binSize=binSize,capping=0)
fragment_fControlb = paste(filename_Control_bam,"_fragL",fragment_length,"_bin",binSize,".txt",sep="")
fragment_fCaseb = paste(filename_Case_bam,"_fragL",fragment_length,"_bin",binSize,".txt",sep="")

## 4-3. Reading bin-level data
fileName <- c(fragment_fControlb,fragment_fCaseb)
binTFBS <- readBins(type=c("input","chip"), fileName=fileName)

## 4-4. Fitting the MOSAiCS model
fitTFBS <- mosaicsFit(binTFBS,analysisType="IO",bgEst="rMOM")

## 4-5. Identifying peaks based on the fitted model
peakTFBS <- mosaicsPeak (fitTFBS, signalModel="2S", FDR=0.05, maxgap=200, minsize=50, thres=10)

## 4-6. Exporting result to 'bed' format with the process of elimination of header line
result_peakfile = "sharp_peakTFBS.bed"
export (peakTFBS, type="bed", filename=result_peakfile)
peakfile <- read.table(result_peakfile, header=FALSE)[-1,]
header.eliminated.file = "h.e_peakTFBS.bed"
write.table(peakfile, file=header.eliminated.file, sep="\t", col.names=FALSE, row.names=FALSE)
####################################################################
## Step 5 | Annotation and Visualization Peak calling ## 2 minutes
####################################################################
## 5-1. Preprocessing
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
interesing_chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10"
         ,"chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20"
         ,"chr21","chr22","chrM","chrX","chrY")

## 5-2. Reading peak file and coverage plot
peak <- readPeakFile(header.eliminated.file)
covplot(peak, weightCol="V5", chrs=interesing_chrs)

## 5-3. Get promoter and heatmap of ChIP binding to TSS regions
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)


tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")

par(mfrow = c(3, 5))
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), 
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plot(1:10)

edit(plotAvgProf)

par(mfrow = c(3, 2))
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), facet="column",
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()

layout(matrix(c(1,2,3, 4), 2, 2, byrow = TRUE))
edit(plotAvgProf)

peakAnno <- annotatePeak (header.eliminated.filename, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno)
upsetplot(peakAnno, vennpie=TRUE)

gene <- seq2gene(peak, tssRegion=c(-1000,1000), flankDistance=3000, TxDb=txdb)
pathway <- enrichPathway(gene)
dotplot(pathway)
