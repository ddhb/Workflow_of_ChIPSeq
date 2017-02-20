####################################################################
## Please install the latest release of 'Microsoft R Open'. 
## https://mran.microsoft.com/open/
## Then, entering the commands.
####################################################################

####################################################################
## Step 0 | Install the packages
####################################################################
## Step 0-1. Install bioconductor base program.
source("https://bioconductor.org/biocLite.R")
biocLite()

## Step 0-2. Install each package.
biocLite("BiocParallel")
biocLite("parallel")
biocLite("ShortRead")
biocLite("RSQLite")
biocLite("QuasR")
biocLite("BSgenome")
biocLite("Rsamtools")
biocLite("rtracklayer")
biocLite("GenomicFeatures")
biocLite("Hmisc")
biocLite("Gviz")
biocLite("XML")
biocLite("mosaics")
biocLite("ChIPseeker")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
biocLite("clusterProfiler")
biocLite("ReactomePA")

####################################################################
## Step 1 | Set the working directory and loading packages
####################################################################
## Step 1-1. Set working directory if the path is 'K:/LAB4'. 
## Unless please correct the working directory path.
working_dir="K:/LAB4"
setwd(working_dir)

## Step 1-2. Library loading.
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

####################################################################
## (Optional) Step 2 | Downloading sample data from EBI (2.8GB and 1.5GB respectively) ## 30 minutes
####################################################################
## Step 2-1. Assignment variables for downloading zipped fastq files.
download_filename_Control = "Sample_data_control.fastq.gz"
download_filename_Case = "Sample_data_case.fastq.gz"
sample_control = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR227/SRR227391/SRR227391.fastq.gz"
sample_chip = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR227/SRR227441/SRR227441.fastq.gz"

## Step 2-2. Start download file with internal method from EBI.
download.file(sample_control, destfile=download_filename_Control, method="internal")
download.file(sample_chip, destfile=download_filename_Case, method="internal")

####################################################################
## Step 3 | Align to reference genome # 1.2 hours using 7 cores
##        | Recommend to start at this point
####################################################################
## Step 3-1. Assignment variables for alignment to hg19 reference genome
filename_Control = strsplit(download_filename_Control, ".fastq.gz")[[1]]
filename_Case = strsplit(download_filename_Case, ".fastq.gz")[[1]]
controlFile <- "sample_control.txt"
genomeFile <- "BSgenome.Hsapiens.UCSC.hg19"
caseFile <- "sample_Case.txt"
genomeFile <- "BSgenome.Hsapiens.UCSC.hg19"

## Step 3-2. Making a matrix file sample and control each other to use qAlign function. Please refer the "sample_*.txt" file in GitHub
Matrix_file_control <- matrix(c("FileName",download_filename_Control,"SampleName","Sample1"),nrow=2,ncol=2)
write.table(Matrix_file_control, file="sample_control.txt", sep="\t",row.names=FALSE, col.names=FALSE)
Matrix_file_case <- matrix(c("FileName",download_filename_Case,"SampleName","Sample1"),nrow=2,ncol=2)
write.table(Matrix_file_case, file="sample_Case.txt", sep="\t",row.names=FALSE, col.names=FALSE)

## Step 3-3. Assigning core and aligning to hg19 reference genome
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
Aligned_control<-qAlign(sampleFile=controlFile, genome=genomeFile,
                        clObj=cl, alignmentsDir = working_dir)
Aligned_case<-qAlign(sampleFile=caseFile, genome=genomeFile,
                     clObj=cl, alignmentsDir = working_dir)

## Step 3-4. Post-processing (releasing core and removing temporal matrix files)
stopCluster(cl)
file.remove(controlFile)
file.remove(caseFile)

####################################################################
## Step 4 | Quality check ## 1.2 hours using 7 cores
####################################################################
## Step 4-1. The QC plotting of case and control files
qQCReport(Aligned_control, pdfFilename = paste(filename_Control,"_QCReport.pdf",sep=""))
qQCReport(Aligned_case, pdfFilename = paste(filename_Case,"_QCReport.pdf",sep=""))

## Step 4-2. Renaming aligned output files
file.rename(list.files(pattern=paste("^",filename_Control,"(.*).bam$",sep="")), to=paste(filename_Control,".bam", sep=""))
file.rename(list.files(pattern=paste("^",filename_Control,"(.*).bam.bai$",sep="")), to=paste(filename_Control,".bam.bai", sep=""))
file.rename(list.files(pattern=paste("^",filename_Control,"(.*).bam.txt$",sep="")), to=paste(filename_Control,".bam.txt", sep=""))
file.rename(list.files(pattern=paste("^",filename_Case,"(.*).bam$",sep="")), to=paste(filename_Case,".bam", sep=""))
file.rename(list.files(pattern=paste("^",filename_Case,"(.*).bam.bai$",sep="")), to=paste(filename_Case,".bam.bai", sep=""))
file.rename(list.files(pattern=paste("^",filename_Case,"(.*).bam.txt$",sep="")), to=paste(filename_Case,".bam.txt", sep=""))

####################################################################
## Step 5 | Peak calling ## 13 minutes
####################################################################
## Step 5-1. Preprocessing before calling peaks
filename_Control_bam = paste(filename_Control,".bam",sep="")
filename_Case_bam = paste(filename_Case,".bam",sep="")

## Step 5-2. Setting the fragment length and binsize and contructing bin-level files from bam files
fragment_length = 200
binSize = 200
constructBins(infile=filename_Control_bam, fileFormat="bam", outfileLoc="./",byChr=FALSE, useChrfile=FALSE,
              chrfile=NULL, excludeChr=NULL, PET=FALSE, fragLen=fragment_length,binSize=binSize,capping=0)
constructBins(infile=filename_Case_bam, fileFormat="bam", outfileLoc="./",byChr=FALSE, useChrfile=FALSE,
              chrfile=NULL, excludeChr=NULL, PET=FALSE, fragLen=fragment_length,binSize=binSize,capping=0)
fragment_fControlb = paste(filename_Control_bam,"_fragL",fragment_length,"_bin",binSize,".txt",sep="")
fragment_fCaseb = paste(filename_Case_bam,"_fragL",fragment_length,"_bin",binSize,".txt",sep="")

## Step 5-3. Reading bin-level data
fileName <- c(fragment_fControlb,fragment_fCaseb)
binTFBS <- readBins(type=c("input","chip"), fileName=fileName)

## Step 5-4. Fitting the MOSAiCS model
fitTFBS <- mosaicsFit(binTFBS,analysisType="IO",bgEst="rMOM")

## Step 5-5. Identifying peaks based on the fitted model
peakTFBS <- mosaicsPeak (fitTFBS, signalModel="2S", FDR=0.05, maxgap=200, minsize=50, thres=10)

## Step 5-6. Exporting result to 'bed' format with the process of elimination of header line
result_peakfile = "sharp_peakTFBS.bed"
export (peakTFBS, type="bed", filename=result_peakfile)
peakfile <- read.table(result_peakfile, header=FALSE)[-1,]
header.eliminated.file = "h.e_peakTFBS.bed"
write.table(peakfile, file=header.eliminated.file, sep="\t", col.names=FALSE, row.names=FALSE)

####################################################################
## Step 6 | Annotation and Visualization Peak calling ## 2 minutes
####################################################################
## Step 6-1. Preprocessing
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
interesing_chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10"
         ,"chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20"
         ,"chr21","chr22","chrM","chrX","chrY")

## Step 6-2. Reading peak file and coverage plot
peak <- readPeakFile(header.eliminated.file)
covplot(peak, weightCol="V5", chrs=interesing_chrs)

## Step 6-3. Get promoter and heatmap of ChIP binding to TSS regions
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), 
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

## Step 6-4. Annotate peak and illustrate pie-plot & upset-plot
peakAnno <- annotatePeak (header.eliminated.filename, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno)
upsetplot(peakAnno, vennpie=TRUE)

## Step 6-5. Functional enrichment analysis and plotting
gene <- seq2gene(peak, tssRegion=c(-1000,1000), flankDistance=3000, TxDb=txdb)
pathway <- enrichPathway(gene)
dotplot(pathway)
