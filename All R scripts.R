## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("EnsDb.Hsapiens.v75")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("ensembldb")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("tximport")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("rjson")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("readr")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("statmod")

library(tximport)
dir <-'~/Desktop/quants'
list.files(dir)

samples <- read.table (file.path (dir, "samples.txt"), header = TRUE)
samples

files <- file.path (dir, c("quant_normal.sf", "quant_tumor1.sf", "quant_tumor2.sf", "quant_tumor3.sf" ))
names(files) <- paste0("sample", 1:4)
all(file.exists(files))

#create the tx2gene
library(EnsDb.Hsapiens.v75)
library(Rsamtools)
edb <- EnsDb.Hsapiens.v75
txdf <- transcripts(EnsDb.Hsapiens.v75, return.type="DataFrame")
tx2gene <- as.data.frame(txdf[,c("tx_id","gene_id")])

tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))
head(tx2gene)

#create txi
library(tximport)
library(readr)
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
names(txi)
head(txi$counts)

txi.tx <- tximport(files, type= "salmon", txOut=TRUE, tx2gene=tx2gene)

txi.sum <- summarizeToGene(txi.tx, tx2gene)
all.equal(txi$counts, txi.sum$counts)

files <- file.path (dir, c("quant_normal.sf", "quant_tumor1.sf", "quant_tumor2.sf", "quant_tumor3.sf" ))
names(files) <- paste0("sample", 1:4)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi.salmon$counts)

#Creating a DGEList object
library(edgeR)
cts <- txi$counts
normMat <- txi$length
normMat <- normMat/exp(rowMeans(log(normMat)))
library(edgeR)
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y$offset <- t(t(log(normMat)) + o)
colnames(cts) <- c("Normal", "Tumor1", "Tumor2", "Tumor3")
dim(cts)
head(cts)
y$samples

keep <- rowSums(cpm(y) > 2) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
y$offset <- t(t(log(normMat)) + o)

#filtering
countsPerMillion <- cpm(y)
summary(countsPerMillion)
#'summary' is a useful function of exploring numeric data.

countCheck <- countsPerMillion > 1
head(countCheck)

keep <- which(rowSums(countCheck) >= 2)
dgList <- y[keep,]
summary(cpm(y))

keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]

#Normalization
y <- calcNormFactors(y, method="TMM")

pseudoCounts <- log2(y$counts + 1)
head(pseudoCounts)

normCounts <- cpm(y)
pseudoNormCounts <- cpm(y, log = TRUE, prior.count = 1)
boxplot(pseudoNormCounts, col = "gray", las = 3, cex.names = 1)

#Data Exploration(Multidimensional scaling)
plotMDS(y)
design <- model.matrix(~ 0 + group)
design

#Design the matrix
Group <- factor(c(1,2,2,2))
Tissue <- factor(c("N","T","T","T"))
data.frame(Sample=colnames(y),Group,Tissue)

design <- model.matrix(~Group)
rownames(design) <- colnames(y)
design

#Estimate dispersions
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y$offset <- t(t(log(normMat)) + o)
colnames(cts) <- c("Normal", "Tumor1", "Tumor2", "Tumor3")
dim(cts)
head(cts)

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

y <- estimateGLMTrendedDisp(y, design, robust=TRUE)
y$trended.dispersion

y <- estimateGLMTagwiseDisp(y, design, robust=TRUE)
y$tagwise.dispersion

#Biological coefficient of variation
plotBCV(y)

#differential expression
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
lrt <- glmLRT(fit)
topTags(lrt)

deGenes <- decideTestsDGE(lrt, p=0.001)
deGenes <- rownames(lrt)[as.logical(deGenes)]
plotSmear(lrt, de.tags=deGenes)
abline(h=c(-1, 1), col="2")

