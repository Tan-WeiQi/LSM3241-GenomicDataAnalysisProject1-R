#LSM3241 CA1 Code

#install Bioconductor
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()

#install Bioconductor packages
BiocManager::install('foo')

#install relevant packages needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery", version = "3.8")

BiocManager::install("affy")
BiocManager::install("limma")
BiocManager::install("hgu133plus2.db")
BiocManager::install("org.Hs.eg.db")

#load packages installed
library(affy)
library(limma)
library(hgu133plus2.db)
library(org.Hs.eg.db)

#check if packages are loaded
sessionInfo()

#Calling GSE50697
library(GEOquery)
#downloading the gse file
gse <- getGEO('GSE50697',GSEMatrix = FALSE)

#untar file
#untar("GSE50697_RAW.tar", files = NULL, list = FALSE, exdir = "GSE50697_RAW/",
      #compressed = NA, extras = NULL, verbose = FALSE,
      #restore_times =  TRUE, tar = Sys.getenv("TAR"))

#retrieving whole GSE
names(GSMList(gse))

#Getting the raw data for the series using GEOquery
filePaths <- getGEOSuppFiles('GSE50697')
#generating a vector containing the name of all the CEL files.
list.celfiles('GSE50697_RAW')

#finding info from meta data of gse
names(Meta(gse))

#creating GSM object names
gsm <- GSMList(gse)[[1]]

#find information we need from the metadata of gsm
names(Meta(gsm))

#extract metadata of gsm not in gse
names(Meta(gsm))[!(names(Meta(gsm)) %in% names(Meta(gse)))]
Meta(gsm)[!(names(Meta(gsm)) %in% names(Meta(gse)))]

#From the metedata information, we decidede to use the following elements
#source_name_ch1 and characteristics_ch1

#getting sample growth conditions into a data frame
culture_medium <- function(gsm) { 
  Meta(gsm)[['characteristics_ch1']][2]
}
sapply(GSMList(gse),culture_medium)
pd <- data.frame(culture=as.factor(sapply(GSMList(gse),culture_medium)))
pd

#simplifying our dataframe, convert columns to 2 values only
pd$culture <- as.factor(pd$culture)
levels(pd$culture) <- c("control","miR203")
#to enable the kable function, we need to install the knitr package
install.packages("knitr")
library(knitr)
kable(pd)

#Reading in the CEL files with the phenoData
celfiles <- paste0('GSE50697_RAW/', list.celfiles('GSE50697_RAW/'),'.')
affydata <- read.affybatch(celfiles,phenoData = new("AnnotatedDataFrame",pd))
phenoData(affydata)

#pseudo images of chips
image(data[,1])
image(data[,2])
image(data[,3])
image(data[,4])
image(data[,5])
image(data[,6])

#CEL file densities before normalisation or background correction
plotDensity.AffyBatch(affydata)

#perform RMA normalisation
eset <- rma(affydata)
plotDensity(exprs(eset),xlab='log intensity',main="feature level densities after RMA",lwd=2)

#phenotype data for each samples after rma is retained
pData(eset)


#generate model matrix
model <- model.matrix( ~ 0 + eset$culture)
#rename the model columns to correspond to the different growth conditions
colnames(model) <- levels(eset$culture)
model

#look at when the growth conditions differ, we create contrast
contrasts <- makeContrasts(control - miR203, levels=model)
contrasts


#fit model and contrast matrix into a data
fit <- lmFit(eset, model)
fit
fitted.contrast <- contrasts.fit(fit,contrasts)
fitted.contrast

#calculate test statistics via performing eBayes correction
fitted.ebayes <- eBayes(fitted.contrast)
fitted.ebayes

#extracting differentially expressed genes
topTable(fitted.ebayes)

#Get a limited number of probesets
ps <- rownames(topTable(fitted.ebayes))
ps

#The AnnotationDbi interface
#look at available columns for our chip
columns(hgu133plus2.db)
#Look at which that can be used as keys
keytypes(hgu133plus2.db)
#realise all can be used as keys, choose "PROBEID" as most suitable
head(keys(hgu133plus2.db,keytype="PROBEID"))
AnnotationDbi::select(hgu133plus2.db,ps,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")


#Restrict to upregulated genes
differentially_expressed <- topTable(fitted.ebayes,number=Inf,p.value=0.1,lfc=1)
differentially_expressed
upregulated <- differentially_expressed[differentially_expressed$logFC > 0,]
genes_of_interest <- AnnotationDbi::select(hgu133plus2.db,
                                           keys=rownames(upregulated),
                                           columns=c("SYMBOL","ENTREZID","GENENAME"),
                                           keytype="PROBEID")
genes_of_interest
#save data into a csv file
write.csv(genes_of_interest,'genes_of_interest.csv')


#volcanoplots
#
volcanoplot(fitted.ebayes)
interesting_genes <- topTable(fitted.ebayes,number=Inf,p.value = 0.1,lfc=1)
volcanoplot(fitted.ebayes, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes)))
points(interesting_genes[['logFC']],-log10(interesting_genes[['P.Value']]),col='red')


#heatmaps
#normalise expression value
eset_of_interest <- eset[rownames(interesting_genes),]
heatmap(exprs(eset_of_interest))


#fix and beautify heatmap
install.packages("RColorBrewer")
library(RColorBrewer)
eset_of_interest <- eset[rownames(interesting_genes),]
heatmap(exprs(eset_of_interest),
        labCol=eset$culture,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))


#Results
recommendations <- AnnotationDbi::select(hgu133plus2.db,
                                           keys=rownames(interesting_genes),
                                           columns=c("SYMBOL","ENTREZID","GENENAME"),
                                           keytype="PROBEID")
write.csv(recommendations,'recommendations.csv')
