library(htmltools)
library(arrayQualityMetrics)
#library(DESeq)
library(DESeq2)
library(tidyr)
library(genefilter)
library(cowplot)
library(readr)
library(gplots)
library(reshape2)
library(plyr)
library(dplyr)
library(apeglm)
library(IHW)
library(pheatmap)
library(vsn)
library(ggplot2)
library(ggpubr)
library(ggtext)
library(ape)
library(vegan)
library(RColorBrewer)
library(stringr)
library(KOGMWU)
library(DESeq2)
library(arrayQualityMetrics)
library(Biobase)


#Importing counts data and preparing for analyses
setwd("C:/Users/anaca/OneDrive/Escritorio/VS Nuevo O.arbuscula simbiosis/My O.arbuscula")
coldata<-read.table("data/Samp_data.txt", header=TRUE)
rownames(coldata)<-coldata$Samp
coldata<-coldata[,2:3]

OCAA<-read.table("data/OCAA_coral_sam.counts")
OCBA<-read.table("data/OCBA_coral_sam.counts")
OCCA<-read.table("data/OCCA_coral_sam.counts")
OCDA<-read.table("data/OCDA_coral_sam.counts")
OCFA<-read.table("data/OCFA_coral_sam.counts")

OCAS<-read.table("data/OCAS_coral_sam.counts")
OCBS<-read.table("data/OCBS_coral_sam.counts")
OCCS<-read.table("data/OCCS_coral_sam.counts")
OCDS<-read.table("data/OCDS_coral_sam.counts")
OCFS<-read.table("data/OCFS_coral_sam.counts")

cts<-full_join(OCAA,OCBA, by="V1") # full_join merges all the data even if there's a gene found in A and that isn't in B and vice versa, missing data are given NA) 
cts<-full_join(cts,OCCA, by="V1")
cts<-full_join(cts,OCDA, by="V1")
cts<-full_join(cts,OCFA, by="V1")
cts<-full_join(cts,OCAS, by="V1")
cts<-full_join(cts,OCBS, by="V1")
cts<-full_join(cts,OCCS, by="V1")
cts<-full_join(cts,OCDS, by="V1")
cts<-full_join(cts,OCFS, by="V1")

#Fix Column headers 
colnames(cts)<-c("contig","OCAA", "OCBA", "OCCA", "OCDA", "OCFA", "OCAS", "OCBS", "OCCS", "OCDS", "OCFS")

#Make NAs zeros
cts[is.na(cts)]<-0 


# Los autores llaman al archivo "coral_gene_names.tab" pero no se encuentra en el repositorio de GitHub.
# así que, descargamos el transcriptoma desde https://sites.bu.edu/davieslab/data-code/ 

gg=read.delim("O_arbuscula_isogroup_to_genename.tab",sep="\t", header=FALSE)
colnames(gg)<-c("contig","genename")
# La primera columna es el identificador del contig (“contiguous sequence”), isogroup o transcrito (por ejemplo, “comp100001_c0_seq1”).
# La segunda columna es el nombre del gen anotado (por ejemplo, “glutamine synthetase”)

# removing A. thaliana genes from counts
cts<-left_join(cts,gg, by="contig")
cts<-cts[grep("thaliana", cts$genename, invert=TRUE),]

#Name columns and remove gene and genename column
rownames(cts)<-cts$contig
cts<-cts[,2:11]

# Verificar que todo está funcionando
View(cts)