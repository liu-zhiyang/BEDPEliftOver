################################################################################
#########convert genome coordinates in bedpe format between assemblies##########
################################################################################

#set work directory where chain file and file to convert locate in
setwd("./")
library(dplyr)
library(rtracklayer)

#read sample data
#here geome coordinates is in hg19 genome, and will be converted to hg38
hg19.bedpe <- read.table("EnhancerPredictions.bedpe")
colnames(hg19.bedpe) <- c("chr_1","start_1","end_1",
                          "chr_2","start_2","end_2",
                          "name","score","strand_1","strand_2")
hg19.1 <- GRanges(seqnames=hg19.bedpe$chr_1,
                  ranges=IRanges(start=hg19.bedpe$start_1,
                                 end=hg19.bedpe$end_1))
hg19.2 <- GRanges(seqnames=hg19.bedpe$chr_2,
                  ranges=IRanges(start=hg19.bedpe$start_2,
                                 end=hg19.bedpe$end_2))
                                     
#import chain file
#hg19 to hg38 chain file can be downloaded from 
#https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
#more chain file can be found from ucsc genome browser
chain.file <- import.chain.file("hg19ToHg38.over.chain.file")

#liftOver each part of bedpe file to hg38 genome
hg38.1 <- liftOver(hg19.1,chain.file)
hg38.2 <- liftOver(hg19.2,chain.file)

#check and discard inconsistent coordinates
ex <- vector()
for (i in 1:length(hg38.1)) {
    if(length(hg38.1[[i]]) != 1){
        ex <- c(ex,i)
    }
}
hg38.1 <- liftOver(hg19.1[-ex],chain.file) %>% unlist()
hg38.2 <- liftOver(hg19.2[-ex],chain.file) %>% unlist()

#build and write the successfully coverted bedpe file
hg38.bedpe <- cbind(data.frame(hg38.1)[,1:3],
                    data.frame(hg38.2)[,1:3],
                    hg19.bedpe[-ex,7:10])
write.table(hg38.bedpe,"hg38liftOver.EnhancerPredictions.bedpe",
            col.names=FALSE, row.names=FALSE, sep="\t", quote=F)
            
#write coordinates that failed to convert
write.table(hg19.bedpe,"unliftOver.EnhancerPredictions.bedpe",
            col.names=FALSE, row.names=FALSE, sep="\t", quote=F)
