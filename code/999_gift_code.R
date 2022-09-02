library("RJSONIO")
library("modules")
library(dada2); packageVersion("dada2") ## Dada2 v.1.10.1
library(ShortRead); packageVersion("ShortRead") ## ShortRead 1.40.0
library(phyloseq); packageVersion("phyloseq") ## phyloseq 1.26.0
library(dplyr)

setwd("~/BaseSpace/IrsiCaixa-Run16S-360489476/FASTQ_Generation_2022-08-14_20_40_38Z-596888292/")
fnFs <- sort(list.files( pattern="_R1_001.fastq.gz|_R1.fastq.gz"))
fnRs <- sort(list.files( pattern="_R2_001.fastq.gz|_R2.fastq.gz"))
print("Forward Fiels are")
print(fnFs[1:5])
print("Reverse Files are:")
print(fnRs[1:5])

sample.files.name <- sapply(strsplit(fnFs, "_R1"), `[`, 1)
sample.files.name <- sapply(strsplit(sample.files.name, "_S"), `[`, 1)
print("File names from fastq files are:")
print(sample.files.name[1:5])
filter_metadata<-TRUE

fnFs <- file.path( fnFs)
print(fnFs[1:5])
fnRs <- file.path( fnRs)
print(fnRs[1:5])

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
filt_path <- file.path(".", "filtered")
filtFs <- file.path(filt_path, paste0(sample.files.name, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.files.name, "_R_filt.fastq.gz"))
print(filtFs[1:5])
print(filtRs[1:5])
system(paste0("gunzip -c ",fnFs[1]," | head -n 4"))
system(paste0("gunzip -c ",fnRs[1]," | head -n 4"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,260),
                     maxN=0, maxEE=c(4,10), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=T,verbose=T)

errF <- learnErrors(filtFs, multithread=TRUE,nreads=500000,randomize=T,verbose=T)
errR <- learnErrors(filtRs, multithread=TRUE,nreads=500000,randomize=T,verbose=T)

derepFs <- derepFastq(filtFs, verbose=TRUE,n=1000000)
derepRs <- derepFastq(filtRs, verbose=TRUE,n=1000000)

names(derepFs) <- sample.files.name
names(derepRs) <- sample.files.name

dadaFs <- dada(derepFs, err=errF, multithread=T)
dadaRs <- dada(derepRs, err=errR, multithread=T)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE,minOverlap = 6)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.files.name
track<-data.frame(track)
track$SampleID<-rownames(track)
head(track)
summary(track)


track %>%
  dplyr::summarise(input_pct=100*(input/input),filtered=100*(filtered/input),denoised=100*(denoised/input),
                   merged=100*(merged/input),nonchim=100*(nonchim/input)) %>%
  dplyr::rename(input=input_pct) %>%
  tidyr::pivot_longer(cols = c("input", "filtered", "denoised", "merged", "nonchim")) %>%
  dplyr::mutate(name=forcats::fct_relevel(name,c("input", "filtered", "denoised", "merged","nonchim"))) %>%
  ggstatsplot::ggwithinstats(.,name,value,pairwise.display ="none")

