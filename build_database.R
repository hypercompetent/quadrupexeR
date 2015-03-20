library(rtracklayer)
library(dplyr)

## Read the refGene table for hg19
ref <- read.table("raw_data/refGene.txt",header=F,stringsAsFactors=F)
names(ref) <- c("bin","id","chr","strand","txn.start","txn.end","cds.start","cds.end",
                "n.exons","exon.starts","exon.ends","score","symbol","x","xx","exon.frame")

# Remove non-standard chromosomes
ref <- ref %>% filter(!grepl("_",ref$chr))

# Select TSS locations for plus and minus-strand genes
tss.p <- ref %>% filter(strand == "+") %>% select(chr,txn.start,id,score,strand)
tss.m <- ref %>% filter(strand == "-") %>% select(chr,txn.end,id,score,strand)

# Make regions 500 bp around tss in bed-like format
tss.p.500 <- data.frame(chr=tss.p$chr,
                        start=tss.p$txn.start - 500,
                        end=tss.p$txn.start - 500 + 1000,
                        name=tss.p$id,
                        score=tss.p$score,
                        strand=tss.p$strand)
tss.m.500 <- data.frame(chr=tss.m$chr,
                        start=tss.m$txn.end - 500,
                        end=tss.m$txn.end -500 + 1000,
                        name=tss.m$id,
                        score=tss.m$score,
                        strand=tss.m$strand)

# Convert to  GRanges format for comparison to G4 data
tss.p.gr <- GRanges(seqnames=tss.p.500$chr,
                  ranges=IRanges(start=tss.p.500$start,end=tss.p.500$end,),
                  strand=tss.p.500$strand,
                  mcols=data.frame(name=tss.p.500$name))
tss.m.gr <- GRanges(seqnames=tss.m.500$chr,
                    ranges=IRanges(start=tss.m.500$start,end=tss.m.500$end,),
                    strand=tss.m.500$strand,
                    mcols=data.frame(name=tss.m.500$name))

# Read in the G4 data and overlap with tss regions to get nt and ts g4s
g4.tss.p.nt <- import(con="raw_data/g4-12_plus.bw",format="bigWig",which=tss.p.gr)
g4.tss.p.ts <- import(con="raw_data/g4-12_minus.bw",format="bigWig",which=tss.p.gr)
g4.tss.m.nt <- import(con="raw_data/g4-12_minus.bw",format="bigWig",which=tss.m.gr)
g4.tss.m.ts <- import(con="raw_data/g4-12_plus.bw",format="bigWig",which=tss.m.gr)

# Build overlap matrices
ol.tss.p.nt <- as.data.frame(findOverlaps(tss.p.gr,g4.tss.p.nt))
ol.tss.p.ts <- as.data.frame(findOverlaps(tss.p.gr,g4.tss.p.ts))
ol.tss.m.nt <- as.data.frame(findOverlaps(tss.m.gr,g4.tss.m.nt))
ol.tss.m.ts <- as.data.frame(findOverlaps(tss.m.gr,g4.tss.m.ts))

