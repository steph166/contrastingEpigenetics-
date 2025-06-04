# Contrasting Epigenetics of Ixodes scapularis Populations
#### Stephanie Guzman-Valencia, Jacob Cassens, Perot Saelao, Abagail Leal, Elizabeth Lohstroh, Cristina Harvey, Brenda Galvan-Leal, Cross Chambers, Crys Wright, Sydney Orsborn, Tietjen Mackenzie, Tammi Johnson, Nicole Mehta, Michael Golding, Raul Medina, Danielle M. Tufts, Jianmin Zhong, Christopher Faulk, Jonathan Oliver, and Adela Oliva Chavez
### Suplememntary material
## Table of contents
Characterization of 5-mC patterns in Texas and Minnesota I. scapularis populations using WGBS
1. Methylation call for TX and MN reads using Bismark
2. Hyper/hypormethylated genes identification with methylKit
3. Variable sites (SNPs) that interfere with DNA methylation using WGS with FreeBayes
4. Filtering variants with VCFtools and AWK

Characterization of 5-mC patterns in Texas and Minnesota I. scapularis using Oxford nanopore technology
1. Data processing and mapping using Minimap2
2. Differentially methylated regions identification using methylKit
3. Visualization of DMR with karyotypeR

### Characterization of 5-mC patterns in Texas and Minnesota I. scapularis populations using WGBS

#### Methylation call in TX and MN reads

#### Read assessment and quality control

Before trimming, raw sequences were asssessed with FASTQC
```shell
trim_galore --fastqc "${FILE_NAME_1}.fastq.gz" "${FILE_NAME_2}.fastq.gz" -o fastqcRaw
```

Raw sequences were trimmed and filtered out using Trim_galore
```shell
trim_galore --paired -q 30 --illumina --gzip --fastqc "${FILE_NAME_1}.fastq.gz" "${FILE_NAME_2}.fastq.gz" -o trimmed
```

#### Mapping reads to reference genome (GCF_016920785.2) using Bismark. Deduplication and methylation calls were performed with Bismark using reads previously sorted with SAMtools.

```shell
bismark_genome_preparation --verbose preparation  #Before the alignment, reference genome was indexed. Preparation is a folder that contains the reference genome

bismark --multicore 8 --genome alignment/ -1 "alignment/${FILE_NAME_1}_val_1.fastq.gz" -2 "alignment/${FILE_NAME_2}_val_2.fastq.gz" --output_dir alignmentOutput
```

Sequences were sorted with SAMtools
```shell
samtools sort -n "alignmentOutput/${FILE_NAME_1}_val_1_bismark_bt2_pe.bam" -o sorted_output.bam
```

Deduplication and methylation call were perfomed with Bismark
```shell
deduplicate_bismark --output_dir deduplication_output sorted_output.bam

bismark_methylation_extractor --gzip --bedGraph --buffer_size 24G --cytosine_report --genome_folder "$REAL_PATH" deduplication_output/sorted_output.deduplicated.bam -o MethylationCall
```

#### Hyper/hypormethylated genes identification with Methylkit
```R
#Package installation in Rstudio
BiocManager::install("genomation")
BiocManager::install("GRanges")
install.packages("ggplot2")
install.packages("reshape")
install.packages("dplyr")


#Loading libraries
library("genomation")
library("GenomicRanges")
library("methylKit")
library(ggplot2)
library(reshape)
library(dplyr)
library(vegan)

#Some objects were defined previously start the analysis    
gffge <- gffToGRanges("/Users/stephanie/Desktop/CoverageFile/annotationFeaturesGFF/ncbi_dataset/data/GCF_016920785.2/genomic.gff", filter = "gene")

gfftype=as(split(gff, gff$type), "GRangesList")

#Coverage files from Bismark were loaded in Rstudio. TX was the control group and MN was considered a treatment  
myobj <- methRead(list("/Users/stephanie/Desktop/CoverageFile/MN12.cov.gz", "/Users/stephanie/Desktop/CoverageFile/MN23.cov.gz", "/Users/stephanie/Desktop/CoverageFile/MN31.cov.gz","/Users/stephanie/Desktop/CoverageFile/MN46.cov.gz", "/Users/stephanie/Desktop/CoverageFile/MN47.cov.gz","/Users/stephanie/Desktop/CoverageFile/TX1.cov.gz","/Users/stephanie/Desktop/CoverageFile/TX2.cov.gz","/Users/stephanie/Desktop/CoverageFile/TX3.cov.gz","/Users/stephanie/Desktop/CoverageFile/TX5.cov.gz","/Users/stephanie/Desktop/CoverageFile/TX7.cov.gz", "/Users/stephanie/Desktop/CoverageFile/TX8.cov.gz"
),
sample.id = list("MN12","MN23","MN31","MN46","MN47","TX1","TX2","TX3","TX5","TX7","TX8"),
assembly = "tick",
pipeline = "bismarkCoverage",
treatment = c(1,1,1,1,1,0,0,0,0,0,0),
mincov = 10)

##Methylation statistics were obtained from objects defined above##

#Percentage methylation distribution and CpG coverage histogram were obtained from all samples
getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj[[1]],plot=TRUE,both.strands=FALSE)

###Comparative analysis###

#High and low read coverages were discarded to increase the power of statistical tests
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)

#Samples were merged in one object
meth=unite(myobj, destrand=FALSE)
head(meth)

#Samples were hierarchical clustering based on methylation data
#Hierarchical clustering of methylation data was assessed using correlation and ward method
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)

#customize colors in dendrogram
hc <- clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
plot(hc, main="Cluster Dendrogram", xlab="", sub="", cex=2)
dend <- as.dendrogram(hc) 
labels_colors(dend) <- c("blue","blue","blue","blue","blue","blue","red", "red","red","red","red")
plot(dend, main="CpG methylation clustering ", xlab="", sub="", cex=2)

#Methylation data was tested using Principal Component Analysis
PCASamples(meth, screeplot=TRUE)
PCASamples(meth)

#Default colors of the labels were customized using Vegan

dataPCA <- PCASamples(meth, obj.return = TRUE)
abcd1234 <- ordiplot(dataPCA, type = "none") 
points(abcd1234, "sites", col= "red", select=c(1,2,3,4,5))
points(abcd1234, "sites", col= "blue", select=c(6,7,8,9,10,11))
text(abcd1234, "sites", pos = 4, col="red", cex=0.9, select=c(1,2,3,4,5))
text(abcd1234, "sites", pos = 4, col="blue", cex=0.9, select=c(6,7,8,11))
text(abcd1234, "sites", pos = 3, col="blue", cex=0.9, select=c(9,10))


##Finding hyper/hypomethylated base among TX and MN ticks##

#Differential methylation were calculated from TX and MN samples
myDiff=calculateDiffMeth(meth)

#HYPERMETHYLATION. Significant hypermethylated bases are retained with at least methylation differences of 25% with q-value of equal or major of 0.01
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")

#Hypermethylated data was writen into a table
mydiffhypergr<- as(myDiff25p.hyper, "GRanges")  #New object was converted to GRange object
tt <- as.data.frame(mydiffhypergr)              #GRange object was convert to data frame
write.table(tt, "/Users/stephanie/Desktop/CoverageFile/MNhyper.txt", quote = FALSE, row.names = FALSE, sep="\t", col.names = TRUE)

#HYPOMETHYLATION. Significant hypomethylated bases are retained with at least methylation differences of 25% with q-value of equal or major of 0.01
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")

#Hypomethylated data was writen into a table   
mydiffhypogr<- as(myDiff25p.hypo, "GRanges")  #New object was converted to GRange object
tf <- as.data.frame(mydiffhypogr)  #GRange objet was reshaped to data frame
write.table(tf, "/Users/stephanie/Desktop/CoverageFile/MNhypo.txt", quote = FALSE, row.names = FALSE, sep="\t", col.names = TRUE)

#Differential methylated sites were calculated retaining at least methylation differences of 25% with q-value of equal or major of 0.01
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

#Hyper/hypomethylated sites were annotated and plotted by features

#Differentially hypermethylated sites were intersected with annotated features in GRanges object
diffannhyper <- annotateWithFeatures(as(myDiff25p.hyper, "GRanges"), gfftype)

#Percentage of target elements that overlappping with features and percentage of fectures overlapping with target were calculated in hyper object
diffannhyper

#Differentially hypomethylated sites were intersected with annotated features in GRanges object
diffannhypo <- annotateWithFeatures(as(myDiff25p.hypo, "GRanges"), gfftype)

#Percentage of target elements that overlappping with features and percentage of fectures overlapping with target were calculated in hypo object
diffannhypo

#Plotting the differentially hyper and hypho methylated according with features only target elements overlapping with features
allbar <- data.frame(feature = c("gene" ,"mRNA"  , "exon"  , "CDS"  ,  "pseudogene"   ,  "tRNA" ,"lnc_RNA"  ,  "snRNA", "transcript" ,"rRNA",     "cDNA_match"  ,     "snoRNA"), hyper = c(12.99 ,11.60 ,9.74 ,7.19 ,5.10 ,0.00 ,0.93 ,0.00 ,0.00 ,0.46 ,0.00 ,0.00), hypo = c(48.19   ,    44.04    ,  7.52    ,  4.88     ,  0.93    ,   0.00     ,  4.46   ,  0.00   , 0.76 , 0.00  ,  0.00  , 0.00))

tt <- melt(allbar, id.vars = "feature")  # Melt the data frame to long format

#Create a bar plot using ggplot2
ggplot(data=tt, aes(x = feature, y = value, fill = variable)) + 
  geom_bar(stat ="identity", position = position_dodge()) + 
  geom_text(aes(label = value), vjust = -0.2, color = "black", 
            position = position_dodge(0.9), size = 3.5) +
  theme_minimal() + 
  labs(x = "Feature Type", y = "Value", fill = "Variable") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#Intersection between hypermethylated sites and genes. Hypomethylated overlapping genes are not showing in this script
head(mydiffhypergr)
head(gffge)
print(mydiffhypergr)
print(gffge)

#Overlapping were found by the following GRanges objects
overlaps <- findOverlaps(mydiffhypergr,gffge)

#Overlapping sites from hypermethylated and genomic features were extracted
overlapping_ranges_Hyper <- mydiffhypergr[queryHits(overlaps)]
overlapping_ranges_gffge <- gffge[subjectHits(overlaps)]

df_Hyper <- as.data.frame(overlapping_ranges_Hyper)  #Extraction files were coverted to data frame
df_gffge <- as.data.frame(overlapping_ranges_gffge)

head(df_Hyper)
head(df_gffge)

#Some columns were extracted in each data frame: df_gffge and df_Hyper
cols_to_extract <- c(1,10,11,17, 24,25)
df_gffge <- df_gffge[, cols_to_extract]

cols_to_extract <- c(1,2,3,7,8)
df_Hyper <- df_Hyper[, cols_to_extract]

head(df_gffge)
head(df_Hyper)

#Data frames were combined in joined_ranges
joined_ranges <- cbind(df_Hyper, df_gffge)

head(joined_ranges)

#Before writing the data to csv file, all columns were tested if columns are converted to character type  
joined_ranges <- apply(joined_ranges,2,as.character)

#Hypermethylated genes were writen a table   
write.csv(joined_ranges, file = "/Users/stephanie/Desktop/CoverageFile/MNhyper_geneInfo.csv", row.names = FALSE, quote = FALSE)
```
#### Hyper-and Hypomethylated density displayed in three different scaffolds of I. scapularis
```R
#Package installation in Rstudio
BiocManager::install("karyoploteR")
BiocManager::install("readr")
BiocManager::install("seqTools")
BiocManager::install ("GenomicRanges")
BiocManager::install("genomation")


#Loading libraries
library(readr)
library(GenomicRanges)
library(readxl)
library(rtracklayer)
library(regioneR)
library(karyoploteR)

#Hypermethylated bases detected using WGBS were loaded in Rstudio 
dmr_sigHY <- read_excel("/Users/stephanie/Desktop/PlottingKaryoploteR/MNHyperExcel.xlsx", sheet=1)
#Convert significant DMRs to GRanges
dmr_sigHY_granges <- makeGRangesFromDataFrame(dmr_sigHY, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = TRUE)

#Hypomethylated bases detected using WGBS were loaded in Rstudio
dmr_sigHYPo <- read_excel("/Users/stephanie/Desktop/PlottingKaryoploteR/MNhypoExcel.xlsx", sheet=1)
#Convert significant DMRs to GRanges
dmr_sigHYPo_granges <- makeGRangesFromDataFrame(dmr_sigHYPo, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = TRUE)

#I. scapularis genome were loaded and interested columns were extracted 
genome_sizes <- read.table("/Users/stephanie/Desktop/PlottingKaryoploteR/GCF_016920785.2_ASM1692078v2_genomic.fa.fai", header = FALSE, stringsAsFactors = FALSE)
head(genome_sizes)
genome_sizes$start <- 1  
colnames(genome_sizes) <- c("chr", "end", "junk1", "junk2", "junk3","start")
genome_sizes <- genome_sizes[, c("chr","start", "end")]

#Three scaffolds were selected for plotting and genome in GRanges format is stored in custome_genome
list_of_chromosoma= c("NW_024609835.1", "NW_024609846.1", "NW_024609857.1")
custom_genome <- toGRanges(genome_sizes)

#Settings and drawing the scaffolds
kp <- plotKaryotype( genome = genome_sizes,
                     chromosomes = list_of_chromosoma, 
                     plot.type = 2,
                     ideogram.plotter = NULL,
                     cex=1.5)

kpAddCytobandsAsLine(kp)

#formatting the labes in the plot
kpText(kp, chr=genesGRange@seqnames,
       x=genesGRange$x,
       y=0.3,
       labels=genesGRange$gene,
       data.panel = 2,
       cex=1.2,
       srt=-45)

#Gene annotation from I. scapularis were loaded and gene column were used for the label 
gff_file <- "/Users/stephanie/Desktop/PlottingKaryoploteR/genomic.gff"
gffRangedData<-import.gff(gff_file)

gffRangedData <- gffRangedData[
  gffRangedData$type == "gene", c("source","type", "gene")
] 

# Some hyper- and hypomethylated genes were filtered together with their gene annotation for being plotted
gffRangedData <- gffRangedData[
  
  ##1st block of hypermethylated genes
  
    gffRangedData$gene == "LOC8023308" |
    gffRangedData$gene == "LOC115311576" |
    gffRangedData$gene == "LOC115317095" |
    
    
  ###2nd block hypomethylated genes
    
    gffRangedData$gene == "LOC8027741" | 
    gffRangedData$gene == "LOC8024530" |
    gffRangedData$gene == "LOC8030101" |
    gffRangedData$gene == "LOC8025626" |
    gffRangedData$gene == "LOC8037629" |
    gffRangedData$gene == "LOC8033621" |
    gffRangedData$gene == "LOC8035386" |
    gffRangedData$gene == "LOC8025554" |
    gffRangedData$gene == "LOC80224204" |
    gffRangedData$gene == "LOC8041333" |
    gffRangedData$gene == "LOC121833955" 
]

genesGRange<-as(gffRangedData, "GRanges")
genesGRange$x <- genesGRange@ranges@start

#Hyper- and hypomethylated density were plot in the panel 1 and 2 with the following settings 
kpAddBaseNumbers(kp)
kpPlotDensity(kp, data=dmr_sigHY_granges, data.panel = 2, col="#E41A1C" , r0=-0.15, r1=1.5)
kpPlotDensity(kp, data=dmr_sigHYPo_granges, data.panel = 1, col="#27AEF9", r0=-0.15, r1=1.5)
```

#### Variable sites (SNPs) that interfere with DNA methylation using WGS with FreeBayes

Raw sequences were trimmed using bbduk
```shell
./bbmap/bbduk.sh in1=Minnesota_Male_4_S4_R1_001.fastq.gz in2=Minnesota_Male_4_S4_R2_001.fastq.gz out1=MNmaleclean1.fq.gz out2=MNmaleclean2.fq.gz -Xmx38g ktrim=rl k=23 mink=11 hdist=1 tpe tbo minlength=25 ftm=5 qtrim=rl trimq=10 ziplevel=9
```

Reads were mapped to I. scapularis reference genome (GCF_016920785.2) with BWA

```shell
./bwa/bwa index -p myIndexGeno -a bwtsw GCF_016920785.2_ASM1692078v2_genomic.fa  #Index the genome of I. scapularis

./bwa/bwa mem -t 8 -v 2 myIndexGeno MNFemaleclean1.fq.gz MNFemaleclean2.fq.gz | gzip -3 > out.sam.gz
```

Aligments were processed using SAMtools and Picard
```shell
samtools view -b /staging/guzmanvalenc/out.sam.gz > /staging/guzmanvalenc/out.bam  #SAM files conversion and sorting in SAMtools
samtools sort -l 9 -o out.srt.bam out.bam
picard MarkDuplicates --INPUT sortedMNfem.bam --OUTPUT marked_duplMNfem.bam --METRICS_FILE marked_dup_metrics.txt --COMPRESSION_LEVEL 9  #Mark duplicates and index in BAM files using Picard
picard BuildBamIndex -I sortedMNfem.bam

```

Variable sites (SNPs) were identied with FreeBayes

```shell
freebayes-parallel <(fasta_generate_regions.py GCF_016920785.2_ASM1692078v2_genomic.fa.fai 100000) 36 -f GCF_016920785.2_ASM1692078v2_genomic.fa marked_duplMNfem.bam marked_duplMNmale.bam marked_duplTXFem.bam marked_duplTXMale.bam >var.vcf  #Calling the variants using Freebayes for one file

```

#### Filtering variants with VCFtools and AWK

```shell
vcftools --vcf var2.vcf --remove-indels --max-alleles 2 --min-alleles 2 --recode --recode-INFO-all --out filtered_biallelic_snps  #Removes indel sites and keep only biallelic SNPs
vcftools --vcf filtered_biallelic_snps.recode.vcf --minDP 4 --recode --recode-INFO-all --out filtered_depth #Filter genotypes with read depth less than 4 
vcftools --vcf filtered_depth.recode.vcf --max-missing 0.9 --recode --recode-INFO-all --out filtered_missing #Retaining site where at least 90% of sample 
vcftools --vcf filtered_missing.recode.vcf --minQ 20 --recode --recode-INFO-all --out filtered_quality  #Retain variants with a quality score of at least 20
vcftools --vcf filtered_quality.recode.vcf --mac 2 --recode --recode-INFO-all --out filtered_mac #Remove variants with a minor allele count less than 2
#Remove SNPs with Excessive coverage
vcftools --vcf filtered_mac.recode.vcf --site-mean-depth --out depth_summary  #fist calculate the mean coverage per sit
awk '{sum += $3} END {print "Mean Depth:", sum / NR}' depth_summary.ldepth.mean  #Calculate the mean of depth using awk
awk '{sum += $3} END {mean = sum / NR; print "Low Threshold:", mean * 0.5, "\nHigh Threshold:", mean * 2}' depth_summary.ldepth.mean  #Calculate the low and high threshold in awk
vcftools --vcf filtered_mac.recode.vcf --max-meanDP 43.3874 --recode --recode-INFO-all --out filtered_coverage ##Remove SNPs with coverage exceeding the high threshold in VCFtools
```


Retrieve variables sites that are present in genes of I. scapularis

```R
#Loading libraries
library("GenomicRanges")
library(ggplot2)
library(reshape)
library(rtracklayer)   
library(VariantAnnotation)
library(tidyr)

gffge <- gffToGRanges("genomic.gff", filter = "gene")   #Define the gff file that contains only genes
cols_to_extract <- c(1,5,6,12,19,20)  #Remove some columns that are no needed from the analysis
gffge <- gffge[, cols_to_extract]
gffge <- expand(gffge, colnames="Dbxref")  #Unlist the list-like  "Dbxref"
vcf_file <- "filtered_coverage.recode.vcf"  #Using VariantAnnotation::readVcf() to read the VCF file
vcf <- readVcf(vcf_file, genome = "IxodesSca")
vcf_gr <- rowRanges(vcf)  #Convert VCF data to a GRanges object
vcf_gr <- expand(vcf_gr, colnames="ALT")  #Unlist the list-like  "ALT"
#Use the findOverlaps() function to identify overlaps between the VCF and GFF regions
overlaps <- findOverlaps(vcf_gr, gffge)
#Retrieve the overlapping entries
vcf_overlapping <- vcf_gr[queryHits(overlaps)]
gff_overlapping <- gffge[subjectHits(overlaps)]
#Creating data frame
df_vcf <- as.data.frame(vcf_overlapping, row.names = NULL, optional = TRUE)
df_gffge <- as.data.frame(gff_overlapping, row.names = NULL, optional = TRUE)
#joining the dataframes and writen a table
df_joined_ranges <- cbind(df_vcf, df_gffge, stringsAsFactors = FALSE)
write.csv(df_joined_ranges, file = "joined_rangesSNPS.csv")
```
### Characterization of 5-mC patterns in Texas and Minnesota I. scapularis using Oxford nanopore technology

#### Data processing and mapping using Minimap2

#### Methylation 

Concatenate fastq files
As a pre-processing step, all fastq files generated from the sequencing experiment were concatenated into a single fastq file for each individual. Downstream analyses use the concatenated fastq file
```shell
cat *.fastq.gz > iscap_MN_meth.fastq.gz
```

Mapping reads to reference genome
Raw reads were mapped to the PalLabHifi reference genome using Minimap2. Mapped sam files were converted to bam files, sorted, and indexed using samtools
```shell
# Align to reference
minimap2 -ax map-ont pal_reference.fasta iscap_MN_meth.fastq.gz -t 32 > iscap_MN_align_meth_pal.sam
# Convert to sam to bam file 
samtools view -S -b iscap_MN_align_meth_pal.sam > iscap_MN_meth_mapped_pal.bam
# Sort mapped bam file to removed unmapped reads
samtools view -b -F 4 iscap_MN_meth_mapped_pal.bam > iscap_MN_meth_aligned_pal.bam
```

Converting to bedMethyl format
The mapped BAM file containing modified base information was converted to bedMethyl format using Modkit to extract methylation counts
```shell
modkit pileup --log-filepath modkit.log -t 28 -r pal_reference.fasta --cpg --combine-strands --only-tabs --prefix iscap_MN-HP --partition-tag HP iscap_MN_meth_aligned_pal.bam modkit_longphase
```
modkit output contains both 5mC and 5hmC. We extracted 5mC counts with awk -v OFS=“ ‘$4 ==“m”’ HP.bed

#### Differentially methylated regions identification using methylKit

#### Methylkit
Filtering and normalization
The following steps provide an example for one comparison group (MN and PA versus TX), although this analysis was repeated for three other comparison groups

```shell
# Load packages
library(tidyverse)
library(methylKit)
## 5mC analysis comparing Minnesota and Pennsylvania individuals to the Texas individual
# List the modkit files
file.list.IS <- list("/Volumes/Xtreme/aim2/methylation/input_files/IS_MN_modkit_5mc.bed.gz",
                      "/Volumes/Xtreme/aim2/methylation/input_files/IS_PA_modkit_5mc.bed.gz",
                      "/Volumes/Xtreme/aim2/methylation/input_files/IS_TX_modkit_5mc.bed.gz")
sample.id.IS <- list("IS_MN", "IS_PA", "IS_TX")
# Specify the modkit column IDs
modkit_cols <- list(
  fraction = FALSE,
  chr.col = 1,
  start.col = 2,
  end.col = 3,
  coverage.col = 5,
  strand.col = 6,
  freqC.col = 11
)
# Read in methylation data
methylationData.IS <- methRead(file.list.IS,
                                sample.id = sample.id.IS,
                                assembly = "ASM1692078v2",  # Specify your reference genome assembly
                                treatment = c(1, 1, 0),
                                pipeline = modkit_cols,
                                header = F,
                                dbtype = "tabix",
                                dbdir = "/Volumes/Xtreme/aim2/methylation/methylkit_out_IS/tabix",
                                context = "CpG",  # Context of methylation
                                mincov = 10)  # Minimum coverage to filter sites

# Merging and filtering
mk_meth.IS <- methylationData.IS |>
  filterByCoverage(lo.count = 10,
                   hi.perc = 99.9,
                   save.db = F) |>
  normalizeCoverage(save.db = F) |>
  unite(save.db = F)

```

Call differentially methylated regions (DMRs)
Count filtering, tiling into 100 bp windows, and differential methylation analysis were performed with MethylKit according to Flack et al., 2024.

```R
# DMRs
mk_meth.all <- methylationData.all |>
  unite(min.per.group = 1L, save.db = FALSE)
tile_meth.all <- tileMethylCounts(mk_meth.all,
                              win.size = 100,
                              step.size = 100,
                              mc.cores = 8,
                              cov.bases = 10,
                              save.db = F)
tile_diff.all <- calculateDiffMeth(tile_meth.all,
                               slim = F,
                               save.db = F)
tile_delta.all <- getMethylDiff(tile_diff.all,
                            difference = 50,
                            qvalue = 0.05)
# Summary statistics
dmr_summary.all <- getData(tile_delta.all) |>
  summarize(n_DMR = n(),
            mean_abs_delta = mean(abs(meth.diff)),
            min_delta = min(meth.diff),
            max_delta = max(meth.diff)) |>
  mutate(across(everything(), ~ round(.x, digits = 1))) |>
  bind_cols(data.frame(windows_tested = nrow(getData(tile_diff.all)),
                       window_size = 100)) |>
  t() |>
  data.frame() %>%
  rename_with(function(x){x = "value"}) |>
  rownames_to_column(var = "stat")
# Write table with DMR information
write.table(tile_delta.IS_all[, c("chr", "start", "end")],
            file = "dmr_data_IS.bed",
            sep = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

```
#### Visualization of DMR with karyotypeR

#### Visualization

Differential methylation was visualized by the fourteen longest scaffolds in the PalLabHifi Ixodes scapularis genome using karyotypeR. Once again, this visualization is for one of the comparison groups, and the analysis was repeated for three other comparison groups

```R
# Load packages
library(karyoploteR)
# Read in raw dmr data
dmr_all <- read_excel("/Volumes/Xtreme/aim2/nuc_genome/excel/5mc_dmr_analysis.xlsx", sheet="dmr_data_IS_all")
# Convert raw dmr data to GRanges
dmr_granges <- makeGRangesFromDataFrame(
  dmr_all,
  seqnames.field = "chr",  # Column with chromosome/scaffold names
  start.field = "start",   # Column with start positions
  end.field = "end",       # Column with end positions
  keep.extra.columns = TRUE  # Keep other columns if needed
)
# Read in significant DMRs
dmr_sig <- read_excel("/Volumes/Xtreme/aim2/nuc_genome/excel/5mc_dmr_analysis.xlsx", sheet="dmr_data_IS_sig")
# Convert significant DMRs to GRanges
dmr_sig_granges <- makeGRangesFromDataFrame(
  dmr_sig,
  seqnames.field = "chr",  # Column with chromosome/scaffold names
  start.field = "start",   # Column with start positions
  end.field = "end",       # Column with end positions
  keep.extra.columns = TRUE  # Keep other columns if needed
)
# Read in gene information
iscap_genes <- read_excel("/Volumes/Xtreme/aim2/nuc_genome/excel/high_impact_variants.xlsx", sheet="PalLabHifi_gene_annotation_IDs")
# Convert gene information to GRanges
iscap_genes_granges <- makeGRangesFromDataFrame(
  iscap_genes,
  seqnames.field = "CHROM",  # Column with chromosome/scaffold names
  start.field = "BIN_START",   # Column with start positions
  end.field = "BIN_END",       # Column with end positions
  keep.extra.columns = TRUE  # Keep other columns if needed
)
# Create the plot
kp <- plotKaryotype(genome = iscap.genome, chromosomes = c("NW_024609835", "NW_024609846", "NW_024609857", "NW_024609868", "NW_024609879",
                                                           "NW_024609880", "NW_024609881", "NW_024609883", "NW_024609882", "NW_024609884",
                                                           "NW_024609836", "NW_024609839", "NW_024609837", "NW_024609838"), plot.type = 2)
kpAddBaseNumbers(kp,
                 tick.dist = 50000000,  # Distance between ticks (adjust based on your genome size)
                 cex = 0.6,             # Text size for labels
                 add.units = TRUE,      # Adds 'Mb' units to the labels
                 digits = 2,            # Number of digits after the decimal in labels
                 minor.ticks = 0,
                 r0=0, r0=1)  
kpPlotHorizon(kp, data = dmr_granges, 
              y = dmr_granges$`meth.diff.MN-PAvTX`,
              y0 = 0,
              col = c("#27AEF9", "#E41A1C"), # You can adjust the color scheme
              border = NA,
              data.panel = 1)
kpPlotDensity(kp, data=dmr_sig_granges, data.panel = 2, col="#FFB000")
kpPlotDensity(kp, data=iscap_genes_granges, data.panel = 2, col="#AA88FF")

```

