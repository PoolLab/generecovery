#############################################
#### 1_Genome annotation pre-processing #####
#############################################


#############################################
#### Input data and statement of purpose ####
#############################################

# Background: Generating a scRNA-seq optimized transcriptomic reference requires optimizing
# the genome annotation ("xxx.gtf") file that transcriptomic references are based on.  The 
# following three aspects of genome annotations need to be optimized:
# A) Resolving gene overlap derived read loss;
# B) Recovering intergenic reads from 3’ un-annotated exons; and
# C) Recovering intronic reads.

# After optimizing and assembling the genome annotation, you can use "cellranger mkref"
# pipeline to assemble the optimized transcriptomic reference for mapping sequencing read
# data and compiling gene-cell matrices with the "cellranger count" (or other) pipeline.

# Purpose: this script generates the input data for automated genome annotation assembly
# ("2_Optimized_annotation_assembler.R").

# Required input data:
# 1. Genome annotation file ("xxx.gtf" file, from 10x Genomics provided reference transcriptome "gene" folder: can be downloaded at "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest" or Ensembl.org if wish to customize more)
# 2. Cell Ranger aligned scRNA-seq sequencing data (.bam file)
# 3. Optional: Refseq genome annotation for rapid extension of 3' gene ends in the 10x Genomics/Ensembl genome annotation. Can be accessed from UCSC genome browser: http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/ for mice and https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/ where I downloaded hg38.ncbiRefSeq.gtf.gz for humans.


###################################################
#### A. Resolve gene overlap derived read loss ####
###################################################

#### Statement of purpose: Identify all overlapping genes based on the ENSEMBL/10x Genomics
# default genome annotation file (GTF), rank-order them according to the # of gene overlaps.
# Prioritize this gene list for manual curation focusing on exonically overlapping genes.
# Resolve gene overlaps by manual inspection of transcript and gene structures and eliminating
# rare gene-expression obscuring transcripts and pseudogenes/gene models.



#### A1: Generate a rank-ordered gene list of same-strand overlapping genes ####
################################################################################

library("rtracklayer")
library("GenomicRanges")

exonic_df<- import(con = '<location_of_genome_annotation>/name.gtf', format = "gtf") # Import the original exonic genome annotation file
exonic_df = as.data.frame(exonic_df)
genes_df = exonic_df[exonic_df$type == "gene",1:13] # Extract all "gene" entries in the genome annotation to a new variable
row.names(genes_df) = 1:nrow(genes_df)
gene_names = genes_df$gene_name
genes_df = makeGRangesFromDataFrame(genes_df, keep.extra.columns=T) # convert into granges object

overlapper = rep(FALSE, length(gene_names))
number_of_overlaps = rep(0, length(gene_names))
overlapping_genes = rep("", length(gene_names))

for (i in 1:length(gene_names)){
  a = sum(countOverlaps(genes_df, genes_df[i]))
  if (a>1){
    overlapper[i] = TRUE
    number_of_overlaps[i] = a-1
    conflict_genes = gene_names[as.logical(countOverlaps(genes_df, genes_df[i]))]
    conflict_genes = setdiff(conflict_genes, gene_names[i])
    overlapping_genes[i] = paste(conflict_genes, collapse = ', ')
  }
}

overlapping_gene_list = as.data.frame(cbind(gene_names, number_of_overlaps, overlapping_genes))[overlapper,]
colnames(overlapping_gene_list) = c("Gene", "Number_of_gene_overlaps", "Overlapping_genes")
overlapping_gene_list$Number_of_gene_overlaps = as.integer(overlapping_gene_list$Number_of_gene_overlaps)

o = order(overlapping_gene_list$Number_of_gene_overlaps, decreasing = TRUE) # Rank order genes by the number of gene overlaps
overlapping_gene_list = overlapping_gene_list[o,]
row.names(overlapping_gene_list) = 1:nrow(overlapping_gene_list)

write.csv(overlapping_gene_list, "overlapping_gene_list.csv")

#### A2: Prioritize gene list by classifying overlapping into recommended action categories ####
################################################################################################

# Optional: R-implementation coming shortly (current version in Python 3)

#### A3: Manually curate gene overlap list ####
###############################################

# Inspect gene overlaps and transcript structure in Ensembl
#(https://nov2020.archive.ensembl.org/Mus_musculus/Info/Index for mice and 
# https://useast.ensembl.org/Homo_sapiens/Info/Index for humans) to determine which
# readthrough and premature start transcripts to eliminate (mark under "transcripts_for_deletion"
# column in "overlapping_gene_list.csv"). Also, for pseudogenes/poorly supported gene models
# that obscure protein coding genes, examine evidence for its existence (are they reported
# by other annotation consortia - e.g. Refseq). You can also examine if any reads map to it
# with regular exonic reference and determine whether to keep or eliminate it from the annotation.
# To mark gene for deletion, mark its status as "delete" under column "final_classification"
# in the overlapping_gene_list.csv.


################################################################
#### B. Recover intergenic reads from 3’ un-annotated exons ####
################################################################


#### Instantiate libraries and define inputs

library("GenomicFeatures")
library("GenomicRanges")
library("rtracklayer") # For coercing data to bed format
library("IRanges")

## Inputs:
# 1. Genome annotation file ("xxx.gtf" file, from 10x Genomics provided reference transcriptome "gene" folder: can be downloaded at "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest" or Ensembl.org if wish to customize more)
# 2. Regular reference transcriptome aligned sequencing data (note, needs to be aligned with Cell Ranger 4.00 or later to take advantage of Cell Ranger's classification of sequencing read in the bam tag "RE"). Also, note that different scRNA-seq
# datasets may have different sets of genes with unincorporated intergenic sequencing reads 3' of known ends. 
# 3. Optional: Refseq genome annotation for rapid extension of 3' gene ends in the 10x Genomics/Ensembl genome annotation. Can be accessed from UCSC genome browser: http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/ for mice and https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/ where I downloaded hg38.ncbiRefSeq.gtf.gz for humans

## Outputs:
# A. Gene_3'_intergenic_read_summary

## Steps:
# B1: Isolate intergenic sequencing reads
# B2: Create gene ranges file for linking reads to genes
# B3: Identify candidate genes for extension with excess 3' intergenic reads
# B4: Add Ensembl/10x Genomics and Refseq 3' gene ends to gene extension decision file
# B5: Manual curation of candidate genes for 3' gene extension

#### B1: Isolate intergenic sequencing reads ####
#################################################

# Intergenic reads are extracted from Cell Ranger aligned bam file. Intergenic reads can be identified by two features: their read identity tag RE = "I" (for intergenic) OR
# their RE=E (for exonic) with AN = <some gene>. The latter reads are in fact intergenic reads since Cell Ranger wrongly classifies reads mapping antisense to an exon as exonic (i.e. RE="E").
# The false exonic reads can be recognized and captured as proper intergenic reads by extracting two kinds of reads (RE=I and RE=E & AN=<something else than NA).
# Also, removing duplicates command in GenomicAlignments package does not work for intergenic (nor for intronic) reads. Duplicate and corrupt read removal has to be done manually
# (i.e. make sure cellular and molecular barcodes have specified lengths and duplicate barcodes removed).

## Load sequencing data from transcriptome aligned .bam file (note, needs to be aligned with Cell Ranger!).

bamfile = file.path("<file_location>/possorted_genome_bam.bam")
indexfile = file.path("<file_location>/possorted_genome_bam.bam") #Note, you don't have to specify ".bai" extension here.
seq_data = readGAlignments(bamfile, index=indexfile, param = ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE), tag = c("GN", "RE", "CB", "UB", "AN"), what = "flag", tagFilter = list("RE"=c("I", "E")))) # Extract all non-duplicate intergenic and exonic sequening reads with the following bam tags: "GN" - aligned gene; "RE" - read classification into E-exonic, N-intronic, I-intergenic; "CB" - corrected cellular barcode; "UB" - corrected UMI/molecular barcode, "AN" - antisense gene.
seq_data = data.frame(seq_data)

## Keep only intergenic reads by removing all true exonic reads (i.e. remove exonic reads that lack antisense gene mapping: AN tag =NA)
intergenic_reads = seq_data$RE=="I"
false_exonic_reads = !is.na(seq_data$AN)
all_intergenic_reads = as.logical(intergenic_reads + false_exonic_reads)

seq_data = seq_data[all_intergenic_reads,] # remaining dataframe contains only intergenic reads. Note, that we are assuming that all false exonic reads are intergenic, which slightly overestimates intergenic read count. This is since some will likely also end up being intronic.

## Remove all duplicate reads and reads with corrupt barcodes (i.e. keep reads with 16 nucleotide cellular barcodes and 10 nucleotide molecular barcodes). Note that duplicate removal is required since Cell Ranger does not automatically flag duplicates for intronically and intergenically classified reads.

library(stringr)
seq_data$CB = str_sub(seq_data$CB, end=-3) # Remove last two elements of the cell barcode. This is an artifact ("-1") added by Cell Ranger software.
seq_data$barcodes = paste(seq_data$CB, seq_data$UB, sep="") # Assemble the cell barcode / molecular barcode list. Each read included in the gene_cell matrix will have a unique index comprised of the two.

a = nchar(seq_data$barcodes)==26 # logical vector for selecting reads with non-corrupt barcodes
seq_data = seq_data[a,] # exclude all reads that don't have an intact full cellular and molecular barcodes
length(unique(seq_data$barcodes)) # Determine # of unique intergenic reads
seq_data = seq_data[!duplicated(seq_data$barcodes),] # exclude all duplicated intergenic reads

## Save extracted intergenic reads as a separate file
gr_seq_data = makeGRangesFromDataFrame(seq_data) # coerce to granges object as that makes it possible to save it as a bedfile that bedtools can parse

ga_seq_data = as(gr_seq_data, "GAlignments") # Coerces GRanges object into a GAlignments object, that can be saved as a bed file. Required for bedtools to link reads to closest 3' gene end.
m = asBED(ga_seq_data)# converts GAlignments object into the bed format
export.bed(m, con = "<save_location>/intergenic_reads.bed")


#### B2: Create gene ranges file for linking reads to genes (note: partially in BASH/Terminal) ####
###################################################################################################

#### Purpose: Make a bed file with gene boundaries, which is required for assigning intergenic reads to a specific gene. This step is in bash with one step in R.

#### Input: Genome annotation file ("xxx.gtf" file, from 10x Genomics provided reference transcriptome "gene" folder: can be downloaded at "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest" or Ensembl.org if wish to customize more)

#### Notes:
# Part of this step should be run in linux Terminal in Bash.
# This step requires bedtools: https://bedtools.readthedocs.io/en/latest/content/installation.html. Put it in PATH after installing.


## In bash/linux terminal: Create bed file with gene boundaries from the gene annotation ("xxx.gtf") file (perform in linux Terminal in bash!)

# In linux terminal, navigate to folder with the genome annotation of interest. Assuming that it is named "genes.gtf" per 10x Genomics convention.

grep -P '\tgene\t' genes.gtf > gene_ranges.gtf # extracts all "gene" entries from genome annotation and saves to a new file 

## In bash/linux terminal: Add "transcript_id """ column to the gtf file to make it compatible with bedtools format
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' gene_ranges.gtf > gene_ranges1.gtf # can discard the original, can discard
rm gene_ranges.gtf # removes intermediate file

## In bash: Convert reference gtf into bed with bedtools
export PATH=<path_to_bedtools_bin_folder>:$PATH  # Include bedtools location to path.
gtf2bed < gene_ranges1.gtf > gene_ranges.bed # Creates a bed file with gene boundaries

## In R: Replace final column with gene name. Make sure you navigate to same folder in R.

library(stringr)
gene_ranges = read.table("gene_ranges.bed", sep = "\t")

for (i in 1:dim(gene_ranges)[1])
{
  a = gene_ranges[i,10]
  res <- str_match(a, "gene_name\\s*(.*?)\\s*;")
  b = res[,2]
  gene_ranges[i, 10] = b
}

## In R: save outcome
write.table(gene_ranges, "gene_ranges.bed", sep="\t",row.names=FALSE, col.names=FALSE, quote = FALSE)


#### B3: Identify candidate genes for extension with excess 3' intergenic reads ####
####################################################################################

## In bash/linux terminal: Make sure bedtools is in PATH.

export PATH=<location of bedtools bin folder>:$PATH # Places bedtools in path.

sortBed -i intergenic_reads.bed > intergenic_reads1.bed # sort intergenic reads file
sortBed -i gene_ranges.bed > gene_ranges1.bed  # sort gene ranges file

bedtools closest -a intergenic_reads1.bed -b gene_ranges1.bed -s -D a -fu > results.txt # resulting file contains sequencing reads with distance data from closest 3' gene identity and end

## In R: Save a rank ordered list of genes with highest-to-lowest number of intergenic reads within 10kb of its known gene end.

summary_data = read.table("results.txt", sep = "\t")

summary_data = summary_data[summary_data$V23>-10000,] # retain only sequencing reads within 10kb of known gene ends. Change to more or less stringent as desired.
summary_data = summary_data[summary_data$V23<0,] # retain only sequencing reads within 10kb of known gene ends. Change to more or less stringent as desired.
hist(summary_data$V23) # plot histogram of intergenic sequencing reads as a function of distance from 3' gene ends.

summary_data_genes = table(summary_data$V22) # Summarizes # of intergenic reads within 10kb of known gene ends for each gene.
o = order(summary_data_genes, decreasing = TRUE) # Rank order the gene list
length(summary_data_genes)
summary_data_genes = summary_data_genes[o]
length(summary_data_genes[summary_data_genes>10]) # Determine number of genes with more than 10 intergenic reads within 10kb of known gene end
summary_data_genes = summary_data_genes[summary_data_genes>10] # Threshold gene list based on the amount of intergenic gene loading.
summary_data_genes = data.frame(summary_data_genes)
dim(summary_data_genes)
summary_data_genes[1:40,]

write.csv(summary_data_genes, "gene_extension_candidates.csv") # Saves a rank ordered list of genes as a function of 3' intergenic read mapping within 10kb of known gene end. You can use this as a prioritized gene list for gene extension to examine in Integrated Genomics Viewer.

summary_data_genes["<gene_of_interest"] # display 3' intergenic read mapping for a gene of interest



#### B4: Add Ensembl/10x Genomics and Refseq 3' gene ends to gene extension decision file ####
##############################################################################################

#### Acquire ensembl gene boundaries from 10x Genomics/Ensembl genome annotation. The genome annotations can be acquired from 10x Genomics provided reference transcriptome "gene" folder: can be downloaded at "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest"

ensembl_db <- import(con ="<location of genome annotation>/genes.gtf", format = "gtf")

ensembl_db = as.data.frame(ensembl_db)
ensembl_db = ensembl_db[,c(13,2,3,5,7)] # Extract relevant information from genome annotation (gene, start, end, strand and feature)
colnames(ensembl_db) = c("gene", "start", "end", "strand", "feature")
ensembl_db = ensembl_db[ensembl_db$feature == "gene",] # Keep only "gene" features from the annotation 
head(ensembl_db)
dim(ensembl_db)

## Deal with gene duplications in ensembl genome annotation (weirdly the annotation contains several entries for some genes)

ensembl_duplicates=unique(ensembl_db$gene[duplicated(ensembl_db$gene)])
start = rep(0, length(ensembl_duplicates))
end = rep(0, length(ensembl_duplicates))
strand = as.factor(rep(c("+", "-"), length(ensembl_duplicates)/2))
feature = rep("gene", length(ensembl_duplicates))

ensembl_fix = data.frame(cbind(ensembl_duplicates, start, end, strand, feature))
ensembl_fix$strand = strand
colnames(ensembl_fix) = c("gene", "start", "end", "strand", "feature")
row.names(ensembl_fix) = ensembl_fix$gene

for (i in 1:length(ensembl_duplicates)){
  gene=ensembl_duplicates[i]
  a = ensembl_db[ensembl_db$gene == gene,]
  ensembl_fix[gene,2] = min(a$start)
  ensembl_fix[gene,3] = max(a$end)
  ensembl_fix[gene,4] = as.factor(a[1,4])
}

ensembl_db = ensembl_db[!(ensembl_db$gene %in% ensembl_duplicates), ] # kick out all the duplicate gene entries to append them later at the end

ensembl_db = rbind(ensembl_db, ensembl_fix)
ensembl_db = data.frame(ensembl_db[,-1], row.names = ensembl_db[,1])



#### Get Refseq gene boundaries (for many genes Refseq has annotated 3' gene boundaries to extend significantly further than in the Ensembl/Gencode references. This allows for rapid extension of hundreds of genes). Refseq genome annotations can be accessed from UCSC genome browser: http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/ for mice and https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/ where you can download hg38.ncbiRefSeq.gtf.gz for humans.

refseq_db = import(con ="C:/Users/allan/OneDrive/Desktop/Helen Allan Computational Paper/Data/Refseq_GTFs/mm10.ncbiRefSeq.gtf", format = "gtf")
refseq_db = as.data.frame(refseq_db)
refseq_db = refseq_db[,c(12,2,3,5,7)] # Note, these positions can change with successive updates of the format.
colnames(refseq_db) = c("gene", "start", "end", "strand", "feature")
refseq_db = refseq_db[refseq_db$feature == "transcript",]
dim(refseq_db)
head(refseq_db)


## Deal with gene duplications in refseq genome annotation (similarly to ensembl reference, the Refseq reference somewhat unexpectedly contains several gene entries for some genes)

refseq_duplicates=unique(refseq_db$gene[duplicated(refseq_db$gene)])
start = rep(0, length(refseq_duplicates))
end = rep(0, length(refseq_duplicates))
strand = as.factor(rep(c("+", "-"), length(refseq_duplicates)/2))
feature = rep("gene", length(refseq_duplicates))

refseq_fix = data.frame(cbind(refseq_duplicates, start, end, strand, feature))
refseq_fix$strand = strand
colnames(refseq_fix) = c("gene", "start", "end", "strand", "feature")
row.names(refseq_fix) = refseq_fix$gene

for (i in 1:length(refseq_duplicates)){
  gene=refseq_duplicates[i]
  a = refseq_db[refseq_db$gene == gene,]
  refseq_fix[gene,2] = min(a$start)
  refseq_fix[gene,3] = max(a$end)
  refseq_fix[gene,4] = as.factor(a[1,4])
}

refseq_db = refseq_db[!(refseq_db$gene %in% refseq_duplicates), ] # kick out all the duplicate gene entries to append them later at the end

refseq_db = rbind(refseq_db, refseq_fix)
refseq_db = data.frame(refseq_db[,-1], row.names = refseq_db[,1])


#### Assemble summary file for gene extension curation

strand = as.factor(rep(c("+", "-"), 10000))[1:dim(summary_data_genes)[1]]
ensembl_start = rep(0, dim(summary_data_genes)[1])
ensembl_end = rep(0, dim(summary_data_genes)[1])
refseq_term_left = rep(0, dim(summary_data_genes)[1])
refseq_term_right = rep(0, dim(summary_data_genes)[1])
ref_seq_overhang = rep(0, dim(summary_data_genes)[1])

for (i in 1:dim(summary_data_genes)[1]){
  gene = row.names(summary_data_genes)[i]
  ensembl_start[i] = as.numeric(ensembl_db[gene,]$start)
  ensembl_end[i] = as.numeric(ensembl_db[gene,]$end)
  strand[i] = ensembl_db[gene,]$strand
}

summary_data_genes = data.frame(cbind(summary_data_genes, strand, ensembl_start, ensembl_end, refseq_term_left, refseq_term_right, ref_seq_overhang))
# summary_data$refseq_term_right = as.numeric(summary_data$refseq_term_right)

ref_seq_genes = intersect(rownames(summary_data_genes), row.names(refseq_db))

for (i in ref_seq_genes){
  if (summary_data_genes[i,2] == "+"){
    summary_data_genes[i, 6] = as.numeric(refseq_db[i, 2])
    summary_data_genes[i, 7] = summary_data_genes[i, 6] - summary_data_genes[i, 4]
  } else {
    summary_data_genes[i, 5] = as.numeric(refseq_db[i, 1])
    summary_data_genes[i, 7] = summary_data_genes[i, 3] - summary_data_genes[i, 5]
  }
}

write.csv(summary_data_genes, "gene_extension_candidates.csv") # Final decision file for gene extension curation



#### B5: Manual curation of candidate genes for 3' gene extension ####
######################################################################

# Download and install https://software.broadinstitute.org/software/igv/
# Load the transcriptome aligned sequencing read data ("xxx.bam" file) to the igv and choose appropriate species/genome.
# Go through the rank ordered candidate gene list in "gene_extension_candidates.csv" and examine evidence for unannotated 3' exons and UTRs:
#  a) extensive splicing between annotated and candidate unannotated regions at the genetic locus as visualized by igv.
#  b) continuous read mapping from known gene end (delayed cutoff of sequencing reads posterior from known gene end)
#  c) examine external evidence (e.g. Allen in situ brain atlas, Human Protein Atlas etc.)

# Based on the latter, provide new gene boundaries under two new columns: "update_start" and "update_end" depending on whether gene is on "-" or "+" strands, respectively. Use excel formulas, if necessary, to scale curation, and save file as "gene_extension_candidates.xlsx".
# Resulting manually curated "gene_extension_candidates.xlsx" will be used as input in the Optimized annotation assembly stage.

###################################
#### C. Recover intronic reads ####
###################################

# Note: several methods are available for recovering intronic reads from scRNA-seq data. These
# include the traditional approach of redefining all genes/transcripts as exons, which enables
# incorporation of intronic reads to downstream analysis. This approach works well but has
# several problems including poorer performance at mapping spliced reads (annotation loses
# splicing boundaries) and eliminates hundreds of genes due to gene overlaps. There are also
# specific modes in the latest versions of sequence aligners (e.g. --include-introns parameter
# in Cell Ranger) for incorporating intronic reads that should work better. However, we found
# that the latter systematically fails to retrieve intronic reads for many genes that the
# traditional premrna approach performs well on. We therefore adopted a hybrid pre-mRNA
# reference approach where we supplement normal gene annotation entries by traditional
# pre-mRNA entries where transcripts have been redefined as exons and map in the
# --include-introns mode to retrieve most of available intronic reads.


# Input: original genome annotation file (.gtf)
# Output: pre-mRNA genome annotation (.gtf), which contains transcripts redefined as exons. This is used as input in the final optimized genome annotation assembly.

#### Generate a basic pre-mRNA reference (define transcripts as exons).

library("rtracklayer")

exonic_df<- import(con = '<location_of_genome_annotation>/name.gtf', format = "gtf") # Import the original exonic genome annotation file
exonic_df = as.data.frame(exonic_df)
premrna_df = exonic_df[exonic_df$type == "transcript",] # Extract all "transcript" entries in the genome annotation to a new variable
premrna_df = premrna_df$feature = rep("exon", nrow(premrna_df)) # Rename all "feature" 

premrna_df = makeGRangesFromDataFrame(premrna_df, keep.extra.columns=TRUE)
rtracklayer::export(premrna_df, "premrna.gtf", format = "gtf")

