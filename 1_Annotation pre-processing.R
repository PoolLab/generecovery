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
# B) Recovering intergenic reads from 3' un-annotated exons; and
# C) Optimally recovering intronic reads.

# After optimizing and assembling the genome annotation, you can use "cellranger mkref"
# pipeline to assemble the optimized transcriptomic reference. The latter can be used  for
# mapping sequencing read data and compiling gene-cell matrices with the "cellranger count"
# (or other) pipeline.

# Purpose: this script generates the input data for automated genome annotation assembly
# ("2_Optimized_annotation_assembler.R").

# Required input data:
# 1. Genome annotation file ("xxx.gtf" file, from 10x Genomics provided reference transcriptome "gene" folder: can be downloaded at "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest" or Ensembl.org if wish to customize more)
# 2. Cell Ranger aligned scRNA-seq sequencing data (.bam file). This is required to identify genes with excessive 3' intergenic reads that may indicate the presence of an unannotated 3' exon.
# 3. Optional: Refseq genome annotation for rapid extension of 3' gene ends in the 10x Genomics/Ensembl genome annotation. Can be accessed from UCSC genome browser: http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/ for mice and https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/ where I downloaded hg38.ncbiRefSeq.gtf.gz for humans.

# Required software: This code relies on two software packages that are only available for linux/mac OS paltforms. Therefore it is highly advisable you run this analysis
# on a linux or mac OS machine.
# 1. BEDOPS: available for install at https://bedops.readthedocs.io/en/latest/content/installation.html 
# 2. Bedtools: available for install at https://bedtools.readthedocs.io/en/latest/content/installation.html

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

install.packages("BiocManager")
BiocManager::install("rtracklayer")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomicAlignments")
BiocManager::install("IRanges")
install.packages("stringr")

library("rtracklayer")
library("GenomicRanges")
library("stringr")


exonic_df<- import(con = '<location of unoptimized genome annotation>/genes.gtf', format = "gtf") # Import the original exonic genome annotation file
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
colnames(overlapping_gene_list) = c("genes", "number_of_gene_overlaps", "overlapping_genes")
overlapping_gene_list$number_of_gene_overlaps = as.integer(overlapping_gene_list$number_of_gene_overlaps)

o = order(overlapping_gene_list$number_of_gene_overlaps, decreasing = TRUE) # Rank order genes by the number of gene overlaps
overlapping_gene_list = overlapping_gene_list[o,]
row.names(overlapping_gene_list) = 1:nrow(overlapping_gene_list)

final_classification = rep("", nrow(overlapping_gene_list))
transcripts_for_deletion = rep("", nrow(overlapping_gene_list))
comments = rep("", nrow(overlapping_gene_list))

overlapping_gene_list = cbind(overlapping_gene_list, final_classification, transcripts_for_deletion, comments)

head(overlapping_gene_list)

write.csv(overlapping_gene_list, "overlapping_gene_list.csv")

#### A2 (OPTIONAL): Prioritize gene list by classifying overlapping into recommended action categories ####
###########################################################################################################

# Packages to import
library("rtracklayer")
library("stringr")

# Import data
data = read.csv("./overlapping_gene_list.csv", header=T)

# Import genes.gtf file for gene overlap scenario determination
genes_gtf <- import(con = '<location of the unoptimized genome annotation>/genes.gtf', format = "gtf") # Import the original exonic genome annotation file
genes_gtf = as.data.frame(genes_gtf)

# Import overlapping genes list from A1
data$overlapping_genes <- strsplit(data$overlapping_genes, ", ") 

# FUNCTIONS
return_exons <- function(gene_name){
  exon_subset <- subset(gene_name, type == 'exon')
  return(data.frame(exon_subset['start'], exon_subset['end']))
}

both_pseudo <- function(key, overlapping){
  return((str_sub(key, 1, 2) == 'Gm' | str_sub(key, - 3, - 1) == 'Rik') & (str_sub(overlapping, 1, 2) == 'Gm' | str_sub(overlapping, - 3, - 1) == 'Rik'))
}

exon_overlap <- function(gene_A_exons, gene_B_exons){
  for(row_exonA in 1:nrow(gene_A_exons)){
    for(row_exonB in 1:nrow(gene_B_exons)){ 
      x = seq(from = gene_A_exons[row_exonA,1], to = gene_A_exons[row_exonA,2]-1, by = 1)
      y = seq(from = gene_B_exons[row_exonB,1], to = gene_B_exons[row_exonB,2]-1, by = 1)
      
      if(length(intersect(x,y))!=0){
        return (TRUE)
      }
    }
  }
  return (FALSE)
}

pseudo_overlap <- function(key, overlapping, gene_A_exons, gene_B_exons){
  # Check for exon overlap
  if(str_sub(key, 1, 2) == 'Gm' | str_sub(key, - 3, - 1) == 'Rik' | str_sub(overlapping, 1, 2) == 'Gm' | str_sub(overlapping, - 3, - 1) == 'Rik'){
    if(exon_overlap(gene_A_exons, gene_B_exons) == TRUE){
      # Check if gene_A is a pseudogene
      if(str_sub(key, 1, 2) == 'Gm' | str_sub(key, - 3, - 1) == 'Rik'){
        return(key)
      }
      else{
        return(overlapping)
      }
    }
    else{
      return('exonic')
    }
  }
  
  else{
    return('empty')
  }
}

readthrough_or_premature <- function(upstream_name, downstream_name, upstream_trx, downstream_trx){
  max_u = 0
  for(trx_u in 1:length(upstream_trx[[1]])){
    count_u = 0
    for (trx_d in 1:length(downstream_trx[[1]])){
      x = seq(downstream_trx[[1]][trx_d], downstream_trx[[2]][trx_d]-1)
      y = seq(upstream_trx[[1]][trx_u], upstream_trx[[2]][trx_u]-1)
      
      if(length(intersect(x,y)) > 0){
        count_u = count_u + 1
      }
      
      if(count_u > max_u){
        max_u = count_u
      }
    }
  }
  
  max_d = 0
  for(trx_d in 1:length(downstream_trx[[1]])){
    count_d = 0
    for(trx_u in 1:length(upstream_trx[[1]])){
      x = seq(downstream_trx[[1]][trx_d], downstream_trx[[2]][trx_d]-1)
      y = seq(upstream_trx[[1]][trx_u], upstream_trx[[2]][trx_u]-1)
      if(length(intersect(x,y)) > 0){
        count_d = count_d + 1
      }
    }
    if(count_d > max_d){
      max_d = count_d
    }
  }
  
  if(max_u > max_d){
    result <- list(upstream_name, downstream_name, "readthrough")
    return(result)
  }
  else if(max_d > max_u){
    result <- list(upstream_name, downstream_name, "premature")
    return(result)
  }
  else if(min(length(upstream_trx), length(downstream_trx)) == 1 & max(length(upstream_trx), length(downstream_trx)) != 1){
    print("MAX & MIN")
    result <- list(upstream_name, downstream_name, "manual")
    return(result)
  }
  else if(max_u == 1 & max_d == 1){
    result <- list(upstream_name, downstream_name, "manual")
    return(result)
  }
  else if(max_u == max_d){
    if(length(upstream_trx[[1]]) > length(downstream_trx[[1]])){
      result <- list(upstream_name, downstream_name, "readthrough")
      return(result)
    }
    else{
      result <- list(upstream_name, downstream_name, "premature")
      return(result)
    }
  }
  else{
    result <- list(upstream_name, downstream_name, "manual")
    return(result)
  }
}

readthrough_or_premature_plus <- function(name_A, gene_A, name_B, gene_B, gene_A_exons, gene_B_exons){
  select_gene_A <- gene_A[gene_A$type=='gene', ]
  gene_A_start <-  min(select_gene_A[,'start'])
  gene_A_end <- max(select_gene_A[,'end'])
  select_gene_B <- gene_B[gene_B$type=='gene', ]
  gene_B_start <-  min(select_gene_B[,'start'])
  gene_B_end <- max(select_gene_B[,'end'])
  
  if(gene_A_start < gene_B_start){
    upstream_name <- name_A
    upstream <- gene_A
    downstream_name <- name_B
    downstream <- gene_B
  }
  
  else if(gene_B_start < gene_A_start){
    upstream_name <- name_B
    upstream <- gene_B
    downstream_name <- name_A
    downstream <- gene_A
  }
  
  else if(gene_A_end < gene_B_end){
    upstream_name <- name_A
    upstream <- gene_A
    downstream_name <- name_B
    downstream <- gene_B
  }
  
  else{
    upstream_name <- name_B
    upstream <- gene_B
    downstream_name <- name_A
    downstream <- gene_A
  }
  
  upstream_trx <- list(upstream[upstream$type == 'transcript',][,'start'], upstream[upstream$type == 'transcript',][,'end'])
  downstream_trx = list(downstream[downstream$type == 'transcript',][,'start'], downstream[downstream$type == 'transcript',][,'end'])
  
  return(readthrough_or_premature(upstream_name, downstream_name, upstream_trx, downstream_trx))
}

readthrough_or_premature_min <- function(name_A, gene_A, name_B, gene_B, gene_A_exons, gene_B_exons){
  select_gene_A <- gene_A[gene_A$type=='gene', ]
  gene_A_start <-  min(select_gene_A[,'start'])
  gene_A_end <- max(select_gene_A[,'end'])
  select_gene_B <- gene_B[gene_B$type=='gene', ]
  gene_B_start <-  min(select_gene_B[,'start'])
  gene_B_end <- max(select_gene_B[,'end'])
  
  if(gene_B_end > gene_A_end){
    upstream_name <- name_B
    upstream <- gene_B
    downstream_name <- name_A
    downstream <- gene_A
  }
  
  else if(gene_A_end > gene_B_end){
    upstream_name <- name_A
    upstream <- gene_A
    downstream_name <- name_B
    downstream <- gene_B
  }
  
  else if(gene_B_start > gene_A_start){
    upstream_name <- name_B
    upstream <- gene_B
    downstream_name <- name_A
    downstream <- gene_A
  }
  
  else{
    upstream_name <- name_A
    upstream <- gene_A
    downstream_name <- name_B
    downstream <- gene_B
  }
  
  upstream_trx <- list(upstream[upstream$type == 'transcript',][,'start'], upstream[upstream$type == 'transcript',][,'end'])
  downstream_trx = list(downstream[downstream$type == 'transcript',][,'start'], downstream[downstream$type == 'transcript',][,'end'])
  
  return(readthrough_or_premature(upstream_name, downstream_name, upstream_trx, downstream_trx))
}

# CLASSIFICATION
for(key in (rownames(data))){
  
  # Check that the gene is not classified already
  if(is.na(data[key,'automatic_classification'])){
    gene_A <- subset(genes_gtf, gene_name == key)
    
    if(data[key,'number_of_gene_overlaps'] > 1){
      overlaps <- data[key,'overlapping_genes']
      for(item in overlaps[[1]]){
        gene_B = genes_gtf[genes_gtf['gene_name'] == item,]
        
        gene_A_exons = return_exons(gene_A)
        gene_B_exons = return_exons(gene_B)
        
        if(exon_overlap(gene_A_exons, gene_B_exons) == TRUE){
          data[item, 'automatic_classification'] = 'Manual inspection'
          
          if(is.na(data[key,'automatic_classification']) | data[key,'automatic_classification'] != 'Manual inspection'){
            data[key,'automatic_classification'] = 'Manual inspection'
          }
        }
        else{
          if(is.na(data[key,'automatic_classification'])){
            data[key,'automatic_classification'] = 'Keep as is'
            
            if(data[item,'number_of_gene_overlaps'] > 1){
              data[item,'automatic_classification'] = 'Manual inspection'
            }
            else{
              data[item,'automatic_classification'] = 'Keep as is'
            }
          }
        }
      } 
    }
    
    if(data[key,'number_of_gene_overlaps'] == 1){
      overlapping <- data[key,'overlapping_genes'][[1]]
      gene_B <- subset(genes_gtf, gene_name == overlapping)
      strand <- gene_A[1,'strand']
      
      gene_A_exons = return_exons(gene_A)
      gene_B_exons = return_exons(gene_B)
      
      # Check if both - key and overlapping gene - are pseudogenes
      if(both_pseudo(key, overlapping) == TRUE){
        data[key, 'automatic_classification'] = 'Manual inspection'
        data[overlapping[[1]], 'automatic_classification'] = 'Manual inspection' 
      }
      
      #  Check for pseudogene
      else if(pseudo_overlap(key, overlapping, gene_A_exons, gene_B_exons) == key){
        data[key, 'automatic_classification'] = 'Delete'
        data[overlapping[[1]], 'automatic_classification'] = 'Keep as is'
      }
      
      else if(pseudo_overlap(key, overlapping, gene_A_exons, gene_B_exons) == overlapping){
        data[key, 'automatic_classification'] = 'Keep as is' 
        data[overlapping[[1]], 'automatic_classification'] = 'Delete'
      }
      
      else if(pseudo_overlap(key, overlapping, gene_A_exons, gene_B_exons) == 'exonic'){
        data[key, 'automatic_classification'] = 'Keep as is'
        data[overlapping[[1]], 'automatic_classification'] = 'Keep as is' 
      }
      
      # Check for readthrough
      else if(exon_overlap(gene_A_exons, gene_B_exons) == TRUE){
        if(strand == '+'){
          name_A = key
          name_B = overlapping
          result = readthrough_or_premature_plus(name_A, gene_A, name_B, gene_B, gene_A_exons, gene_B_exons)
          
          if(result[[3]] == 'readthrough'){
            data[result[[1]],'automatic_classification'] = 'Readthrough transcript deletion'
            data[result[[2]],'automatic_classification'] = 'Keep as is' 
          }
          else if(result[[3]] == 'premature'){
            data[result[[1]],'automatic_classification'] = 'Keep as is'
            data[result[[2]],'automatic_classification'] = 'Premature transcript deletion' 
          }
          else if(result[[3]] == 'manual'){
            data[result[[1]],'automatic_classification'] = 'Manual inspection' 
            data[result[[2]],'automatic_classification'] = 'Manual inspection'
          }
        }
        
        else if(strand == '-'){
          name_A = key
          name_B = overlapping
          result = readthrough_or_premature_min(name_A, gene_A, name_B, gene_B, gene_A_exons, gene_B_exons)
          
          if(result[[3]] == 'readthrough'){
            data[result[[1]],'automatic_classification'] = 'Readthrough transcript deletion'
            data[result[[2]],'automatic_classification'] = 'Keep as is' 
          }
          else if(result[[3]] == 'premature'){
            data[result[[1]],'automatic_classification'] = 'Keep as is'
            data[result[[2]],'automatic_classification'] = 'Premature transcript deletion'
          }
          else if(result[[3]] == 'manual'){
            data[result[[1]],'automatic_classification'] = 'Manual inspection' 
            data[result[[2]],'automatic_classification'] = 'Manual inspection'
          }
        }
      }
      
      else if(exon_overlap(gene_A_exons, gene_B_exons) == FALSE){
        data[key,'automatic_classification'] = 'Keep as is' 
        data[overlapping,'automatic_classification'] = 'Keep as is'
      }
    }    
  }
}

# Convert data$overlapping_genes from list to string
string_overlapping = rep("", length(data$overlapping_genes))

for (i in 1:length(data$overlapping_genes)){
  a = unlist(data$overlapping_genes[i])
  a = paste(a, collapse = " ")
  string_overlapping[i] = a
}

data$overlapping_genes = string_overlapping


# SAVE DATA
write.csv(data,"./overlapping_gene_list.csv")


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
# To mark gene for deletion, mark its status as "Delete" under column "final_classification"
# in the overlapping_gene_list.csv.


# For genes that cannot be unentangled due to overlap of terminal exons, delete one of them and rename the other to reflect that reads
# could orignate from both. To this end, mark one of the genes for deletion (add "Delete" under "overlapping_gene_list.csv" file $final_classification)
# field. Also, create a new csv file "rename_genes.csv" where under column "old_names" you can list gene names that need to be renamed and under
# column "new_names" list the new hybrid gene name. See https://github.com/PoolLab/generecovery/tree/main/mouse_mm10_input_files for formating example.

################################################################
#### B. Recover intergenic reads from 3' un-annotated exons ####
################################################################


#### Instantiate libraries and define inputs

library("GenomicAlignments")
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
# A. gene_extension_candidates.csv

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

bamfile = file.path("<location of Cell Ranger generated bam file with aligned sequencing reads>/possorted_genome_bam.bam") # Add the location of the Cell Ranger generated bam file containing transcriptome aligned reads with the unoptimized reference
indexfile = file.path("<location of Cell Ranger generated bam file with aligned sequencing reads>/possorted_genome_bam.bam") #Note, you don't have to specify ".bai" extension here.
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
export.bed(m, con = "./intergenic_reads.bed")


#### B2: Create gene ranges file for linking reads to genes ####
################################################################

#### Purpose: Make a bed file with gene boundaries, which is required for assigning intergenic reads to a specific gene. This step requires R to be able to access a linux terminal and a software package that is only available on linux / mac OS.

#### Input: Genome annotation file ("xxx.gtf" file, from 10x Genomics provided reference transcriptome "gene" folder: can be downloaded at "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest" or Ensembl.org if wish to customize more)

#### Notes:
# This step requires BEDOPS (https://bedops.readthedocs.io/en/latest/index.html) which is only available for linux based environments. which is only available for linux environments.
# We will invoke the bash scripts from R for which you need to have BEDOPS installed on your linux/mac OS system.
# Alternatively, one can run the BEDOPS and relevant bash scripts straight from a linux terminal if desired. In this case, make sure BEDOPS is set in PATH.


# Extract all gene entries from genome annotation. Assuming that input unoptimzied genome reference is named "genes.gtf" per 10x Genomics convention.

gene_ranges_df<- import(con = "./genes.gtf", format = "gtf") # Import the original unoptimized genome annotation file
gene_ranges_df = as.data.frame(gene_ranges_df)
gene_ranges_df = gene_ranges_df[gene_ranges_df$type == "gene",] # Extract all "gene" entries in the genome annotation to a new variable
gene_ranges_df = makeGRangesFromDataFrame(gene_ranges_df, keep.extra.columns=TRUE)
rtracklayer::export(gene_ranges_df, "gene_ranges.gtf", format = "gtf")

## Add "transcript_id """ column to the gtf file to make it compatible with BEDOPS expectation of bed format

system('awk \'{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }\' gene_ranges.gtf > gene_ranges1.gtf')

## Convert reference gtf into bed with BEDOPS command gtf2bed. Note, that R needs to have access to linux terminal and uses BEDOPS software that is only available for linux / mac OS.
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "/usr/bin/bedops", sep = ":")) # modify BEDOPS location folder if different
system("gtf2bed < gene_ranges1.gtf > gene_ranges.bed") # Creates a bed file with gene boundaries

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

file.remove("./gene_ranges.gtf")
file.remove("./gene_ranges1.gtf")

## In R: save outcome
write.table(gene_ranges, "gene_ranges.bed", sep="\t",row.names=FALSE, col.names=FALSE, quote = FALSE)

#### B3: Identify candidate genes for extension with excess 3' intergenic reads ####
####################################################################################

# Sort bed files
system("sort -k 1,1 -k2,2n gene_ranges.bed > gene_ranges_sorted.bed")
system("sort -k 1,1 -k2,2n intergenic_reads.bed > intergenic_reads_sorted.bed")

## This step requires bedtools (https://bedtools.readthedocs.io/en/latest/content/installation.html) which is only available for linux environments. We will invoke the bash scripts from R for which you need to have bedtools installed on your linux/mac system.
# Alternatively, one can run the bedtools commands and relevant bash scripts straight from a linux terminal. In this case, make sure bedtools is set in PATH.

old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "/opt/bedtools2/bin", sep = ":")) # Add bedtools to PATH. Modify bedtools location folder if different

system("bedtools closest -a intergenic_reads_sorted.bed -b gene_ranges_sorted.bed -s -D a -fu > results.txt") # resulting file contains sequencing reads with distance data from closest 3' gene identity and end

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

colnames(summary_data_genes) =  c("genes", "3primcount")
write.csv(summary_data_genes, "gene_extension_candidates.csv") # Saves a rank ordered list of genes as a function of 3' intergenic read mapping within 10kb of known gene end. You can use this as a prioritized gene list for gene extension to examine in Integrated Genomics Viewer.

summary_data_genes["<gene_of_interest"] # display 3' intergenic read mapping for a gene of interest



#### B4: Add Ensembl 3' gene ends to gene extension decision file ####
##############################################################################################

#### Acquire ensembl gene boundaries from 10x Genomics/Ensembl genome annotation. The genome annotations can be acquired from 10x Genomics provided reference transcriptome "gene" folder: can be downloaded at "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest"

ensembl_db <- import(con ="./genes.gtf", format = "gtf")

ensembl_db = as.data.frame(ensembl_db)
ensembl_db = ensembl_db[,c(13,2,3,5,7)] # Extract relevant information from genome annotation (gene, start, end, strand and feature)
colnames(ensembl_db) = c("gene", "start", "end", "strand", "type")
ensembl_db = ensembl_db[ensembl_db$type == "gene",] # Keep only "gene" features from the annotation 
head(ensembl_db)
dim(ensembl_db)

## Deal with gene duplications in ensembl genome annotation (weirdly the annotation contains several entries for some genes)

ensembl_duplicates=unique(ensembl_db$gene[duplicated(ensembl_db$gene)])
start = rep(0, length(ensembl_duplicates))
end = rep(0, length(ensembl_duplicates))
strand = as.factor(rep(c("+", "-"), length(ensembl_duplicates)/2))
type = rep("gene", length(ensembl_duplicates))

ensembl_fix = data.frame(cbind(ensembl_duplicates, start, end, strand, type))
ensembl_fix$strand = strand
colnames(ensembl_fix) = c("gene", "start", "end", "strand", "type")
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

#### OPTIONAL: Get Refseq gene boundaries (for many genes Refseq has annotated 3' gene boundaries to extend significantly further than in the Ensembl/Gencode references. This allows for rapid extension of hundreds of genes). Refseq genome annotations can be accessed from UCSC genome browser: http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/ for mice and https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/ where you can download hg38.ncbiRefSeq.gtf.gz for humans.

refseq_db = import(con ="C:/Users/allan/OneDrive/Desktop/Helen Allan Computational Paper/Data/Refseq_GTFs/mm10.ncbiRefSeq.gtf", format = "gtf")
refseq_db = as.data.frame(refseq_db)
refseq_db = refseq_db[,c(12,2,3,5,7)] # Note, these positions can change with successive updates of the format.
colnames(refseq_db) = c("gene", "start", "end", "strand", "type")
refseq_db = refseq_db[refseq_db$feature == "transcript",]
dim(refseq_db)
head(refseq_db)


#### OPTIONAL: Deal with gene duplications in refseq genome annotation (similarly to ensembl reference, the Refseq reference somewhat unexpectedly contains several gene entries for some genes)

refseq_duplicates=unique(refseq_db$gene[duplicated(refseq_db$gene)])
start = rep(0, length(refseq_duplicates))
end = rep(0, length(refseq_duplicates))
strand = as.factor(rep(c("+", "-"), length(refseq_duplicates)/2))
feature = rep("gene", length(refseq_duplicates))

refseq_fix = data.frame(cbind(refseq_duplicates, start, end, strand, type))
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

strand = as.factor(rep(c("+", "-"), 10000))[1:nrow(summary_data_genes)]
ensembl_start = rep(0, nrow(summary_data_genes))
ensembl_end = rep(0, nrow(summary_data_genes))
# refseq_term_left = rep(0, nrow(summary_data_genes))
# refseq_term_right = rep(0, nrow(summary_data_genes))
# ref_seq_overhang = rep(0, nrow(summary_data_genes))
update_start = as.numeric(rep(NA, nrow(summary_data_genes)))
update_end = as.numeric(rep(NA, nrow(summary_data_genes)))

for (i in 1:dim(summary_data_genes)[1]){
  gene = summary_data_genes$genes[i]
  ensembl_start[i] = as.numeric(ensembl_db[gene,]$start)
  ensembl_end[i] = as.numeric(ensembl_db[gene,]$end)
  strand[i] = ensembl_db[gene,]$strand
}

summary_data_genes = data.frame(cbind(summary_data_genes, strand, ensembl_start, ensembl_end, update_start, update_end))
# summary_data_genes = data.frame(cbind(summary_data_genes, strand, ensembl_start, ensembl_end, refseq_term_left, refseq_term_right, ref_seq_overhang)) # if refseq data is used


# ref_seq_genes = intersect(rownames(summary_data_genes), row.names(refseq_db))
#for (i in ref_seq_genes){
#  if (summary_data_genes[i,2] == "+"){
#    summary_data_genes[i, 6] = as.numeric(refseq_db[i, 2])
#    summary_data_genes[i, 7] = summary_data_genes[i, 6] - summary_data_genes[i, 4]
#  } else {
#    summary_data_genes[i, 5] = as.numeric(refseq_db[i, 1])
#    summary_data_genes[i, 7] = summary_data_genes[i, 3] - summary_data_genes[i, 5]
#  }
#}

write.csv(summary_data_genes, "gene_extension_candidates.csv") # Final decision file for gene extension curation
file.remove("./gene_ranges.bed")
file.remove("./gene_ranges.gtf")
file.remove("./gene_ranges1.gtf")
file.remove("./gene_ranges_sorted.bed")
file.remove("./intergenic_reads.bed")
file.remove("./intergenic_reads_sorted.bed")
file.remove("./results.txt")

#### B5: Manual curation of candidate genes for 3' gene extension ####
######################################################################

# Download and install https://software.broadinstitute.org/software/igv/
# Load the transcriptome aligned sequencing read data ("xxx.bam" file) to the igv and choose appropriate species/genome.
# Go through the rank ordered candidate gene list in "gene_extension_candidates.csv" and examine evidence for unannotated 3' exons and UTRs:
#  a) extensive splicing between annotated and candidate unannotated regions at the genetic locus as visualized by igv.
#  b) continuous read mapping from known gene end (delayed cutoff of sequencing reads posterior from known gene end)
#  c) examine external evidence (e.g. Allen in situ brain atlas, Human Protein Atlas etc.)

# Based on the latter, provide new gene boundaries under two new columns: "update_start" and "update_end" depending on whether gene is on "-" or "+" strands, respectively. Use excel formulas, if necessary, to scale curation, and save file as "gene_extension_candidates.xlsx".
# Resulting manually curated "gene_extension_candidates.csv" will be used as input in the Optimized annotation assembly stage.
