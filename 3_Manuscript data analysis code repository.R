#####################################################
#### 3. Manuscript data analysis code repository ####
#####################################################


#################################################################################
#### 1. Quantify Exonic, Intronic and Intergenic Read Proportions (Figure 1) ####
#################################################################################

# Input file: Cell Ranger (4.00 or newer) aligned sequencing data (.bam file).

BiocManager::install("GenomicAlignments")
library("GenomicAlignments")

bamfile = file.path("<path_to_data>/possorted_genome_bam.bam")
indexfile = file.path("<path_to_data>/possorted_genome_bam.bam")
import_data = readGAlignments(bamfile, index=indexfile, param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = FALSE), tag = c("GN", "RE", "AN", "xf", "NH"), what = "flag")) # Imports all all aligned reads (except secondary alignments). These include all duplicate reads.
import_data = as.data.frame(import_data)
import_data = import_data[import_data$NH==1, ] # keep only uniquely mapped reads (i.e. single genomic location).
a = table(import_data$RE) # RE tag contains Cell Ranger generated read identity (E=exonic, N=intronic, I=intergenic). Note that Cell Ranger (up to 6.1.2 unless they fix it later), wrongly classifies reads antisense to exons as exonic (10x Genomics confirmed). These false exonic reads can be identified by looking up entries in AN tag (lists genes antisense to read's mapping location).
intronic_reads = a["N"] # final intronic read count

b = table(import_data$RE[is.na(import_data$AN)]) # Classification of all non-antisense reads.
preliminary_intergenic_reads = b["I"] # Intergenic reads that do not include false exonic reads.
corrected_exonic_reads = b["E"] # this is true exonic read count (RE=E and AN=NA) which excludes false exonic reads (RE=E and AN!=NA). 

false_exonic_reads = nrow(import_data[!is.na(import_data$AN),]) # Identifies number of false exonic reads (RE=E and AN!=NA).
corrected_intergenic_reads = preliminary_intergenic_reads + false_exonic_reads

corrected_exonic_reads
intronic_reads
corrected_intergenic_reads

