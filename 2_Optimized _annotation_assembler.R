###############################################
#### 2_Automated genome reference assembly ####
###############################################


#### Input data and statement of purpose ####
#############################################

# 0. Load data and libraries
#    - genome annotation file (.gtf) to be optimized
#    - "overlapping_gene_list.csv" file specifying how to resolve gene overlap derived issues. "Delete" entries in $final_classification field mark genes for deletion. Transcript names in $transcripts_for_deletion mark specific transcripts for deletion.
#    - "gene_extension_candidates.csv" specifying updated gene boundaries for incorporating intergenic reads
#    - "rename_genes.csv" specifying gene names to be replaced and new names (under $old_names and $new_names fields, respectively)
# 1. Creates pre-mRNA genome annotation from input genome annotation. This step extracts all transcript entries from the genome annotation and defines them as full length exons with new transcript IDs and corresponding transcripts. This allows to capture many intronically mapped reads that otherwise get discarded.
# 2. Gene deletion step: Deletes all annotation entries for genes destined for deletion (has "Delete" entry in $final_classification field of "overlapping_gene_list.csv"
# 3. Transcript deletion step: Deletes all transcripts destined for deletion (transcript names listed in the "transcripts_for_deletion" column in ""overlapping_gene_list.csv"
# 4. Gene coordinate adjustment step: replace the left most or right most coordinate of the first exon of a gene in genome annotation if there is a coordinate in columns $new_left or $new_right in the "gene_extension_candidates.csv".
# 5. Add pre-mRNA reads to all genes not in the gene overlap list.
# 6. Rename genes to avoid discarding expression data with near perfect terminal exon overlap.
# 7. Save the optimized genome annotation in a new gtf file


#### 0. Load libraries and import data ####
###########################################
install.packages("BiocManager")
install.packages("gdata")
install.packages("stringr")
BiocManager::install("rtracklayer")
BiocManager::install("GenomicRanges")


library("rtracklayer")
library("gdata")
library("stringr")
library("rtracklayer")
library("GenomicRanges")

exonic_gtf<- import(con = 'genes.gtf', format = "gtf")
exonic_df = as.data.frame(exonic_gtf)

overlap_df = read.csv("overlapping_gene_list.csv", header=T)

boundary_fix = read.csv("gene_extension_candidates.csv", header=T)

rename_genes = read.csv("rename_genes.csv", header=T)

new_df = exonic_df
rm(exonic_gtf)

####  1. Create premRNA genome annotation from input gtf that defines transcripts as exons ####
###############################################################################################

# Note: several methods are available for recovering intronic reads from scRNA-seq data. These
# include the early approach of redefining all genes/transcripts as exons, which enables
# incorporation of intronic reads to downstream analysis. This approach works well but has
# several problems including poorer performance at mapping spliced reads (annotation loses
# splicing boundaries) and eliminates hundreds of genes due to gene overlaps. There are also
# specific modes in the latest versions of sequence aligners (e.g. --include-introns parameter
# in Cell Ranger 6) for incorporating intronic reads that work better. However, we found
# that the latter systematically fails to retrieve intronic reads for many genes that the
# traditional premrna approach performs well on. We therefore adopted a hybrid pre-mRNA
# reference approach where we supplement normal gene annotation entries by traditional
# pre-mRNA entries where transcripts have been redefined as exons and map in the
# --include-introns mode in Cell Ranger 6 or default mode in later ierations to retrieve 
# most of available intronic reads.

transcripts_df = exonic_df[exonic_df$type == "transcript",]
exons_df = transcripts_df # Create new dataframe to contain premrna exons
exons_df$type = rep("exon", nrow(exons_df)) # rename "type" from transcripts to exon

premrna_df = interleave(transcripts_df, exons_df) # interleave transript entries with exon entries
premrna_df$transcript_id = gsub("0000", "8888", premrna_df$transcript_id)

####  2. Delete select genes ####
#################################

genes_to_delete = overlap_df$genes[overlap_df$final_classification == "Delete"]
new_df = new_df[!new_df$gene_name %in% genes_to_delete,]

####  3. Delete select transcripts ####
#######################################

transcripts_to_delete = overlap_df$transcripts_for_deletion
transcripts_to_delete <- transcripts_to_delete[transcripts_to_delete!=""]

transcripts_to_delete_final = transcripts_to_delete[!str_detect(transcripts_to_delete, ", ")]

for (i in 1:length(transcripts_to_delete)){
  a = transcripts_to_delete[i]
  if (str_detect(a, ", ")){
    split_elements <- unlist(str_split(a, ", "))
    transcripts_to_delete_final = c(transcripts_to_delete_final, split_elements)
  }
}

transcripts_to_delete = transcripts_to_delete_final

new_df = new_df[!new_df$transcript_name %in% transcripts_to_delete,]


####  4. Adjust gene coordinates ####
#####################################

left_genes = as.data.frame(cbind(boundary_fix$genes[!is.na(boundary_fix$update_start)], boundary_fix$update_start[!is.na(boundary_fix$update_start)]))
colnames(left_genes) = c("genes", "update_start")
left_genes$update_start = as.numeric(left_genes$update_start)
right_genes = as.data.frame(cbind(boundary_fix$genes[!is.na(boundary_fix$update_end)], boundary_fix$update_end[!is.na(boundary_fix$update_end)]))
colnames(right_genes) = c("genes", "update_end")
right_genes$update_end = as.numeric(right_genes$update_end)

left_exon_difs = rep(0, length(left_genes)) # for troubleshooting
right_exon_difs = rep(0, length(right_genes))

for (i in 1:dim(left_genes)[1]){
  gene_entries = which(new_df$gene_name == left_genes[i, 1])
  type_entries = new_df$type[gene_entries]
  first_gene_exon = head(gene_entries[type_entries == "exon"], 1)
  new_df[first_gene_exon, 2] = left_genes[i, 2]
  left_exon_difs[i] = new_df[first_gene_exon, 3] - new_df[first_gene_exon, 2]
}

for (i in 1:dim(right_genes)[1]){
  gene_entries = which(new_df$gene_name == right_genes[i, 1])
  type_entries = new_df$type[gene_entries]
  last_gene_exon = tail(gene_entries[type_entries == "exon"], 1)
  new_df[last_gene_exon, 3] = right_genes[i, 2]
  right_exon_difs[i] = new_df[last_gene_exon, 3] - new_df[last_gene_exon, 2]
}

#### 5. Add pre-mRNA transcripts to genes not in the gene overlap list ####
############################################################################

# Explanation: Cellranger --include-introns mode unfortunately does not pick up on many intronic reads (unclear why despite lengthy correspondence with their support). I can pick those up however if I add the pre-mRNA transcripts to respective genes as exons with new transcript_id values.

## Genes to modify

# Note that all genes in the overlap gene list (overlap_df$genes) will be excluded from this step to avoid the introduction of large exonic overlaps between genes

genes_to_append = unique(new_df$gene_name)
genes_to_append = setdiff(genes_to_append, overlap_df$genes)

## Reformat the final gtf dataframe such that we can add premrna data to it

new_df = new_df[, colnames(premrna_df)]

## Append premrna transcript to the end of the gene

genes_to_append = genes_to_append[1:(length(genes_to_append)-1)]

for (i in genes_to_append){
  insert = premrna_df[premrna_df$gene_name %in% i,]
  first_section = new_df[0:tail(which(new_df$gene_name == i), 1),]
  last_section = new_df[(tail(which(new_df$gene_name == i), 1)+1):dim(new_df)[1],]
  new_df = rbind(first_section, insert, last_section)
}

#### 6. Rename desired genes ####
#################################

# Rename desired genes (example from mouse genome): "Cers1"==>"Cers1_Gdf1" // "Chtf8" ==> "Chtf8_Derpc" // "Insl3" ==> "Insl3_Jak3" // "Pcdhga1" ==> "Pcdhg_all" // "Pcdha1" ==> "Pcdha_all" // "Ugt1a10" ==> "Ugt1a_all" // "4933427D14Rik" ==> "4933427D14Rik_Gm43951" // "Mkks" ==> "Mkks_plus"

old_names = rename_genes$old_names
new_names = rename_genes$new_names

for (i in 1:length(old_names)){
  new_df$transcript_name = str_replace_all(new_df$transcript_name, old_names[i], new_names[i])
  new_df$gene_name = str_replace_all(new_df$gene_name, old_names[i], new_names[i])
} 

sum(str_detect(new_df$gene_name, "Cers1"))
sum(str_detect(new_df$gene_name, "Cers1-Gdf1"))

sum(str_detect(new_df$transcript_name[!is.na(new_df$transcript_name)], "Cers1"))
sum(str_detect(new_df$transcript_name[!is.na(new_df$transcript_name)], "Cers1-Gdf1"))


#### 7. Save the optimized genome annotation in a new gtf file ####
###################################################################

new_gtf = makeGRangesFromDataFrame(new_df, keep.extra.columns=TRUE)

rtracklayer::export(new_gtf, "optimized_reference.gtf", format = "gtf")
