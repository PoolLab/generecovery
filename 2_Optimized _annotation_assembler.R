###############################################
#### 2_Automated genome reference assembly ####
###############################################


#### Input data and statement of purpose ####
#############################################

# 0. Load data and libraries
# 1. In mouse_exonic_genes.gtf delete all rows for genes classified as "Delete" in "mouse_overlap_gene_summary.xlsx" and save as new_genes.gtf
# 2. In new_genes.gtf delete all rows (transcripts not genes: "transcript_name") containing specific transcripts in the "Transcripts_for_deletion" column in "mouse_overlap_gene_summary.xlsx" and save as new_genes.gtf
# 3. Gene coordinate adjustment step: replace the left most or right most coordinate of the first exon of a gene in new_genes.gtf if there is a coordinate in columns new_left or new_right, respectively and save as new.gtf
# 4. Add pre-mRNA reads to all genes not in the gene overlap list (and change the transcript names in the former).
# 5. Rename genes to avoid discarding expression data with near perfect terminal exon overlap.


#### 0. Load libraries and import data ####
###########################################

library("rtracklayer")

exonic_gtf<- import(con = '<location_of_basic_exonic_reference>/genes.gtf', format = "gtf")
premrna_gtf <- import(con = '<location_of_assembled_premRNA_reference>/premrna.gtf', format = "gtf") # Assembled in "1_Annotation pre-processing.R"

exonic_df = as.data.frame(exonic_gtf)
premrna_df = as.data.frame(premrna_gtf)

library(readxl)

overlap_data = read_excel("<location of gene overlap resolving file>/overlapping_gene_list.xlsx")
overlap_df = data.frame(overlap_data)

new_df = exonic_df
rm(exonic_df)
rm(exonic_gtf)
rm(premrna_gtf)

####  1. Delete select genes ####
#################################

genes_to_delete = overlap_df$genes[overlap_df$final_classification == "Delete"]
new_df = new_df[!new_df$gene_name %in% genes_to_delete,]


####  2. Delete select transcripts ####
#######################################

library(stringr)

transcripts_to_delete = overlap_df$transcripts_for_deletion
transcripts_to_delete <- transcripts_to_delete[!is.na(transcripts_to_delete)]

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


####  3. Adjust gene coordinates ####
#####################################

boundary_fix = read_excel("<location of 3'UTR gene extension file>/gene_extension_candidates.xlsx")
boundary_df = data.frame(boundary_fix)

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


#### 4. Add pre-mRNA transcripts to genes not in the gene overlap list ####
############################################################################

# Explanation: Cellranger --include-introns mode unfortunately does not pick up on many intronic reads (unclear why despite lengthy correspondence with their support). I can pick those up however if I add the pre-mRNA transcripts to respective genes as exons with new transcript_id values.

## Genes to modify

overlap_df$genes # genes to exclude

genes_to_append = unique(new_df$gene_name)
genes_to_append = setdiff(genes_to_append, overlap_df$genes)

## Give new transcript_ids to everything in the pre-mRNA gtf

for (i in 1:dim(premrna_df)[1]){
  premrna_df$transcript_id[i] = as.character(i)
}

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

#### 5. Rename desired genes ####
#################################

# Rename desired genes (example from mouse genome): "Cers1"==>"Cers1_Gdf1" // "Chtf8" ==> "Chtf8_Derpc" // "Insl3" ==> "Insl3_Jak3" // "Pcdhga1" ==> "Pcdhg_all" // "Pcdha1" ==> "Pcdha_all" // "Ugt1a10" ==> "Ugt1a_all" // "4933427D14Rik" ==> "4933427D14Rik_Gm43951" // "Mkks" ==> "Mkks_plus"

old_names = c("Cers1", "Chtf8", "Insl3", "Pcdhga1", "Pcdha1", "Ugt1a10", "4933427D14Rik", "Mkks")
new_names = c("Cers1-Gdf1", "Chtf8-Derpc", "Insl3-Jak3", "Pcdhg-all", "Pcdha-all", "Ugt1a-all", "4933427D14Rik-Gm43951", "Mkks-plus")

for (i in 1:length(old_names)){
  new_df$transcript_name = str_replace_all(new_df$transcript_name, old_names[i], new_names[i])
  new_df$gene_name = str_replace_all(new_df$gene_name, old_names[i], new_names[i])
} 

sum(str_detect(new_df$gene_name, "Cers1"))
sum(str_detect(new_df$gene_name, "Cers1-Gdf1"))

sum(str_detect(new_df$transcript_name[!is.na(new_df$transcript_name)], "Cers1"))
sum(str_detect(new_df$transcript_name[!is.na(new_df$transcript_name)], "Cers1-Gdf1"))


#### 6. Export object to gtf file ####
######################################

new_gtf = makeGRangesFromDataFrame(new_df, keep.extra.columns=TRUE)

rtracklayer::export(new_gtf, "optimized_reference.gtf", format = "gtf")

