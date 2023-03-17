This is the code repository for the bioRxiv manuscript "Enhanced recovery of single-cell RNA-sequencing reads for missing gene expression data" available at https://www.biorxiv.org/content/10.1101/2022.04.26.489449v1. If you are interested in the optimized mouse and human reference transcriptomes for scRNA-seq, you can find the latest versions at https://www.thepoollab.org/. These are optimized versions of the latest mouse and human transcriptomic references provided by 10x Genomics at "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest" (for mouse GRCm38/mm10 and human GRCh38 genome assemblies). You can use the optimized references in lieu of the basic 10x Genomics transcriptomic references by downloading and untaring them (tar -zxvf <xxx.tar.gz file> in terminal).



If you are interested in optimizing your own transcriptomic reference with your own sequencing data, the code for doing that can be found in the following to R files:
* "1_Genome annotation pre-processing.R"
* "2_Optimized _annotation_assembler.R"

We have also generated an R-package "ReferenceEnhancer" for optimizing any transcriptomic reference for scRNA-seq data analysis which is available here: https://github.com/PoolLab/ReferenceEnhancer


Finally, if you are interested in reproducing analysis in this paper, you will find the raw data for the paper in the Gene Expression Omnibus Repository GSE198528 and the code for analysis and figures in the folder "Data_analysis".

For any questions, comments or suggestions please reach out to allan-hermann.pool@utsouthwestern.edu.
