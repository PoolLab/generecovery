####################################################################
#### ScRNA-seq analysis of exonically mapped MnPO neuronal data ####
####################################################################

library("Seurat")

# Load data

MnPO.data = read.csv("./MnPO_neurons_clean_exon.csv", header = T, row.names = 1)
MnPO_neuro_dbase <- CreateSeuratObject(MnPO.data, min.features = 0, min.cells = 0, project = "MnPO_neuro")

# Plot QC data

VlnPlot(object = MnPO_neuro_dbase, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"))

# Print metrics

length(Cells(MnPO_neuro_dbase)) # number of cells
median(MnPO_neuro_dbase$nFeature_RNA) # Median resoultion in genes/cell

# Data normalization

MnPO_neuro_dbase <- NormalizeData(object = MnPO_neuro_dbase, normalization.method = "LogNormalize", scale.factor = 10000)

# Detect variable genes across the single cells

MnPO_neuro_dbase <- FindVariableFeatures(object = MnPO_neuro_dbase, selection.method = "vst", nfeatures = 850, loess.span = 0.3, clip.max = "auto")
VariableFeaturePlot(object = MnPO_neuro_dbase)

# Scale data

MnPO_neuro_dbase <- ScaleData(object = MnPO_neuro_dbase, do.scale = TRUE, do.center = TRUE) # vars.to.regress = c("percent.mito", "nCount_RNA")

#### Carry out PCA and evaluate PC dimensions

MnPO_neuro_dbase <- RunPCA(object = MnPO_neuro_dbase, features = VariableFeatures(object=MnPO_neuro_dbase), verbose = TRUE, ndims.print = 1:5, nfeatures.print = 5, npcs=40, seed.use = 42)

# Determining PCs to use
ElbowPlot(object = MnPO_neuro_dbase, reduction = "pca", ndims = 40)

# Find clusters (save.SNN = T saves the SNN so that the clustering algorithm can be rerun using the same graph but with a different resolution value (see docs for full details)

MnPO_neuro_dbase <- FindNeighbors(object = MnPO_neuro_dbase, reduction = "pca", dims = 1:23, compute.SNN = TRUE, verbose = TRUE)
MnPO_neuro_dbase <- FindClusters(object = MnPO_neuro_dbase, modularity.fxn = 1, resolution = 2, algorithm = 1, n.start = 10, n.iter = 10, random.seed = 1, verbose = TRUE)

#Run Non-linear dimensional reduction (tSNE)
MnPO_neuro_dbase <- RunTSNE(object = MnPO_neuro_dbase, reduction = "pca", dims = 1:23, tsne.method = "Rtsne", seed.use = 1, dim.embed = 2, perplexity = 15)
DimPlot(object = MnPO_neuro_dbase, reduction = "tsne", label = T, pt.size = 1.15) + NoLegend()
