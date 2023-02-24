
###########################################################################
#### ScRNA-seq analysis of optimal reference mapped MnPO neuronal data ####
###########################################################################

library("Seurat")

# Load cells 

MnPO.data = read.csv(./MnPO_neurons_clean_optim.csv", header = T, row.names = 1) # without doublets
MnPO_neuro_dbase3 <- CreateSeuratObject(MnPO.data, min.features = 0, min.cells = 0, project = "MnPO_neuro")

# Plot QC data

VlnPlot(object = MnPO_neuro_dbase2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"))

# Print metrics

length(Cells(MnPO_neuro_dbase2)) # number of cells
median(MnPO_neuro_dbase2$nFeature_RNA) # Median resoultion in genes/cell


# Data normalization

MnPO_neuro_dbase2 <- NormalizeData(object = MnPO_neuro_dbase2, normalization.method = "LogNormalize", scale.factor = 10000)

# Detect variable genes across the single cells
MnPO_neuro_dbase3 <- FindVariableFeatures(object = MnPO_neuro_dbase2, selection.method = "vst", nfeatures = 850, loess.span = 0.3, clip.max = "auto")
VariableFeaturePlot(object = MnPO_neuro_dbase2)

# Scale the data

MnPO_neuro_dbase2 <- ScaleData(object = MnPO_neuro_dbase2, do.scale = TRUE, do.center = TRUE)

# Carry out PCA

MnPO_neuro_dbase2 <- RunPCA(object = MnPO_neuro_dbase2, features = VariableFeatures(object=MnPO_neuro_dbase2), verbose = TRUE, ndims.print = 1:5, nfeatures.print = 5, npcs=40, seed.use = 42)

# Determining which PCs to use is to look at a plot of the standard deviations of the principle components and draw your cutoff where there is a clear elbow in the graph

ElbowPlot(object = MnPO_neuro_dbase2, reduction = "pca", ndims = 40)

# Find clusters

MnPO_neuro_dbase2 <- FindNeighbors(object = MnPO_neuro_dbase2, reduction = "pca", dims = 1:23, compute.SNN = TRUE, verbose = TRUE)
MnPO_neuro_dbase2 <- FindClusters(object = MnPO_neuro_dbase2, resolution = 2, modularity.fxn = 1, algorithm = 1, n.start = 10, n.iter = 10, random.seed = 1, verbose = TRUE)

#Run Non-linear dimensional reduction (tSNE)
MnPO_neuro_dbase3 <- RunTSNE(object = MnPO_neuro_dbase2, reduction = "pca", dims = 1:23, tsne.method = "Rtsne", seed.use = 4, dim.embed = 2, perplexity = 15)
DimPlot(object = MnPO_neuro_dbase2, reduction = "tsne", label = TRUE, pt.size = 1.15) + NoLegend()
