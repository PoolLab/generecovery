#### Analysis of human T-lymphocytes from PBMC data mapped to exonic scRNA-seq reference

PBMC.data = read.csv("./T_lymphocytes_exonic_mapping.csv", header = T, row.names = 1) # change folder if data not saved in the same location
PBMC_lymph_dbase <- CreateSeuratObject(PBMC.data, min.features = 0, min.cells = 0, project = "PBMC_lymph")

# Plot QC data and probe basic data metrics

VlnPlot(object = PBMC_lymph_dbase, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"))
length(Cells(PBMC_lymph_dbase)) # number of cells
median(PBMC_lymph_dbase$nFeature_RNA) # Median resoultion in genes/cell

# Data normalization

PBMC_lymph_dbase <- NormalizeData(object = PBMC_lymph_dbase, normalization.method = "LogNormalize", scale.factor = 10000)

# Detect variable genes
PBMC_lymph_dbase <- FindVariableFeatures(object = PBMC_lymph_dbase, selection.method = "vst", nfeatures = 450, loess.span = 0.3, clip.max = "auto")

# Scale the data
PBMC_lymph_dbase <- ScaleData(object = PBMC_lymph_dbase, do.scale = TRUE, do.center = TRUE) # vars.to.regress = c("percent.mito", "nCount_RNA")
# Carry out PCA and evaluate PC dimensions

PBMC_lymph_dbase <- RunPCA(object = PBMC_lymph_dbase, features = VariableFeatures(object=PBMC_lymph_dbase), verbose = TRUE, ndims.print = 1:5, nfeatures.print = 5, npcs=40, seed.use = 42)
ElbowPlot(object = PBMC_lymph_dbase, reduction = "pca", ndims = 40)

# Find clusters (save.SNN = T saves the SNN so that the clustering algorithm can be rerun using the same graph but with a different resolution value (see docs for full details)

PBMC_lymph_dbase <- FindNeighbors(object = PBMC_lymph_dbase, dims = 1:11, reduction = "pca", compute.SNN = TRUE, verbose = TRUE)
PBMC_lymph_dbase <- FindClusters(object = PBMC_lymph_dbase, resolution = 1.2, modularity.fxn = 1, algorithm = 1, n.start = 10, n.iter = 10, random.seed = 1, verbose = TRUE)

# Run Non-linear dimensional reduction (tSNE)
PBMC_lymph_dbase <- RunTSNE(object = PBMC_lymph_dbase, dims = 1:11, reduction = "pca", tsne.method = "Rtsne", seed.use = 8, dim.embed = 2, perplexity = 25)
DimPlot(object = PBMC_lymph_dbase, reduction = "tsne", label = TRUE, pt.size = 1.5, label.size = 3) + NoLegend() # group.by = "orig.ident",  no.legend = F, label.size = 5
