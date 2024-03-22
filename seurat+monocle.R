set.seed(4180)
#load
sc_data <- Read10X(data.dir = "D:/sc_data/raw_feature_bc_matrix")
sc_data <- CreateSeuratObject(counts = sc_data, project = "sc_data", min.cells =2 ,min.features = 200)
sc_data$group <- "PE"
#sc_data$stim <- "Normal"
sc_data <- NormalizeData(sc_data, verbose = FALSE) 
sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 2000)
#QC
sc_data[["percent.mt"]] <- PercentageFeatureSet(sc_data, pattern = "^MT-")
HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") # ????ѪҺ??????ϸ??????
HB_m <- match(HB.genes_total,rownames(sc_data@assays$RNA))
HB.genes <- rownames(sc_data@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
sc_data[["percent.HB"]]<-PercentageFeatureSet(sc_data,features=HB.genes)
sc_data <- subset(sc_data,nFeature_RNA > 500 & nFeature_RNA <5000& percent.mt <20 &nCount_RNA<30000)
VlnPlot(sc_data,features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.HB"), ncol = 2,pt.size = 0)

#??????ȡϸ??3000??
immune.anchors <- FindIntegrationAnchors(object.list = list()
                                         , dims = 1:20)
gc()
allsample <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
gc()
DefaultAssay(allsample) <- "integrated"
# Run the standard workflow for visualization and clustering
allsample <- ScaleData(allsample, verbose = FALSE)
allsample <- RunPCA(allsample, npcs = 30, verbose = FALSE)
#allsample <- JackStraw(allsample, num.replicate = 100)
#ElbowPlot(allsample)
# t-SNE and Clustering
allsample <- RunUMAP(allsample, reduction = "pca", dims = 1:20)
allsample <- RunTSNE(allsample, reduction = "pca", dims = 1:20)
allsample <- FindNeighbors(allsample, reduction = "pca", dims = 1:20)
allsample <- FindClusters(allsample, resolution = 0.4)
library(Seurat)
library(ggplot2)
library(dplyr)
DefaultAssay(allsample) <- "RNA"
#yanse
cols<-c('firebrick','goldenrod1','seagreen','slateblue')
pal<-colorRampPalette(cols)
DimPlot(allsample, reduction = "umap",label = T,raster=F,raster.dpi = c(1000,1000))+
  theme(legend.position = "right",legend.text = element_text(size =10),
        axis.title.x = element_text(color = 'black',size = 12),
        axis.title.y = element_text(color = 'black',size = 12))+scale_color_manual( values = pal())
#Dotplot
markers.to.plot <- c(markergene$gene)
DotPlot(subset(allsample), features = rev(c(unique(unlist(PBMC_gene_symbol_的副本)),'PECAM1','KDR')),cols = c("white", "firebrick"), dot.scale = 5) + RotatedAxis()+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) #??ת??????
#UAMPgene
FeaturePlot(subset(allsample), features = c(markers.to.plot),cols = c("gray95","darkblue"),
           ncol = 1,max.cutoff = 4,raster = F,min.cutoff = 0)
#Vl
VlnPlot(subset(allsample),
                 features = c(markers.to.plot ) ,cols = pal(2),
                 pt.size = 0,ncol = 2,combine = F)
#newid
new.cluster.ids <- c()
names(new.cluster.ids) <- levels(allsample)
allsample <- RenameIdents(allsample, new.cluster.ids)
allsample$celltype <- allsample@active.ident
#markers
all.markers <- FindAllMarkers(subset(allsample),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
markers <- FindMarkers(subset(allsample),ident.1='PE',
                       only.pos = F, min.pct = 0.25, logfc.threshold = 0.5,group.by = "group")
DoHeatmap(subset(allsample),features =top10$gene ,angle = 0,raster = T,size = 3,hjust = 0.5)
#subset
VSMC <- subset(allsample, celltype %in% 'VSMC')
#############monocle###############
library(monocle3)
library(dplyr)
VSMC <- subset(VSMC)
DefaultAssay(VSMC) <- "RNA"
data<-GetAssayData(VSMC,assay ='RNA',slot ='counts')
cell_metadata <-Tcell@meta.data
gene_annotation <-data.frame(gene_short_name =rownames(Tcell),row.names = rownames(Tcell))
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 50)
## Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, alignment_group = "orig.ident")
## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

## Step 4: Cluster the cells
cds <- cluster_cells(cds)
## Step 5: Learn a graph
cds <- learn_graph(cds)

## Step 6: Order cells
cds <- order_cells(cds)
cols <- c('#AD7B9D',"#92CFD0","#F6CF80","#EE816A",'#D35A7F')
# nolint
pal <- colorRampPalette(cols)
plot_cells(cds,color_cells_by = 'seurat_clusters',label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)+scale_color_manual( values = pal(13))
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=10)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
plot_cells(cds, genes=c(),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
########
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$seurat_clusters)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")
plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(18,43, 14, 11)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
