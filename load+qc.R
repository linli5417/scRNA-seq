set.seed(4180)
#load
Y13_sc <- Read10X(data.dir = "D:/epcPBMC/Y13_sc/raw_feature_bc_matrix")
Y13_sc <- CreateSeuratObject(counts = Y13_sc, project = "Y13_sc", min.cells =2 ,min.features = 200)
Y13_sc$week <- "24"
Y13_sc$group <- "PE"
#Y13_sc$stim <- "Normal"
Y13_sc <- NormalizeData(Y13_sc, verbose = FALSE) 
Y13_sc <- FindVariableFeatures(Y13_sc, selection.method = "vst", nfeatures = 2000)
#QC
Y13_sc[["percent.mt"]] <- PercentageFeatureSet(Y13_sc, pattern = "^MT-")
HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") # ????ѪҺ??????ϸ??????
HB_m <- match(HB.genes_total,rownames(Y13_sc@assays$RNA))
HB.genes <- rownames(Y13_sc@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
Y13_sc[["percent.HB"]]<-PercentageFeatureSet(Y13_sc,features=HB.genes)
Y13_sc <- subset(Y13_sc,nFeature_RNA > 500 & nFeature_RNA <5000& percent.mt <20 &nCount_RNA<30000)
VlnPlot(Y13_sc,features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.HB"), ncol = 2,pt.size = 0)

#??????ȡϸ??3000??
index <- sample(ncol(Y13_sc),3000,replace=FALSE)
Y13_sc <- Y13_sc[,index]