library(Seurat)
library(ggplot2)
library(dplyr)
DefaultAssay(allsample) <- "RNA"
#yanse
cols<-c('firebrick','goldenrod1','seagreen','slateblue')

pal<-colorRampPalette(cols)

#umap
p1 <- DimPlot(allsample, reduction = "umap",label = T,raster=F,raster.dpi = c(1000,1000))+
  theme(legend.position = "right",legend.text = element_text(size =10),
        axis.title.x = element_text(color = 'black',size = 12),
        axis.title.y = element_text(color = 'black',size = 12))+scale_color_manual( values = pal(17))
p2 <- DimPlot(allsample, reduction = "umap",label = TRUE,split.by = "allgroup",raster=TRUE,repel = TRUE,pt.size = 3)+
  scale_color_manual( values = pal(6))

p1<-DimPlot(subset(allsample, sample=="Decidual"),reduction = "umap",label = T,raster = T,group.by = "celltype")+
        theme(legend.position = "right",legend.text = element_text(size =10),
        axis.title.x = element_text(color = 'black',size = 12),
        axis.title.y = element_text(color = 'black',size = 12))+
        scale_color_manual( values = pal(13))
p2<-DimPlot(subset(allsample, sample=="Decidual"),reduction = "tsne",label = T,raster = T,group.by = "celltype")+
  theme(legend.position = "right",legend.text = element_text(size =10),
        axis.title.x = element_text(color = 'black',size = 12),
        axis.title.y = element_text(color = 'black',size = 12))+
  scale_color_manual( values = pal(13))

p3 <- DimPlot(allsample, reduction = "umap",  group.by = "group",cols = c("tan1","orange4"),label = FALSE)+theme(legend.position = "bottom",legend.text = element_text(size = 12),
                                                                                                                 axis.title.x = element_text(color = 'black',size = 12),
                                                                                                                 axis.title.y = element_text(color = 'black',size = 12))
p4 <- DimPlot(allsample, reduction = "umap", group.by = "sample",cols = c("steelblue1","dodgerblue4"),label = FALSE)+theme(legend.position = "bottom",legend.text = element_text(size = 12),
                                                                                                                           axis.title.x = element_text(color = 'black',size = 12),
                                                                                                                           axis.title.y = element_text(color = 'black',size = 12))
CombinePlots(plots = list(p1,p2), ncol =2)
DimPlot(allsample, reduction = "umap", group.by = "celltype")
#Dotplot
markers.to.plot <- c("ECM1","THY1","COL3A1","COL1A1","IGFBP1","PRL","COL4A5","COL14A1","COL4A6","ARC","GEM","BDKRB1","PLIN2",
                     "SPP1","MYH11","ACTA2","CD4","CD3E","CD8A","FOXP3","CD79A","CD19","CD38","CD160","KLRB1","ANXA1","CD39",
                     "ITGA1","CD9","ENTPD1","NCR1","CD86","CSF1R","CD53","IL3RA","CTSS","LYZ","MRC1","CD163","CD209","HBEGF",
                     "CD14","CD52","CD83","ITGAX","HLA-G","JAG1","PECAM1","MMP2","PARP1","GATA3","ITGA6","EGFR","CGB","CSH1",
                     "KRT8","CD34","VEGFA","TGFB1","STMN1","TAC3","PRG2","SERPINE1","ENG","JAM2","HBA1","CD73")
markers.to.plot <- c()

markers.to.plot <- c('CD163','CD14','CSF1R','CD68','IGFBP5','C1R','TAGLN','CCDC80',"ECM1",'CD3D','CD3G','CD3E','CCL5',
                     'S100A8','IL1B','VEGFA','ITGAX','NKG7','KLRC1','GZMB','KLRB1','ACTA2','COL3A1','COL1A1',
                     'CCL21','CAVIN2','GNG11','PECAM1','CDH5','AOC1','PAPPA2','MFAP5','HTRA4','DIO2','HLA-G','XAGE3',
                     'PARP1','SLC27A2','SLC22A11','KRT7',
                     'IGFBP1','IGFBP2','PRL','IGHG1','CD19','IGKC','CD79A','HBA1'
                     ,'HBB','HBM', 'VWF','ACKR1','ID3','ENG','SLCO2A1',
                     'TPSB2','TPSAB1','GATA2','KIT',
                     'CGA','CSH2','CYP19A1','S100P','HSD3B1')
markers.to.plot <- c('CD163','CD14','CSF1R','CD68','CD3D','CD4','CD8A','NKG7','KLRC1','GZMB','S100A8',
                     'PECAM1','IGHG1','CD19','IGKC','CD79A',"SDC1")
markers.to.plot <- c('CD163','CD14','CSF1R','CD68','IGFBP5','C1R','TAGLN','CCDC80',"ECM1",'CD3D','CD3G','CD3E','CCL5',
                     'S100A8','IL1B','VEGFA','ITGAX','NKG7','KLRC1','GZMB','KLRB1','ACTA2','COL3A1','COL1A1',
                     'CCL21','CAVIN2','GNG11','PECAM1','CDH5','AOC1','PAPPA2','MFAP5','HTRA4','DIO2','HLA-G','XAGE3',
                     'PARP1','SLC27A2','SLC22A11','KRT7',
                     'IGFBP1','IGFBP2','PRL','IGHG1','CD19','IGKC','CD79A'
                     ,'ID3','ENG','SLCO2A1',
                     'TPSB2','TPSAB1','GATA2','KIT',
                     'CGA','CSH2','CYP19A1','S100P','HSD3B1')
markers.to.plot <- c("ACO1","IREB2","FTH1","FTL","TF","TFRC","TFR2","SLC11A2","SLC39A8","SLC39A14","SLC7A11","SLC3A2",
                     "SLC40A1","HAMP","GPX4","STEAP4","STEAP3","LPCAT3","ACSL4","AIFM2","HJV","PTGS2","HFE","BLVRB",
                     "ACSF2","NFE2L2","TP53","GSS")
markers.to.plot <- c('TFRC','PTPRC')
DotPlot(subset(allsample), features = rev(c(unique(unlist(PBMC_gene_symbol_的副本)),'PECAM1','KDR')),cols = c("white", "firebrick"), dot.scale = 5) + RotatedAxis()+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) #??ת??????

DotPlot(subset(allsample), features = rev(markers.to.plot ), cols = c("skyblue2", "coral2"), dot.scale = 6) + RotatedAxis()
#UAMPgene
markers.to.plot <- c('CD34','ACKR1', 'PCDH17','LYVE1','CDH5',"ENG","PDPN")
markers.to.plot <- c("PECAM1","KDR","CD34","ENG" ,"ACKR1", "PCDH17","FAM167B", "AQP1")
markers.to.plot <- c('ADM','ANGPT2','CDH5','EPHA2','ENG','EPAS1','FLT4','HSPG2','KDR','MMP2','NOTCH1','SHC1','TIE1','NRP2','STAB1','EGFL7','ROBO4','HSPG2','CTGF')
markers.to.plot <- c('ANXA1','CEBPD','IFI35','JAK1','EIF2AK2','SP100','STAT1','TGFB1','IL32',
                     'ISG15','DDX58','IFIT5','HERC5','TNIP2','RSAD2','CCL2','BIRC3','HES1','ICAM1')
markers.to.plot <- c('SNAI1', 'SNA2','TWIST1','CDH1','CDH2','VIM','ZEB1','ZEB2','CLDN1' )
FeaturePlot(subset(allsample,sample %in% "Placenta"), features = c("MME"),cols = c("gray95","darkblue"),
           ncol = 1,max.cutoff = 4,raster = F,min.cutoff = 0)
#RidgePlot
plots <-  RidgePlot(subset(allsample) ,group.by = "group",
          features = c("VWF") ,
          combine = FALSE)
#Vl
markers.to.plot <- keggSet[["KEGG_PPAR_SIGNALING_PATHWAY"]]
plots <-VlnPlot(subset(allsample,celltype1%in%c('EC','SMC')&group%in% 'Normal'),
                 features = c('TNFRSF10A','TNFRSF10B','TNFRSF10C','TNFRSF10D') ,cols = pal(2),
                 pt.size = 0,ncol = 2,combine = F) #x????ǩƫת45?㣬???½?0.5
 plots <- VlnPlot(subset(allsample) ,split.by = "group",group.by = "celltype1",
                 features = c("ANXA2"),split.plot = TRUE,cols = pal(2),
                 pt.size = 0, combine = TRUE)
plots <- VlnPlot(subset(allsample,group=="Low altitude") ,split.by = "stim",group.by = "celltype",
                 features = c("SLC40A1","SLC6A6","NR1D2","TFRC","RXRB","ESR1","SLC7A2","RXRA",
                              "SLC3A2","SLC29A3","SLC38A2","RARG","HCAR1","SCARB1","PTGER2","NR1D1","SLC44A2","PPARA",
                              "TRPV6","TSPO2"),split.plot = TRUE,cols = pal(2),
                 pt.size = 0, combine = TRUE)

plots <- VlnPlot(subset(allsample,stim%in%"PE"&celltype==c("enEVT") ),group.by = "group",
                 features = c("NFKBIA","NFKBIZ","JUNB","MAFF","IER2","CCL2","CEBPD","JUND",
                              "DUSP1","C2CD4B","CXCL2","GADD45B",'JUN','SOCS3','DNAJB1',"HSPA1B","KLF4","HSPA1A","FAM118A") , pt.size = 0, combine = FALSE)

plots <- VlnPlot(subset(allsample,group%in%"Plain"&celltype %in% c("CTB","enEVT","EVT","STB")) ,group.by = "celltype",
                 features = c("YAP1","PPP2R2B"), split.plot = TRUE,split.by = "stim",
                 pt.size = 0, combine = FALSE)

CombinePlots(plots = plots, ncol =2)
#newid
new.cluster.ids <- c("Macrophages_1(0)","FBs_1(1)","T cells_1(2)","DCs_1(3)","Macrophages_2(4)",
                     "NK_1(5)","DCs_2(6)","SMC_1(7)","T cells_2(8)","Macrophages_3(9)",
                     "T cells_3(10)","Macrophages_4(11)","EC(12)","Macrophages_5(13)","NK_2(14)","EVT(15)","CTB_1(16)",
                     "DSC_1(17)","DSC_2(18)","B cells(19)","T cells_4(20)","CTB_2(21)","Macrophages_5(22)",
                     "Erythrocytes(23)","enEVT(24)","DCs_3(25)","T cells_5(26)","Stem cells(27)","STB(28)",
                     "FBs_2(29)","Macrophages_6(30)","SMC_2(31)")
new.cluster.ids <- c("Macrophages","FBs","T cells","DCs","Macrophages",
                       "NK","DCs","SMC","T cells","Macrophages",
                       "T cells","Macrophages","EC","Macrophages","NK","EVT","CTB",
                       "DSC","DSC","B cells","T cells","CTB","Macrophages",
                       "Erythrocytes","enEVT","DCs","T cells","Stem cells","STB",
                       "FBs","Macrophages","SMC")
new.cluster.ids <- c("Macrophages","T cells","Fibroblasts","Macrophages","Macrophages","DCs","T cells","DCs","NK cells","VSMCs" 
                     ,"T cells","NK cells","Macrophages","Macrophages","Endothelia cells","EVT","CTB","DSC","DSC"
                     ,"B cells","Macrophages","T cells","Erythrocytes","CTB","Undefined(24)","enEVT",
                     "STB","Stem cells")

new.cluster.ids <- c("DSC(0)","Macrophage(1)","Macrophage(2)","DC(3)","DC(4)","T cell(5)","Macrophage(6)","NK(7)","T cell(8)","NK(9)" 
                     ,"SMC(10)","EC(11)","EVT(12)","T cell(13)","CTB(14)","B cell(15)","STB(16)","Macrophage(17)","EC(18)"
                     ,"Stem cell(19)","FB(20)")

new.cluster.ids <- c("DSC","Macrophage","Macrophage","DC","DC","T cell","Macrophage","NK","T cell","NK" 
                     ,"SMC","EC","EVT","T cell","CTB","B cell","STB","Macrophage","EC"
                     ,"Stem cell","FB")
new.cluster.ids <- c()
new.cluster.ids <- c("VSMC_0","VSMC_1","VSMC_2","VSMC_3","VSMC_4","VSMC_5","VSMC_6","VSMC_7")

names(new.cluster.ids) <- levels(allsample)
allsample <- RenameIdents(allsample, new.cluster.ids)
allsample$celltype <- allsample@active.ident
#markers
all.markers <- FindAllMarkers(subset(allsample),only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
markers <- FindMarkers(subset(allsample,group %in% "Low altitude"),ident.1='PE',
                       only.pos = F, min.pct = 0.25, logfc.threshold = 0.5,group.by = "stim")
DoHeatmap(subset(allsample),features =top10$gene ,angle = 0,raster = T,size = 3,hjust = 0.5)
#ratio
metaallsample <- allsample@meta.data
ratio <- metaallsample
ggplot(ratio,aes(x = group,fill=celltype1)) + #x???ķ???Ϊclarity????????ɫΪcolor??J??H??
  geom_bar(position = position_fill()) + 
  scale_fill_manual(values = sample(pal(15),15))+ #??????ɫ?? 
        theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 14),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 14,vjust = -1.2),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 10),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 10,angle = 0,vjust = 0), #x????ǩƫת45?㣬???½?0.5
        panel.grid = element_blank(),
        legend.position = 'top',
        legend.key.height = unit(0.6,'cm'),#????ͼ????ɫ???ĸ߶?
        legend.text = element_text(face = 'bold',color = 'black',size = 8),
        panel.background=element_blank()) + #????????
  labs(y = 'ratio',x='') + #????y????Ϊ??Percent??
  coord_flip() #??ת??????

ggplot(ratio,aes(x =allgroup , fill = celltype)) + #x???ķ???Ϊclarity????????ɫΪcolor??J??H??
  geom_bar(position = position_fill()) + 
  scale_color_manual( values = pal(15))  + #??????ɫ?? #????????
  labs(x='T cell',y = 'frequency')+ #????y????Ϊ??frequency??
  coord_flip() #??ת??????
