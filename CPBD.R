#data prepare
counts <- subset(allsample,celltype %in% c("NK","Decidual cell"))
write.table(as.matrix(counts@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
counts$celltype1 <- paste(counts$allgroup,counts$celltype)
meta_data <-counts@meta.data
meta_data <- cbind(rownames(meta_data),counts@meta.data[,'celltype1', drop=F])
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  ϸ???????в?????NA
write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

DotPlot(subset(data,group %in% c("Y18.5","Y13.5")), features = rev(unique(tiedata$genesymbol)),
        cols =c("white","firebrick"),dot.scale = 5,col.max = 5,col.min = -5,assay = "RNA",scale.by = "size") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


setwd('D:/new/HCJ' )
pbmc='D:/new/HCJ/' 
library(psych)
library(qgraph)
library(igraph)
library(dplyr)
library(ggplot2)
library(readxl)
#????ͼall
mypvals <- read.delim(paste0(pbmc,"pvalues.txt"), check.names = FALSE)
mymeans <- read.delim(paste0(pbmc,"means.txt"), check.names = FALSE)

recepterdata <- read_excel(paste0(pbmc,"enEVT30.xlsx"))
recepter <- c("NPR2","NPR3")
interaction.name <-c("EC|SMC","DSC|SMC","FB|SMC","Macrophage|SMC","DC|SMC","T cell|SMC","B cell|SMC","NK|SMC","EVT|SMC"
                     ,"CTB|SMC","STB|SMC","SMC|EC")
interaction.name <-c("LN EVT|LN Macrophages","LP EVT|LP Macrophages","HN EVT|HN Macrophages","HP EVT|HP Macrophages")


myrecepter <-grep("NPR2|NPR3", 
                  mymeans$interacting_pair,value = T)
mymeans %>% dplyr::filter(interacting_pair %in% myrecepter)%>%
  dplyr::select("interacting_pair",starts_with(interaction.name))  %>%  
  reshape2::melt() -> meansdf
colnames(meansdf)<- c("interacting_pair","CC","means")
mypvals %>% dplyr::filter(interacting_pair %in% myrecepter)%>%
  dplyr::select("interacting_pair",starts_with(interaction.name))  %>%  
  reshape2::melt() -> pvalsdf
colnames(pvalsdf)<- c("interacting_pair","CC","pvals")
pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")
pldf$interacting_pair.x <- factor(pldf$interacting_pair.x, levels = rev(myrecepter))
summary((filter(pldf, means >0.1))$means)
library(ggplot2)
library(Seurat)
plot <- pldf %>% filter( means >0.1)
plot <- plot[(1:300),]

plot1 <- ggplot(pldf,aes(CC.x,interacting_pair.x))+
  geom_point(aes(color=log2(means),size=-log10(pvals+0.001)) ) +
  scale_size_continuous(range = c(1,4))+
  scale_color_gradientn(colours = c("black","darkblue","yellow","red") )+ theme_bw()+ 
  theme(panel.grid=element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))#+geom_vline(xintercept = c(7.5,14.5,21.5),lty=5)

plot1+plot2
  