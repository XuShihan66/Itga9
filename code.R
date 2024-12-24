library(Seurat)
library(multtest)
library(dplyr)
library(ggplot2) 
library(patchwork)
library(SeuratData)
library(tidyverse)
library(devtools)
library(jpeg)
library(plot1cell)
library(stringr)
library(RColorBrewer)

KO <- Read10X(data.dir = "KO.matrix")
HE <- Read10X(data.dir = "NC.matrix")
KO_Pd<- Read10X(data.dir = "KO-Pd.matrix")
HE_Pd <- Read10X(data.dir = "HE-Pd.matrix")
HE<- CreateSeuratObject(counts = HE, project = "HE",min.cells =3,min.features=100)
KO<- CreateSeuratObject(counts =KO, project = "KO",min.cells =3,min.features=100)
HE_Pd<- CreateSeuratObject(counts = HE_Pd, project = "HE_Pd",min.cells =3,min.features=100)
KO_Pd<- CreateSeuratObject(counts =KO_Pd, project = "KO_Pd",min.cells =3,min.features=100)
pbmc<-merge(x=HE,y=c(KO,HE_Pd,KO_Pd),add.cell.ids=c("HE","KO","HE_Pd","KO_Pd"),project="ALL")

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, 
                                             pattern = "^mt-")
pbmc[["percent.hb"]] <- PercentageFeatureSet(pbmc, 
                                             pattern = "^Hb[^(p)]")
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15& percent.hb < 1)

ifnb.list <- SplitObject(pbmc, split.by = "orig.ident")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 4000)
})

features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"


immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident")
p1 

combined <- FindClusters(immune.combined, resolution = 0.5)
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10<-markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file='top10.csv')

Idents(combined)<-"seurat_clusters"
current_ident<- c('0','1',"2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")
new_ident<- c('Epithelial','Epithelial',"Epithelial","Epithelial",'Epithelial','Epithelial',"Epithelial","Epithelial","Immune","Endothelial","Epithelial","Fibroblast","Immune","Immune","Immune","Epithelial","Fibroblast","Epithelial","Fibroblast","Endothelial","Immune","Fibroblast")
combined@active.ident <- plyr::mapvalues(x = combined@active.ident, from = current_ident, to = new_ident)
combined[["new_ORIG"]] <- Idents(object = combined)
table(combined$orig.ident,combined$new_ORIG)##Data for Fig.2a

##Fig.1e,f
DimPlot(combined,reduction="umap",label = TRUE)
saveRDS(combined,file="combined_all.rds")
DefaultAssay(combined) <- "RNA"
FeaturePlot(combined,features = c("Krt14","Krt10","Col1a1","Pecam1","Ptprc","Cd14"),cols = c("grey","blue"),label =F,ncol=3)


immune<-subset(combined,new_ORIG=="Immune")
DefaultAssay(immune) <- "integrated"
immune <- FindClusters(immune, resolution = 0.5)
DimPlot(immune,reduction="umap",label = TRUE)
markers <- FindAllMarkers(immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10<-markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file='top10-immune.csv')

Idents(immune)<-"seurat_clusters"
current_ident<- c('0','1',"2","3","4","5","6","7","8","9","10")
new_ident<- c('MDSCs','Lymphocytes',"Macrophages","Neutrophils","Monocytes","Proliferating monocytes",
              "Neutrophils","Dendritic cells","Lymphocytes",'MDSCs','MDSCs')
immune@active.ident <- plyr::mapvalues(x = immune@active.ident, from = current_ident, to = new_ident)
immune[["new_ORIG"]] <- Idents(object = immune)
table(immune$new_ORIG)
table(immune$orig.ident)
saveRDS(immune,file="All_immune.rds")

##Fig.2b,c
DimPlot(immune,reduction="umap",label = F,split.by = "orig.ident",pt.size=0.1)

Idents(immune)<-"new_ORIG"
list_genes=list(T=c("Cd3e","Cd3g","Cd8a"),
                NK=c("Nkg7","Gzmb", "Cst7"),
                B=c("Cd79a", "Cd79b", "Bank1"),
                proliferatingmonocytes=c("Mki67","Pcna","Kif4","Cd14"),
                Neutrophil=c('Ly6g',"S100a8", "Sod2","G0s2" ),
                DC=c("Cst3","H2-DMa","Ccl22"),
                Macrophage=c("C1qa","C1qb","Selenop","Adgre1"),
                Monocytes=c("Ly6c2","Ccr2","Lyz2"))


p1=DotPlot(immune,
           features=list_genes,
           cols = c(c('#330066','#336699','#66CC66',"#FFCC33","red")),
           cluster.idents = T)+
  RotatedAxis()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5),
        panel.border = element_rect(color="black"), 
        panel.spacing = unit(1, "mm"))+ 
 
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ 
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66',"#FFCC33",'red')) 
p1

table(immune$orig.ident,immune$new_ORIG)##Data for Fig.2d


table(immune@meta.data$orig.ident)
patients_metadata <- data.table::fread("group.csv", header = TRUE)#"group.csv" has been uploaded to GitHub
View(patients_metadata)

metadata <- FetchData(immune,"orig.ident")
metadata$cell_id <- rownames(metadata)
View(metadata)
metadata <- left_join(x= metadata, y= patients_metadata, by = "orig.ident")
rownames(metadata)<-metadata$cell_id
immune<-AddMetaData(immune,metadata =metadata)
table(immune@meta.data$group)

#Fig.2e
DefaultAssay(immune) <- "RNA"
Idents(immune)<-"new_ORIG"
VlnPlot(immune, 
        features = c("Itga9"),split.by="group", pt.size = 0.1)

##Fig.S1
Idents(immune)<-"new_ORIG"
list_genes=c("S100a8","S100a9","Il1b", "Arg1", "Nos2", "Csf1","Tgfb1",
             "Ccl3", "Ccl4", "Cxcl1", "Cxcl2","Cxcl3") 
DotPlot(immune,features = list_genes,cols = c("black","red"),
        assay = 'RNA')+RotatedAxis()


##Fig.4a
rm(list = ls())
combined<-readRDS("combined_all.rds")
immune<-readRDS("All_immune.rds")
iri.integrated <-combined
colnames(iri.integrated@meta.data)  
circ_data <- prepare_circlize_data(iri.integrated, scale = 0.8 )
set.seed(123)

cluster_colors<-rand_color(length(levels(iri.integrated)))
rep_colors<-rand_color(length(names(table(iri.integrated$orig.ident))))
plot_circlize(circ_data,do.label = F, pt.size = 0.1, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 400, repel = T, label.cex = 0.6)
plot_circlize_change <- function (data_plot, do.label = T, contour.levels = c(0.2, 0.3), 
                                  pt.size = 0.5, kde2d.n = 1000, contour.nlevels = 100, bg.color = "#F9F2E4", 
                                  col.use = NULL, label.cex = 0.5, labels.cex = 0.5, circos.cex = 0.5 ,repel = FALSE) 
{
  centers <- data_plot %>% dplyr::group_by(Cluster) %>% summarise(x = median(x = x), 
                                                                  y = median(x = y))
  z <- MASS::kde2d(data_plot$x, data_plot$y, n = kde2d.n)
  celltypes <- names(table(data_plot$Cluster))
  cell_colors <- (scales::hue_pal())(length(celltypes))
  if (!is.null(col.use)) {
    cell_colors = col.use
    col_df <- data.frame(Cluster = celltypes, color2 = col.use)
    cells_order <- rownames(data_plot)
    data_plot <- merge(data_plot, col_df, by = "Cluster")
    rownames(data_plot) <- data_plot$cells
    data_plot <- data_plot[cells_order, ]
    data_plot$Colors <- data_plot$color2
  }
  circos.clear()
  par(bg = bg.color)
  circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0.01, 
                                                            0), track.height = 0.01, gap.degree = c(rep(2, (length(celltypes) - 
                                                                                                              1)), 12), points.overflow.warning = FALSE)
  circos.initialize(sectors = data_plot$Cluster, x = data_plot$x_polar2)
  circos.track(data_plot$Cluster, data_plot$x_polar2, y = data_plot$dim2, 
               bg.border = NA, panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + 
                               mm_y(4), CELL_META$sector.index, cex = labels.cex,
                             col = "black", facing = "bending.inside", niceFacing = T)
                
                 circos.axis(labels.cex = circos.cex, col = "black", labels.col = "black")
               })
  for (i in 1:length(celltypes)) {
    dd <- data_plot[data_plot$Cluster == celltypes[i], ]
    circos.segments(x0 = min(dd$x_polar2), y0 = 0, x1 = max(dd$x_polar2), 
                    y1 = 0, col = cell_colors[i], lwd = 3, sector.index = celltypes[i])
  }
  text(x = 1, y = 0.1, labels = "Cluster", cex = 0.4, col = "black", 
       srt = -90)
  points(data_plot$x, data_plot$y, pch = 19, col = alpha(data_plot$Colors, 
                                                         0.2), cex = pt.size)
  contour(z, drawlabels = F, nlevels = 100, levels = contour.levels, 
          col = "#ae9c76", add = TRUE)
  if (do.label) {
    if (repel) {
      textplot(x = centers$x, y = centers$y, words = centers$Cluster, 
               cex = label.cex, new = F, show.lines = F)
    }
    else {
      text(centers$x, centers$y, labels = centers$Cluster, 
           cex = label.cex, col = "black")
    }
  }
}

plot_circlize_change(circ_data,do.label = T, pt.size = 0.1, 
                     col.use = cluster_colors ,
                     bg.color = 'white', 
                     kde2d.n = 1000, 
                     repel = T, 
                     labels.cex = 0.8, 
                     circos.cex = 0.5,
                     label.cex = 0.8)




add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num = 2) 

Idents(immune)<-"new_ORIG"
iri.immune <-immune

colnames(iri.immune@meta.data)  
circ_data <- prepare_circlize_data(iri.immune, scale = 0.8 )
set.seed(1234)

cluster_colors<-rand_color(length(levels(iri.immune)))
rep_colors<-rand_color(length(names(table(iri.immune$orig.ident))))
plot_circlize(circ_data,do.label = T, pt.size = 0.01, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.6)

plot_circlize_change <- function (data_plot, do.label = T, contour.levels = c(0.2, 0.3), 
                                  pt.size = 0.5, kde2d.n = 1000, contour.nlevels = 100, bg.color = "#F9F2E4", 
                                  col.use = NULL, label.cex = 0.5, labels.cex = 0.5, circos.cex = 0.5 ,repel = FALSE) 
{
  centers <- data_plot %>% dplyr::group_by(Cluster) %>% summarise(x = median(x = x), 
                                                                  y = median(x = y))
  z <- MASS::kde2d(data_plot$x, data_plot$y, n = kde2d.n)
  celltypes <- names(table(data_plot$Cluster))
  cell_colors <- (scales::hue_pal())(length(celltypes))
  if (!is.null(col.use)) {
    cell_colors = col.use
    col_df <- data.frame(Cluster = celltypes, color2 = col.use)
    cells_order <- rownames(data_plot)
    data_plot <- merge(data_plot, col_df, by = "Cluster")
    rownames(data_plot) <- data_plot$cells
    data_plot <- data_plot[cells_order, ]
    data_plot$Colors <- data_plot$color2
  }
  circos.clear()
  par(bg = bg.color)
  circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0.01, 
                                                            0), track.height = 0.01, gap.degree = c(rep(2, (length(celltypes) - 
                                                                                                              1)), 12), points.overflow.warning = FALSE)
  circos.initialize(sectors = data_plot$Cluster, x = data_plot$x_polar2)
  circos.track(data_plot$Cluster, data_plot$x_polar2, y = data_plot$dim2, 
               bg.border = NA, panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + 
                               mm_y(4), CELL_META$sector.index, cex = labels.cex,
                             col = "black", facing = "bending.inside", niceFacing = T)
                 #circos.axis(labels.cex = 0.3, col = "black", labels.col = "black")
                 circos.axis(labels.cex = circos.cex, col = "black", labels.col = "black")
               })
  for (i in 1:length(celltypes)) {
    dd <- data_plot[data_plot$Cluster == celltypes[i], ]
    circos.segments(x0 = min(dd$x_polar2), y0 = 0, x1 = max(dd$x_polar2), 
                    y1 = 0, col = cell_colors[i], lwd = 3, sector.index = celltypes[i])
  }
  text(x = 1, y = 0.1, labels = "Cluster", cex = 0.4, col = "black", 
       srt = -90)
  points(data_plot$x, data_plot$y, pch = 19, col = alpha(data_plot$Colors, 
                                                         0.2), cex = pt.size)
  contour(z, drawlabels = F, nlevels = 100, levels = contour.levels, 
          col = "#ae9c76", add = TRUE)
  if (do.label) {
    if (repel) {
      textplot(x = centers$x, y = centers$y, words = centers$Cluster, 
               cex = label.cex, new = F, show.lines = F)
    }
    else {
      text(centers$x, centers$y, labels = centers$Cluster, 
           cex = label.cex, col = "black")
    }
  }
}

plot_circlize_change(circ_data,do.label = F, pt.size = 0.1, 
                     col.use = cluster_colors ,
                     bg.color = 'white', 
                     kde2d.n = 1000, 
                     repel = T, 
                     labels.cex = 0.01, 
                     circos.cex = 0.01,
                     label.cex = 0.01)




add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num = 2) 

table(immune$orig.ident,immune$new_ORIG)


table(combined$new_ORIG)
table(immune$new_ORIG)
table(immune$orig.ident)
immune$celltype.group <- paste(immune$new_ORIG, immune$orig.ident, sep = "_")
unique(immune$celltype.group)
Idents(immune)<-"celltype.group"

unique(immune$celltype.group)
DefaultAssay(immune) <- "integrated"
Idents(immune)<-"celltype.group"
mydeg <- FindMarkers(immune,ident.1 = "Macrophages_KO",ident.2 = "Macrophages_HE", verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"Macrophages_KOvsHE.csv")
mydeg <- FindMarkers(immune,ident.1 = "Macrophages_KO_Pd",ident.2 = "Macrophages_HE_Pd", verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"Macrophages_KO_PdvsHE_Pd.csv")

mydeg <- FindMarkers(immune,ident.1 = "Macrophages_KO_Pd",ident.2 = "Macrophages_KO", verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"Macrophages_KO_PdvsKONC.csv")
mydeg <- FindMarkers(immune,ident.1 = "Macrophages_HE_Pd",ident.2 = "Macrophages_HE", verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"Macrophages_HE_PdvsHENC.csv")

Macrophages_NC<-read.csv("Macrophages_KOvsHE.csv")
Macrophages_Pd<-read.csv("Macrophages_KO_PdvsHE_Pd.csv")
Macrophages_HE<-read.csv("Macrophages_HE_PdvsHENC.csv")
Macrophages_KO<-read.csv("Macrophages_KO_PdvsKONC.csv")

#Fig.4b
library(ggpubr)
library(ggrepel)

Macrophages_Pd<-read.csv("Macrophages_KO_PdvsHE_Pd.csv")
Macrophages_HE<-read.csv("Macrophages_HE_PdvsHENC.csv")
Macrophages_KO<-read.csv("Macrophages_KO_PdvsKONC.csv")
list <- read.table("inflammatory_genes.txt", header = F, sep = "\t", stringsAsFactors = FALSE)##"inflammatory_genes.csv" has been uploaded to GitHub
View(list)

test=Macrophages_Pd
head(test)[1:6]
class(test)
test$logP <--log10(test$p_val_adj)
test$Group="not-significant"
test$Group[which((test$p_val_adj<0.05)&(test$avg_log2FC> 0.8))]= "up-regulated"
test$Group[which((test$p_val_adj<0.05)&(test$avg_log2FC< -0.8))]= "down-regulated"
table(test$Group)

ggscatter(test,x = "avg_log2FC",y = "logP",xlab="log2FoldChange",ylab="-log10(P.value)",
          alpha=1,size = 3,repel = T,color="Group",
          palette = c("#99FFFF", "#999999", "pink"))+
  geom_hline(yintercept=10,linetype="dashed")+theme_bw()+
  geom_vline(xintercept=c(-0.8,0.8),linetype="dashed")



intersection_values <- intersect(test$X, list$V1)
result <- test[test$X %in% intersection_values, ]
head(result)
result$logP <--log10(result$p_val_adj)
result$Group="not-significant"
result$Group[which((result$p_val_adj<0.05)&(result$avg_log2FC> 0.8))]= "up-regulated"
result$Group[which((result$p_val_adj<0.05)&(result$avg_log2FC< -0.8))]= "down-regulated"
ggscatter(result,x = "avg_log2FC",y = "logP",xlab="log2FoldChange",ylab="-log10(P.value)",
          alpha=1,size = 3,repel = T,color="Group",
          palette = c("blue", "#999999", "red"))+
  geom_hline(yintercept=10,linetype="dashed")+theme_bw()+
  geom_vline(xintercept=c(-0.8,0.8),linetype="dashed")+
  geom_text_repel(
    data = subset(result, result$p_val_adj < 0.05 & abs(result$avg_log2FC) >= 0.8),
    aes(label = X),
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )


test=Macrophages_HE
head(test)[1:6]
class(test)
test$logP <--log10(test$p_val_adj)
test$Group="not-significant"
test$Group[which((test$p_val_adj<0.05)&(test$avg_log2FC> 0.8))]= "up-regulated"
test$Group[which((test$p_val_adj<0.05)&(test$avg_log2FC< -0.8))]= "down-regulated"
table(test$Group)

ggscatter(test,x = "avg_log2FC",y = "logP",xlab="log2FoldChange",ylab="-log10(P.value)",
          alpha=1,size = 3,repel = T,color="Group",
          palette = c("#99FFFF", "#999999", "pink"))+
  geom_hline(yintercept=10,linetype="dashed")+theme_bw()+
  geom_vline(xintercept=c(-0.8,0.8),linetype="dashed")

intersection_values <- intersect(test$X, list$V1)
result <- test[test$X %in% intersection_values, ]
head(result)
result$logP <--log10(result$p_val_adj)
result$Group="not-significant"
result$Group[which((result$p_val_adj<0.05)&(result$avg_log2FC> 0.8))]= "up-regulated"
result$Group[which((result$p_val_adj<0.05)&(result$avg_log2FC< -0.8))]= "down-regulated"
ggscatter(result,x = "avg_log2FC",y = "logP",xlab="log2FoldChange",ylab="-log10(P.value)",
          alpha=1,size = 3,repel = T,color="Group",
          palette = c("blue", "#999999", "red"))+
  geom_hline(yintercept=10,linetype="dashed")+theme_bw()+
  geom_vline(xintercept=c(-0.8,0.8),linetype="dashed")+
  geom_text_repel(
    data = subset(result, result$p_val_adj < 0.05 & abs(result$avg_log2FC) >= 0.8),
    aes(label = X),
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )


#Fig.S3
test=Macrophages_KO
head(test)[1:6]
class(test)
test$logP <--log10(test$p_val_adj)
test$Group="not-significant"
test$Group[which((test$p_val_adj<0.05)&(test$avg_log2FC> 0.8))]= "up-regulated"
test$Group[which((test$p_val_adj<0.05)&(test$avg_log2FC< -0.8))]= "down-regulated"
table(test$Group)


ggscatter(test,x = "avg_log2FC",y = "logP",xlab="log2FoldChange",ylab="-log10(P.value)",
          alpha=1,size = 3,repel = T,color="Group",
          palette = c("#99FFFF", "#999999", "pink"))+
  geom_hline(yintercept=10,linetype="dashed")+theme_bw()+
  geom_vline(xintercept=c(-0.8,0.8),linetype="dashed")

intersection_values <- intersect(test$X, list$V1)
result <- test[test$X %in% intersection_values, ]
head(result)
result$logP <--log10(result$p_val_adj)
result$Group="not-significant"
result$Group[which((result$p_val_adj<0.05)&(result$avg_log2FC> 0.8))]= "up-regulated"
result$Group[which((result$p_val_adj<0.05)&(result$avg_log2FC< -0.8))]= "down-regulated"
ggscatter(result,x = "avg_log2FC",y = "logP",xlab="log2FoldChange",ylab="-log10(P.value)",
          alpha=1,size = 3,repel = T,color="Group",
          palette = c("blue", "#999999", "red"))+
  geom_hline(yintercept=10,linetype="dashed")+theme_bw()+
  geom_vline(xintercept=c(-0.8,0.8),linetype="dashed")+
  geom_text_repel(
    data = subset(result, result$p_val_adj < 0.05 & abs(result$avg_log2FC) >= 0.8),
    aes(label = X),
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )

library(clusterProfiler) 
library(org.Mm.eg.db)
library(data.table)
library(biomaRt)
library(enrichplot)
library(GSEABase)
library(limma)
library(pheatmap)
library(ggsci)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(WGCNA)
library(GSVA)

kegg<-Macrophages_Pd
vector = kegg$p_val_adj < 0.05 & kegg$X !="" 
data_sgni= kegg[vector,]
head(data_sgni)

kegg=as.data.frame(data_sgni)  
deg=kegg
deg$g=ifelse(deg$p_val_adj>0.05,'stable',
             ifelse( deg$avg_log2FC>0.5,'UP',
                     ifelse( deg$avg_log2FC < -0.5,'DOWN','stable') )
)
table(deg$g)

df <- bitr(unique(deg$X), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Mm.eg.db)
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='X')
data_all_sort <- DEG %>% 
  arrange(desc(avg_log2FC))

geneList = data_all_sort$avg_log2FC 
names(geneList) <- data_all_sort$ENTREZID 
head(geneList)

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'mmu',
               nPerm        = 10000,
               minGSSize    = 10,
               maxGSSize    = 200,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none" )
class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
af=as.data.frame(kk2@result)

write.table(af,file="SigMac_PD_kegg.xls",sep="\t",quote=F,col.names=T)



gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Mm.eg.db)
head(gene.df)

##Fig.4c
ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
ego2 <- setReadable(ego2, OrgDb = org.Mm.eg.db)

lineage1_ego <- simplify(ego2, 
                         cutoff=0.5, 
                         by="p.adjust", 
                         select_fun=min)

barplot(lineage1_ego,showCategory=25,font.size =10,label_format = 50)


kegg<-read.csv("Macrophages_KO_PdvsHE_Pd.csv")
vector = kegg$p_val_adj < 0.05 & kegg$X !="" 
data_sgni= kegg[vector,]
head(data_sgni)
kegg=as.data.frame(data_sgni)  
head(kegg)

library(msigdb) 
library(fgsea)   
rt=kegg[order(kegg[,"avg_log2FC"],decreasing=T),]
logFC=rt[,c("X","avg_log2FC")]
list=as.vector(logFC[,"avg_log2FC"])
names(list)=as.vector(logFC[,1])
list[1:10]
gmt=read.gmt("msigdb.v2023.2.Mm.symbols.gmt")
gse=GSEA(list,TERM2GENE = gmt,pvalueCutoff = 1)

kkTab=as.data.frame(gse)
head(kkTab)
write.table(kkTab,file="GSEA.result.xls" ,sep="\t",quote=F,row.names=F)

#Fig.4d
gseaplot2(gse,
          title = "GOBP_BONE_MORPHOGENESIS",  
          "GOBP_BONE_MORPHOGENESIS", 
          color="red", 
          base_size = 10, 
          subplots = 1:8) 

gseaplot2(gse,
          title = "GOBP_POSITIVE_REGULATION_OF_OSTEOBLAST_DIFFERENTIATION",  
          "GOBP_POSITIVE_REGULATION_OF_OSTEOBLAST_DIFFERENTIATION", 
          color="red", 
          base_size = 10, 
          subplots = 1:8) 

##Fig.S4
library(ClusterGVis)
library(org.Mm.eg.db)
library(GOSemSim)

immune<-readRDS("All_immune.rds")
##KO_Pd
KOPD_immune<-subset(immune,orig.ident=="KO_Pd")
table(KOPD_immune$new_ORIG)
Idents(KOPD_immune)<-"new_ORIG"
pbmc.markers.all <- Seurat::FindAllMarkers(KOPD_immune,
                                           only.pos = TRUE,
                                           min.pct = 0.25,
                                           logfc.threshold = 0.25)

pbmc.markers <- pbmc.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)


head(pbmc.markers)

st.data <- prepareDataFromscRNA(object = KOPD_immune,
                                diffData = pbmc.markers.all,
                                showAverage = T)


str(st.data)
saveRDS(st.data, file = "KOPDstdata.rds")


st.data<- lapply(st.data, function(x) {x[is.na(x)] <- 0; return(x)})
str(st.data)


enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Mm.eg.db,
                        type = "BP",
                        organism = "mmu",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)

head(enrich)


markGenes = unique(pbmc.markers$gene)[sample(1:length(unique(pbmc.markers$gene)),40,
                                             replace = F)]


visCluster(object = st.data,
           plot.type = "line")

visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           cluster.order = c(1:7))
dev.off()

visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           
           line.side = "left",
           cluster.order = c(1:7),
           go.col = rep(jjAnno::useMyCol("stallion",n = 7),each = 5),
           add.bar = T)
dev.off()

##HE_Pd
HEPD_immune<-subset(immune,orig.ident=="HE_Pd")
table(HEPD_immune$new_ORIG)
Idents(HEPD_immune)<-"new_ORIG"
pbmc.markers.all <- Seurat::FindAllMarkers(HEPD_immune,
                                           only.pos = TRUE,
                                           min.pct = 0.25,
                                           logfc.threshold = 0.25)


pbmc.markers <- pbmc.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)


head(pbmc.markers)

st.data <- prepareDataFromscRNA(object = HEPD_immune,
                                diffData = pbmc.markers,
                                showAverage = T)

str(st.data)
saveRDS(st.data, file = "HEPDstdata.rds")


st.data<- lapply(st.data, function(x) {x[is.na(x)] <- 0; return(x)})
str(st.data)


enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Mm.eg.db,
                        type = "BP",
                        organism = "mmu",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)

head(enrich)


markGenes = unique(pbmc.markers$gene)[sample(1:length(unique(pbmc.markers$gene)),40,
                                             replace = F)]

visCluster(object = st.data,
           plot.type = "line")

visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           cluster.order = c(1:7))
dev.off()

visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           
           line.side = "left",
           cluster.order = c(1:7),
           go.col = rep(jjAnno::useMyCol("stallion",n = 7),each = 5),
           add.bar = T)
dev.off()

rm(list = ls())
library(Seurat)
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(tidyverse)
library(devtools)
library(clustree)
library(jpeg)
library(CellChat)
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)
library(ggplot2)

combined<-readRDS("combined_all.rds")
combined.list <- SplitObject(combined,split.by = "orig.ident")
HE_Pd <- combined.list$HE_Pd
KO_Pd <- combined.list$KO_Pd
HE <- combined.list$HE
KO <- combined.list$KO

####KO_Pd
pbmc3k.final <- KO_Pd
data.input <- GetAssayData(pbmc3k.final, assay = "RNA", slot = "data")
labels <- Idents(pbmc3k.final)
identity <- data.frame(group = labels, row.names = names(labels)) 
cellchat <- createCellChat(object = data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") 
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"))
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat,features = NULL)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = F)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat.KO_Pd <-cellchat 


####HE_Pd
pbmc3k.final <- HE_Pd
data.input <- GetAssayData(pbmc3k.final, assay = "RNA", slot = "data")
labels <- Idents(pbmc3k.final)
identity <- data.frame(group = labels, row.names = names(labels)) 
cellchat <- createCellChat(object = data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") 
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"))
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat,features = NULL)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = F)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat.HE_Pd <-cellchat 

pbmc3k.final <- KO
data.input <- GetAssayData(pbmc3k.final, assay = "RNA", slot = "data")
labels <- Idents(pbmc3k.final)
identity <- data.frame(group = labels, row.names = names(labels)) 
cellchat <- createCellChat(object = data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") 
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"))
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat,features = NULL)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = F)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat.KO <-cellchat 

####HE_Pd
pbmc3k.final <- HE
data.input <- GetAssayData(pbmc3k.final, assay = "RNA", slot = "data")
labels <- Idents(pbmc3k.final)
identity <- data.frame(group = labels, row.names = names(labels)) 
cellchat <- createCellChat(object = data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") 
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"))
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat,features = NULL)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = F)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat.HE <-cellchat 


levels(cellchat.HE_Pd@idents)
levels(cellchat.KO_Pd@idents)
levels(cellchat.HE@idents)
levels(cellchat.KO@idents)

cellchat.HE_Pd <- netAnalysis_computeCentrality(cellchat.HE_Pd, slot.name = "netP") 
cellchat.HE <- netAnalysis_computeCentrality(cellchat.HE, slot.name = "netP") 
cellchat.KO_Pd <- netAnalysis_computeCentrality(cellchat.KO_Pd, slot.name = "netP") 
cellchat.KO <- netAnalysis_computeCentrality(cellchat.KO, slot.name = "netP") 
saveRDS(cellchat.HE_Pd,file="cellchat.HE_Pd.rds")
saveRDS(cellchat.KO_Pd,file="cellchat.KO_Pd.rds")
saveRDS(cellchat.HE,file="cellchat.HE.rds")
saveRDS(cellchat.KO,file="cellchat.KO.rds")

cellchat.HE_Pd<-readRDS(file="cellchat.HE_Pd.rds")
cellchat.KO_Pd<-readRDS(file="cellchat.KO_Pd.rds")
cellchat.HE<-readRDS(file="cellchat.HE.rds")
cellchat.KO<-readRDS(file="cellchat.KO.rds")

object.list <- list(HE = cellchat.HE, KO = cellchat.KO, HE_Pd = cellchat.HE_Pd, KO_Pd= cellchat.KO_Pd)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

#Fig.5a
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4), measure = "weight")
gg1 + gg2

#Fig.5b
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T,comparison = c(1,3),vertex.label.color = "black",vertex.label.cex = 0.5)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",comparison = c(1,3),vertex.label.color = "black",vertex.label.cex = 0.5)

#Fig.5c
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(3,4))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,comparison = c(3,4))
gg1 + gg2

#Fig.5d
gg1 <- netVisual_heatmap(cellchat,comparison = c(3,4))
gg2 <- netVisual_heatmap(cellchat, measure = "weight",comparison = c(3,4))
gg1 + gg2

#Fig.S5a
gg1 <- netVisual_bubble(cellchat, sources.use = c(2:4), targets.use = 4,  comparison = c(3, 4), max.dataset = 4, title.name = "Increased signaling in KO_Pd", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(2:4), targets.use = 4,  comparison = c(3, 4), max.dataset = 3, title.name = "Decreased signaling in KO_Pd", angle.x = 45, remove.isolate = T)
gg1 + gg2


all<-combined
table(all$new_ORIG)
table(all$orig.ident)
all$celltype.group <- paste(all$new_ORIG, all$orig.ident, sep = "_")
unique(all$celltype.group)
Idents(all)<-"celltype.group"
DefaultAssay(all) <- "integrated"

mydeg <- FindMarkers(all,ident.1 = "Fibroblast_KO_Pd",ident.2 = "Fibroblast_HE_Pd", verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"Fibroblast_KOPdvsHEPd.csv")


library(clusterProfiler) 
library(ggplot2)
library(org.Mm.eg.db)
library(tidyverse)
library(data.table)
library(biomaRt)
library(enrichplot)
library(GSEABase)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(patchwork)
library(WGCNA)
library(GSVA)

kegg<-read.csv("Fibroblast_KOPdvsHEPd.csv")

vector = kegg$p_val_adj < 0.05 & kegg$X !="" 
data_sgni= kegg[vector,]
head(data_sgni)

kegg=as.data.frame(data_sgni)  
head(kegg)

deg=kegg
deg$g=ifelse(deg$p_val_adj>0.05,'stable',
             ifelse( deg$avg_log2FC>0.5,'UP',
                     ifelse( deg$avg_log2FC < -0.5,'DOWN','stable') )
)
table(deg$g)

df <- bitr(unique(deg$X), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Mm.eg.db)
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='X')
data_all_sort <- DEG %>% 
  arrange(desc(avg_log2FC))

geneList = data_all_sort$avg_log2FC 
names(geneList) <- data_all_sort$ENTREZID 
head(geneList)
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'mmu',
               nPerm        = 10000,
               minGSSize    = 10,
               maxGSSize    = 200,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none" )
class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
af=as.data.frame(kk2@result)
write.table(af,file="Fib_pd_kegg.xls",sep="\t",quote=F,col.names=T)

GO_kk_entrez <- gseGO(geneList     = geneList,
                      ont          = "ALL", 
                      OrgDb        = org.Mm.eg.db,
                      keyType      = "ENTREZID",
                      pvalueCutoff = 0.25)   
GO_kk <- DOSE::setReadable(GO_kk_entrez, 
                           OrgDb=org.Mm.eg.db,
                           keyType='ENTREZID')
class(GO_kk)
colnames(GO_kk@result)
GO_result <- as.data.frame(GO_kk)
rownames(GO_kk@result)[head(order(GO_kk@result$enrichmentScore))]
afGo=as.data.frame(GO_kk@result)
write.table(afGo,file="Fib_Pd_go.xls",sep="\t",quote=F,row.names=F)

#Fig.5e
library(aPEAR)
enrichmentNetwork(GO_kk@result,repelLabels = F, drawEllipses = F)
enrichmentNetwork(kk2@result,repelLabels = F, drawEllipses = F)


#Fig.5f
gseaplot2(GO_kk,
          title = "bone morphogenesis",  
          "GO:0060349", 
          color="red", 
          base_size = 20, 
          subplots = 1:3, 
          pvalue_table = T) 
gseaplot2(GO_kk,
          title = "skeletal system development",  
          "GO:0001501", 
          color="red", 
          base_size = 20, 
          subplots = 1:3, 
          pvalue_table = T) 
gseaplot2(GO_kk,
          title = "ossification",  
          "GO:0001503", 
          color="red", 
          base_size = 20, 
          subplots = 1:3, 
          pvalue_table = T) 

#Fig.S5c
kk=af
up_k <- kk[head(order(kk$enrichmentScore,decreasing = T),7),];up_k$group=1
down_k <- kk[tail(order(kk$enrichmentScore,decreasing = T),4),];down_k$group=-1

dat=rbind(up_k,down_k)
colnames(dat)
dat$pvalue = -log10(dat$pvalue)
dat$pvalue=dat$pvalue*dat$group 
dat=dat[order(dat$pvalue,decreasing = F),]
library(ggstatsplot)
p7 <- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group),
             font.size =10,label_format = 20) + 
  geom_bar(stat="identity") + 
  scale_fill_gradient(low="#34bfb5",high="#ff6633",guide = FALSE) + 
  scale_x_discrete(name ="Pathway names",labels = function(x) str_wrap(x, width = 40)) +
  scale_y_continuous(name ="log10P-value") +
  coord_flip() + 
  theme_ggstatsplot()+
  theme(plot.title = element_text(size = 15,hjust = 0.5),  
        axis.text = element_text(size = 10,face = 'bold'),
        panel.grid = element_blank())+
  ggtitle("Pathway Enrichment") 
p7


#Fig.S5b
library(ggpubr)
library(ggrepel)
test<-read.csv("Fibroblasts_KOPdvsHEPd.csv")
head(test)
class(test)

test$logP <--log10(test$p_val_adj)
test$Group="not-significant"
test$Group[which((test$p_val_adj<0.05)&(test$avg_log2FC> 1))]= "up-regulated"
test$Group[which((test$p_val_adj<0.05)&(test$avg_log2FC< -1))]= "down-regulated"

table(test$Group)
ggscatter(test,x = "avg_log2FC",y = "logP",xlab="log2FoldChange",ylab="-log10(P.value)",
          alpha=1,size = 3,repel = T,color="Group", 
          palette = c("#00AFBB", "#999999", "#FC4E07"))+
  geom_hline(yintercept=10,linetype="dashed")+theme_bw()+
  geom_vline(xintercept=c(-1,1),linetype="dashed")+
  geom_text_repel(
    data = subset(test, test$p_val_adj < 0.05 & abs(test$avg_log2FC) >= 1),
    aes(label = X),
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )