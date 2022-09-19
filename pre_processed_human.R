# this is an exsample for human PNETs
library(dplyr)
library(Seurat)
library(patchwork)
library(hdf5r)
library(ggplot2)
library(DoubletFinder)
library(scater)
library(Matrix)
library(scRNAseq)
library(SingleR)
library(tibble)
library(scran)
library(pheatmap)
library(celldex)
library(copykat)
library(CellChat)

#STEP 1: read 10x data
pro.name<-"h66110g"
hpn.data <-Read10X_h5(paste0(pro.name,"/filtered_feature_bc_matrix.h5"))

#STEP 2: QC filter
# create Seurat object
hpn <-CreateSeuratObject(counts = hpn.data, project = pro.name, min.cells = 3, min.features = 200)
# calculate the proportion of mitochondrial reads.
mt.genes <- rownames(hpn)[grep("^MT-",rownames(hpn))]
C<-GetAssayData(object = hpn, slot = "counts")
percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
hpn <- AddMetaData(hpn, percent.mito, col.name = "percent.mito")
# calculate the ribosomal proportion
rb.genes <- rownames(hpn)[grep("^RP[SL]",rownames(hpn))]
percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
hpn <- AddMetaData(hpn, percent.ribo, col.name = "percent.ribo")

## plot QC
### plot"nFeature_RNA", "nCount_RNA", "percent.mito","percent.ribo"
pdf(file = paste0(pro.name,'/pdfs/f1.vlplot_features_count_mito_ribo.pdf'), height = 5, width = 10 )
VlnPlot(hpn, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.ribo"), ncol = 4, pt.size = 0)
dev.off()

## plot each two features for them
p1<-FeatureScatter(hpn, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+NoLegend()
p2<-FeatureScatter(hpn, feature1 = "nFeature_RNA", feature2 = "percent.mito")+NoLegend()
p3<-FeatureScatter(hpn, feature1 = "nCount_RNA", feature2 = "percent.mito")+NoLegend()
p4<-FeatureScatter(hpn, feature1="percent.ribo", feature2="nFeature_RNA")+NoLegend()
pdf(file =  paste0(pro.name,'/pdfs/f2.vs_features_count_mito_ribo.pdf'), height = 10, width = 10 )
p1+p2+p3+p4
dev.off()
## read count vs barcode_rank
pdf(file = paste0(pro.name,'/pdfs/f3.vs_features_count_mito_ribo.pdf'), height = 5, width = 5 )
p<-CalculateBarcodeInflections(hpn,barcode.column = "nCount_RNA",group.column = "orig.ident",threshold.low = NULL, threshold.high = NULL)
BarcodeInflectionsPlot(p)+labs(x="barcode_rank")+ NoLegend()
dev.off()
## histogram
pdf(file = paste0(pro.name,'/pdfs/f4.histo_features_count_mito_ribo.pdf'), height = 7, width = 8 )
par(mfrow=c(2,2))
hist(hpn@meta.data$nCount_RNA, breaks =100, xlab = 'nCount_RNA', main = "")
abline(v = 30000, col = "red")
hist(hpn@meta.data$nFeature_RNA, breaks =100, xlab = 'nFeature_RNA',main = "")
abline(v = 7000, col = "red")
abline(v = 200, col = "red")
hist(hpn@meta.data$percent.mito, breaks =100, xlab = 'percent.mito',main = "")
abline(v = 20, col = "red")
hist(hpn@meta.data$percent.ribo, breaks =100, xlab = 'percent.ribo',main = "")
#abline(v = 5, col = "red")
abline(v = 35, col = "red")
dev.off()

#filter by mitochondrial, nfeatures and ncount
hpn.filt <- subset(hpn, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mito < 20 & nCount_RNA<30000 & percent.ribo<40 & percent.ribo>5)
hpn.filt
pdf(file = paste0(pro.name,'/pdfs/f5.vlplot_filtered_features_count_mito_ribo.pdf'), height = 5, width = 10 )
VlnPlot(hpn.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.ribo"), ncol = 4, pt.size = 0)
dev.off()

## convert Seurat to SingleCellExperiment
sce <- as.SingleCellExperiment(hpn.filt)
# Let’s look at what the top 50 expressed genes are. This can be valuable for detecting genes that are overabundant that may be driving a lot of the variation.
pdf(file = paste0(pro.name,'/pdfs/f6.high_exprs.pdf'), height = 10, width = 10 )
plotHighestExprs(sce, exprs_values = "counts")
dev.off()

## remove MT genes and MALAT1
temp<-hpn.filt
counts <- GetAssayData(temp, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('MALAT1'))),]
counts <- counts[-(which(rownames(counts) %in% mt.genes)),]
temp1 <- subset(temp, features = rownames(counts))
hpn.filt<-temp1
temp1<-NULL
sce <- as.SingleCellExperiment(hpn.filt)
pdf(file = paste0(pro.name,'/pdfs/f7.high_exprs.pdf'), height = 10, width = 10 )
plotHighestExprs(sce, exprs_values = "counts")
dev.off()
sce<-NULL

#STEP 3: doubletfinder 0.031 (3.1% for 4000-5000 cells)

seu_kidney<-hpn.filt
seu_kidney <- NormalizeData(seu_kidney)
seu_kidney <- FindVariableFeatures(seu_kidney, selection.method = "vst", nfeatures = 2000)
seu_kidney <- ScaleData(seu_kidney)
seu_kidney <- RunPCA(seu_kidney)
seu_kidney <- RunUMAP(seu_kidney, dims = 1:20)
seu_kidney <- FindNeighbors(seu_kidney, dims = 1:20)
seu_kidney <- FindClusters(seu_kidney, resolution = 0.5)
annotations <- seu_kidney@meta.data$seurat_clusters

pdf(file = paste0(pro.name,'/pdfs/f8.doubletfinder_1.pdf'), height = 10, width = 10 )
DimPlot(seu_kidney, reduction = "umap", label = TRUE)
dev.off()

#pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:20, sct = FALSE)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.016*nrow(seu_kidney@meta.data))  ## Assuming 1.6% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_1560", sct = FALSE)

colnames(seu_kidney@meta.data)[9]<-c("Doublets")
pdf(file = paste0(pro.name,'/pdfs/f9.doubletfinder.pdf'), height = 10, width = 10 )
DimPlot(seu_kidney, reduction = "umap", group.by="Doublets")
dev.off()

sg_seu_kidney <- subset(seu_kidney, subset = Doublets == 'Singlet')
pdf(file = paste0(pro.name,'/pdfs/f10.doubletfinder_singlet.pdf'), height = 10, width = 10 )
DimPlot(sg_seu_kidney, reduction = "umap",group.by="Doublets", cols="blue")
dev.off()
write.csv(seu_kidney@meta.data, file=paste0(pro.name,'/pdfs/DoubletFinder.csv'))

# removal Doublets
temp<-hpn.filt
temp@meta.data<-cbind(temp@meta.data,seu_kidney@meta.data[,9])
colnames(temp@meta.data)[6]<-c("Doublets")
temp<-subset(temp,subset = Doublets == 'Singlet')
temp@meta.data<-temp@meta.data[,-6]

hpn.filt<-temp
hpn.filt
mm<-seu_kidney@meta.data
colnames(mm)[9]<-"DoubletFinder"

library(tidyverse)
theme_set(theme_bw(16))
p1<-mm %>%
  drop_na()%>%
  ggplot(aes(x=seurat_clusters,
             y=nFeature_RNA,
             fill=DoubletFinder))+
  geom_violin()

p2<-mm %>%
  drop_na()%>%
  ggplot(aes(x=seurat_clusters,
             y=nCount_RNA,
             fill=DoubletFinder))+
  geom_violin()

p3<-mm %>%
  drop_na()%>%
  ggplot(aes(x=seurat_clusters,
             y=percent.mito,
             fill=DoubletFinder))+
  geom_violin()

p4<-mm %>%
  drop_na()%>%
  ggplot(aes(x=seurat_clusters,
             y=percent.ribo,
             fill=DoubletFinder))+
  geom_violin()
pdf(file = paste0(pro.name,'/pdfs/f11.doubletfinder_1.pdf'), height = 10, width = 20 )
p1+p2+p3+p4
dev.off()

saveRDS(hpn.filt, file=paste0(pro.name,'/pdfs/hpn.filt.RDS'))

#STEP 4 SCTransform
#hpn.filt<-DietSeurat(hpn.filt)
hpn.filt <- SCTransform(object = hpn.filt, verbose = FALSE)
hpn.filt <- RunPCA(object = hpn.filt, verbose = FALSE)
hpn.filt <- RunUMAP(object = hpn.filt, dims = 1:20, verbose = FALSE)
hpn.filt <- RunTSNE(object = hpn.filt, dims = 1:20, verbose = FALSE)
hpn.filt <- FindNeighbors(object = hpn.filt, dims = 1:20, verbose = FALSE)
hpn.filt <- FindClusters(object = hpn.filt, resolution = 0.5, verbose = FALSE)
pdf(file = paste0(pro.name,'/pdfs/f12.dimplot.umap.pdf'), height = 5, width = 5 )
DimPlot(object = hpn.filt, reduction = "umap", label = TRUE) + ggtitle("sctransform")
dev.off()

saveRDS(hpn.filt, file=paste0(pro.name,'/pdfs/hpn.filt.SCTr.RDS'))

#STEP 5 cellanotation
#5.1 celltype singleR

#PNET<-readRDS("/Volumes/Disk/CINJ/Trial001006/PNET/PNET_62572/pdfs/final.RDS")
#PNET<-SetIdent(PNET, value=PNET@meta.data$seurat_clusters)
#seurat.obj<-PNET

seurat.obj<-hpn.filt
# multi references
bp <- celldex::BlueprintEncodeData()
mona <- MonacoImmuneData()
dice <- DatabaseImmuneCellExpressionData()
dmap <- NovershternHematopoieticData()
hpca <- HumanPrimaryCellAtlasData()
mpd <- MuraroPancreasData()
gpd <- GrunPancreasData()
mpd <- mpd[,!is.na(mpd$label)]
mpd <-logNormCounts(mpd)
gpd <- gpd[,colSums(counts(gpd)) > 0] # remove libraries with no counts.
gpd<-logNormCounts(gpd)

seurat.obj@meta.data$cell.type <- Idents(seurat.obj)
seurat.obj<-DietSeurat(seurat.obj)
test <- as.SingleCellExperiment(seurat.obj)
#Annotated using seven sets
Anno <- SingleR(test = test, ref = list(BP = bp, MONA = mona, DICE=dice, DMAP=dmap, HPCA=hpca), 
                labels = list(bp$label.main, mona$label.main, dice$label.main, dmap$label.main, hpca$label.main), 
                method = "cluster",cluster = test$cell.type)
Anno$cluster <- rownames(Anno)
fin <- Anno %>% tibble::as_tibble() %>% dplyr::select(cluster,labels)

# assigne cell types
seurat.obj<-hpn.filt
new.cluster.ids <- fin$labels
names(new.cluster.ids) <- levels(seurat.obj)
seurat.obj <- RenameIdents(seurat.obj, new.cluster.ids)

pdf(file = paste0(pro.name,'/pdfs/f13.singleR.dimplot.pdf'), height = 8, width = 8 )
DimPlot(seurat.obj, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

#5.2 mapping using Aziumth: #https://azimuth.hubmapconsortium.org/

predictions <- read.delim(paste0(pro.name,'/azimuth_pred.tsv'), row.names = 1)
seurat.obj <- AddMetaData(object = seurat.obj,metadata = predictions)
cCell <- seurat.obj@meta.data$predicted.annotation.l1
tt<-as.data.frame(table(cCell))
tt$Freq
sum(tt$Freq)

pdf(file = paste0(pro.name,'/pdfs/f14.azimuth.hpn.pdf'), height = 8, width = 8 )
DimPlot(seurat.obj, reduction = "umap", group.by = "predicted.annotation.l1",label = TRUE, pt.size = 0.5)
dev.off()

##mapping immune cells
predictions <- read.delim(paste0(pro.name,'/azimuth_pred_i.tsv'), row.names = 1)
seurat.obj <- AddMetaData(object = seurat.obj,metadata = predictions)
cCell <- seurat.obj@meta.data$predicted.celltype.l2
tt<-as.data.frame(table(cCell))
tt$Freq
sum(tt$Freq)

pdf(file = paste0(pro.name,'/pdfs/f15.azimuth.immune.pdf'), height = 8, width = 8 )
DimPlot(seurat.obj, reduction = "umap", group.by = "predicted.celltype.l2",label = TRUE, pt.size = 0.5)
dev.off()

#5.3 mapping to any other public datasets
#final assignment
seurat.obj<-readRDS(paste0(pro.name,'/pdfs/hpn.filt.SCTr.RDS'))
new.cluster.ids=c("0:endocrine","1:alpha","2:stellate","3:endocrine","4:endocrine","5:NK","6:CD16+Monocytes","7:CD8+T","8:CD14+Monocytes","9:endothelial","10:endothelial","11:B","12:alpha","13:endothelial")
names(new.cluster.ids) <- levels(seurat.obj)
seurat.obj<- RenameIdents(seurat.obj, new.cluster.ids)
pdf(file = paste0(pro.name,'/pdfs/f16.final.dimplot.pdf'), height = 8, width = 8 )
DimPlot(seurat.obj, reduction = "umap", label = TRUE, pt.size = 0.5)+ggtitle(pro.name) #8x12 pdf
dev.off()
seurat.obj@meta.data$celltype.cluster<-seurat.obj@active.ident

seurat.obj<-SetIdent(seurat.obj,value = seurat.obj@meta.data$seurat_clusters)
new.cluster.ids=c("endocrine","alpha","stellate","endocrine","endocrine","NK","CD16+Monocytes","CD8+T","CD14+Monocytes","endothelial","endothelial","B","alpha","endothelial")
names(new.cluster.ids) <- levels(seurat.obj)
seurat.obj<- RenameIdents(seurat.obj, new.cluster.ids)
# plot the cell type
pdf(file = paste0(pro.name,'/pdfs/f17.final.dimplot.celltype_new.pdf'), height = 8, width =8 )
DimPlot(seurat.obj, reduction = "umap", label = TRUE, pt.size = 0.5)+ggtitle(pro.name) #8x12 pdf
dev.off()

seurat.obj@meta.data$celltype<-seurat.obj@active.ident
#count cells in each cell type
tmp<-table(Idents(seurat.obj))
tmp<-as.data.frame(tmp)

write.csv(seurat.obj@meta.data,file=paste0(pro.name,'/pdfs/hpn.celltype.done.csv'))
saveRDS(seurat.obj, file=paste0(pro.name,'/pdfs/hpn.celltype.done.RDS'))

