# this is a KP/PG mouse lung tumor

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
pro.name<-"m8967KPg"
mlung.data <-Read10X_h5(paste0("/Volumes/Disk/CINJ/",pro.name,"/mydata/filtered_feature_bc_matrix.h5"))

#STEP 2: QC filter
# create Seurat object
mlung <-CreateSeuratObject(counts = mlung.data, project = pro.name, min.cells = 3, min.features = 200)
# calculate the proportion of mitochondrial reads.
mt.genes <- rownames(mlung)[grep("^mt-",rownames(mlung))]
C<-GetAssayData(object = mlung, slot = "counts")
percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
mlung <- AddMetaData(mlung, percent.mito, col.name = "percent.mito")
# calculate the ribosomal proportion
rb.genes <- rownames(mlung)[grep("^Rp[sl]",rownames(mlung))]
percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
mlung <- AddMetaData(mlung, percent.ribo, col.name = "percent.ribo")

## plot QC
### plot"nFeature_RNA", "nCount_RNA", "percent.mito","percent.ribo"
pdf(file = paste0(pro.name,'/pdfs/f1.vlplot_features_count_mito_ribo.pdf'), height = 5, width = 10 )
VlnPlot(mlung, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.ribo"), ncol = 4, pt.size = 0)
dev.off()

## plot each two features for them
p1<-FeatureScatter(mlung, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+NoLegend()
p2<-FeatureScatter(mlung, feature1 = "nFeature_RNA", feature2 = "percent.mito")+NoLegend()
p3<-FeatureScatter(mlung, feature1 = "nCount_RNA", feature2 = "percent.mito")+NoLegend()
p4<-FeatureScatter(mlung, feature1="percent.ribo", feature2="nFeature_RNA")+NoLegend()
pdf(file =  paste0(pro.name,'/pdfs/f2.vs_features_count_mito_ribo.pdf'), height = 10, width = 10 )
p1+p2+p3+p4
dev.off()
 ## read count vs barcode_rank
pdf(file = paste0(pro.name,'/pdfs/f3.vs_features_count_mito_ribo.pdf'), height = 5, width = 5 )
p<-CalculateBarcodeInflections(mlung,barcode.column = "nCount_RNA",group.column = "orig.ident",threshold.low = NULL, threshold.high = NULL)
BarcodeInflectionsPlot(p)+labs(x="barcode_rank")+ NoLegend()
dev.off()
## histogram
pdf(file = paste0(pro.name,'/pdfs/f4.histo_features_count_mito_ribo.pdf'), height = 7, width = 8 )
par(mfrow=c(2,2))
hist(mlung@meta.data$nCount_RNA, breaks =100, xlab = 'nCount_RNA', main = "")
abline(v = 50000, col = "red")
hist(mlung@meta.data$nFeature_RNA, breaks =100, xlab = 'nFeature_RNA',main = "")
abline(v = 8000, col = "red")
abline(v = 200, col = "red")
hist(mlung@meta.data$percent.mito, breaks =100, xlab = 'percent.mito',main = "")
abline(v = 25, col = "red")
hist(mlung@meta.data$percent.ribo, breaks =100, xlab = 'percent.ribo',main = "")
#abline(v = 5, col = "red")
abline(v = 40, col = "red")
abline(v = 5, col = "red")
dev.off()

#filter by mitochondrial, nfeatures and ncount
mlung.filt <- subset(mlung, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mito < 25 & nCount_RNA<50000 & percent.ribo<40 & percent.ribo>5)
mlung.filt
pdf(file = paste0(pro.name,'/pdfs/f5.vlplot_filtered_features_count_mito_ribo.pdf'), height = 5, width = 10 )
VlnPlot(mlung.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.ribo"), ncol = 4, pt.size = 0)
dev.off()

## convert Seurat to SingleCellExperiment

sce <- as.SingleCellExperiment(mlung.filt)
# Let’s look at what the top 50 expressed genes are. This can be valuable for detecting genes that are overabundant that may be driving a lot of the variation.
pdf(file = paste0(pro.name,'/pdfs/f6.high_exprs.pdf'), height = 10, width = 10 )
plotHighestExprs(sce, exprs_values = "counts")
dev.off()

## remove MALAT1
temp<-mlung.filt
counts <- GetAssayData(temp, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('Malat1'))),]
#counts <- counts[-(which(rownames(counts) %in% mt.genes)),]
temp1 <- subset(temp, features = rownames(counts))
mlung.filt<-temp1
temp1<-NULL
sce <- as.SingleCellExperiment(mlung.filt)
pdf(file = paste0(pro.name,'/pdfs/f7.high_exprs.pdf'), height = 10, width = 10 )
plotHighestExprs(sce, exprs_values = "counts")
dev.off()
sce<-NULL

#STEP 3: doubletfinder 0.023 (3000-5000 cells)


seu_kidney<-mlung.filt
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
nExp_poi <- round(0.023*nrow(seu_kidney@meta.data))  ## Assuming 1.6% doublet formation rate - tailor for your dataset
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

# removal Doublets in mlung.filt
temp<-mlung.filt
temp@meta.data<-cbind(temp@meta.data,seu_kidney@meta.data[,9])
colnames(temp@meta.data)[6]<-c("Doublets")
temp<-subset(temp,subset = Doublets == 'Singlet')
temp@meta.data<-temp@meta.data[,-6]

mlung.filt<-temp
mlung.filt
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

saveRDS(mlung.filt, file=paste0(pro.name,'/pdfs/mlung.filt.RDS'))

mlung.filt<-readRDS(file=paste0(pro.name,'/pdfs/mlung.filt.RDS'))

#STEP 4 SCTransform
#mlung.filt<-DietSeurat(mlung.filt), for singleR error
mlung.filt <- SCTransform(object = mlung.filt, vst.flavor = "v2",verbose = FALSE) # using SCTransform v2 Jan-11, 2022
mlung.filt <- RunPCA(object = mlung.filt, npcs = 30, verbose = FALSE)
mlung.filt <- RunUMAP(object = mlung.filt, dims = 1:30, verbose = FALSE)
mlung.filt <- RunTSNE(object = mlung.filt, dims = 1:30, verbose = FALSE)
mlung.filt <- FindNeighbors(object = mlung.filt, dims = 1:30, verbose = FALSE)
mlung.filt <- FindClusters(object = mlung.filt, resolution = 1, verbose = FALSE)
pdf(file = paste0(pro.name,'/pdfs/f12.dimplot.umap.pdf'), height = 5, width = 5 )
DimPlot(object = mlung.filt, reduction = "umap", label = TRUE) + ggtitle("sctransform")
dev.off()

pdf(file = paste0(pro.name,'/pdfs/f12.dimplot.tsne.pdf'), height = 5, width = 5 )
DimPlot(object = mlung.filt, reduction = "tsne", label = TRUE) + ggtitle("sctransform")
dev.off()

saveRDS(mlung.filt, file=paste0(pro.name,'/pdfs/mlung.filt.SCTr.RDS'))

#STEP 5 cellanotation
#5.1 celltype singleR


seurat.obj<-mlung.filt
# multi references
immgen.se <- ImmGenData()  # 830 samples
mouserna.se <- MouseRNAseqData() # 358 samples


seurat.obj@meta.data$cell.type <- Idents(seurat.obj)
seurat.obj<-DietSeurat(seurat.obj)
test <- as.SingleCellExperiment(seurat.obj)
#Annotated using seven sets
Anno <- SingleR(test = test, ref = list(IM = immgen.se, MR = mouserna.se), 
                labels = list(immgen.se$label.fine, mouserna.se$label.fine), 
                method = "cluster",cluster = test$cell.type)
Anno$cluster <- rownames(Anno)
fin <- Anno %>% tibble::as_tibble() %>% dplyr::select(cluster,labels)

# assigne cell types
seurat.obj<-mlung.filt
new.cluster.ids <- fin$labels
names(new.cluster.ids) <- levels(seurat.obj)
seurat.obj <- RenameIdents(seurat.obj, new.cluster.ids)

pdf(file = paste0(pro.name,'/pdfs/f13.singleR.dimplot.pdf'), height = 12, width = 12 )
DimPlot(seurat.obj, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

pdf(file = paste0(pro.name,'/pdfs/f13.singleR.dimplot.tsne..pdf'), height = 12, width = 12 )
DimPlot(seurat.obj, reduction = "tsne", label = TRUE, pt.size = 3)
dev.off()

seurat.obj@meta.data$singler<-seurat.obj@active.ident


#5.2 mapping using Aziumth: #https://azimuth.hubmapconsortium.org/ (convert genes to human and use human lung and human PMBC)

predictions <- read.delim(paste0(pro.name,'/azimuth_pred.tsv'), row.names = 1)
seurat.obj <- AddMetaData(object = seurat.obj,metadata = predictions)
cCell <- seurat.obj@meta.data$predicted.ann_level_3
tt<-as.data.frame(table(cCell))
tt$Freq
sum(tt$Freq)

pdf(file = paste0(pro.name,'/pdfs/f14.azimuth.mlung.pdf'), height = 12, width = 12 )
DimPlot(seurat.obj, reduction = "umap", group.by = "predicted.ann_level_3",label = TRUE, pt.size = 0.5)
dev.off()

pdf(file = paste0(pro.name,'/pdfs/f14.azimuth.mlung.tsne.pdf'), height = 12, width = 12 )
DimPlot(seurat.obj, reduction = "tsne", group.by = "predicted.ann_level_3",label = TRUE, pt.size = 0.5)
dev.off()

seurat.obj<-SetIdent(seurat.obj, value = seurat.obj@meta.data$seurat_clusters)

#5.3 mapping to public data-sets

#5.3.1 #https://journals.biologists.com/dev/article-abstract/148/24/dev199512/273783/A-single-cell-atlas-of-mouse-lung-development?redirectedFrom=fulltext
#markers for endothelial subsets
marker_genes <- c("Gpihbp1", "Kit", # Microvascular
                  "Car4", "Kdr", # Car4
                  "Mki67", "Top2a", # Proliferating miEC
                  "Vwf", "Vcam1", #macrovascular
                  "Vegfc", "Prss23", #Venous macrovascular
                  "Cxcl12", "Pcsk5", #Artearial macrovascular
                  "Ephb4", # Vein
                  "Flt4", "Ccl21a" # Lymphatic
)
FeaturePlot(mlung.filt, reduction="tsne",features = marker_genes)

#markers for epithelial subsets
marker_genes <- c("Sftpa1", "Sftpc", # AT2
                  "Hopx", "Aqp5", #AT1
                  "Foxj1", "Dynlrb2", # Ciliated
                  "Mdk", "Mki67", # Primordial
                  "Scgb1a1", "Scgb3a2", # Secretory
                  "Cdkn1a", "Krt8", # Transitional?
                  "Ascl1", "Scg5" #Neuroendocrine
)
FeaturePlot(mlung.filt, reduction="tsne", features = marker_genes)


#5.3.2 mapping immune cells
#predictions <- read.delim(paste0(pro.name,'/azimuth_pred_i.tsv'), row.names = 1)
#seurat.obj <- AddMetaData(object = seurat.obj,metadata = predictions)
#cCell <- seurat.obj@meta.data$predicted.celltype.l2
#tt<-as.data.frame(table(cCell))
#tt$Freq
#sum(tt$Freq)

#pdf(file = paste0(pro.name,'/pdfs/f15.azimuth.immune.pdf'), height = 8, width = 8 )
#DimPlot(seurat.obj, reduction = "umap", group.by = "predicted.celltype.l2",label = TRUE, pt.size = 0.5)
#dev.off()

#check marker genes, signature score

# mapping to mouse atlas https://github.com/ggjlab/scMCA
#520,000 cells
#library(devtools)
#install_github("ggjlab/scMCA")
library(scMCA)

mtx<-as.matrix(mlung@assays$RNA@counts)
write.csv(mtx, file="mtx.csv", quote = FALSE)

mca_result <- scMCA(scdata =mtx, numbers_plot = 3)

tmp<-mca_result$scMCA_probility
tmp1<-distinct(tmp, tmp$Cell, .keep_all = TRUE)
rownames(tmp1)<-tmp1$Cell
tmp1<-tmp1[,-4]
colnames(tmp1)<-c("mca_cell","mca_celltype","mca_score")

seurat.obj<-AddMetaData(seurat.obj, metadata = tmp1)

DimPlot(seurat.obj, reduction = "umap", group.by = "mca_celltype", label = TRUE, pt.size = 5)+
  theme(axis.line.y=element_line(linetype=1,color="black",size=1),axis.line.x=element_line(linetype=1,color="black",size=1))

write.csv(seurat.obj@meta.data, file="all.csv")

tmp<-seurat.obj@meta.data %>% filter(seurat_clusters=="19")
tmp1<-as.matrix(table(tmp$mca_celltype))

# MCA,keep majority based on clusters
new.cluster.ids=c("0:Alveolar macrophage_Pclaf high(Lung)", "1:Alveolar macrophage_Pclaf high(Lung)","2:AT2 Cell(Lung)","3:B","4:AT2 Cell(Lung)",
                  "5:Neutrophil granulocyte(Lung)","6:AT2 Cell(Lung)","7:Conventional dendritic cell_Mgl2 high(Lung)","8:Nuocyte(Lung)/T cell","9:T cell",
                  "10:Interstitial macrophage(Lung)","11:T cell (deviding)","12:Endothelial cell_Tmem100 high(Lung)","13:Endothelial cell_Tmem100 high(Lung)/AT2",
                  "14:Alveolar macrophage_Pclaf high(Lung)/AT2","15:Alveolar macrophage_Pclaf high(Lung)","16:Interstitial macrophage(Lung)/DC",
                  "17:Conventional dendritic cell_Gngt2 high(Lung)","18:AT2 Cell(Lung)/Stromal","19:DC(deviding)")

seurat.obj<-SetIdent(seurat.obj, value = seurat.obj@meta.data$seurat_clusters)
names(new.cluster.ids) <- levels(seurat.obj)
seurat.obj<- RenameIdents(seurat.obj, new.cluster.ids)
seurat.obj@meta.data$mca_majority<-seurat.obj@active.ident

DimPlot(seurat.obj, reduction = "umap", group.by = "mca_majority", label = TRUE, pt.size = 5)+
  theme(axis.line.y=element_line(linetype=1,color="black",size=1),axis.line.x=element_line(linetype=1,color="black",size=1))
write.csv(seurat.obj@meta.data, file="all.csv")

# mapping to mouse atlas MOCA, 1,331,984 high quality cells
#https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/landing

## awk -F'.' '{print $1'\t'$3'\tGene Expression'}'


MOCA.data<-Read10X("MOCA/filtered_feature_bc_matrix/")
MOCA<-CreateSeuratObject(counts=MOCA.data, project="MOCA", min.cells = 3, min.features = 200)

cell_meta<-read.csv("MOCA/cell_annotate.csv")
rownames(cell_meta)<-cell_meta$sample
MOCA<-AddMetaData(MOCA, meta = cell_meta)
saveRDS(MOCA, file="MOCA.seurat.RDS")
saveRDS(seurat.obj, file = "mKP8463g/pdfs/")
## mapping
## go amarel
MOCA<-readRDS("MOCA.seurat.RDS")
MOCA<-SetIdent(MOCA, value=MOCA@meta.data$Main_cell_type)

# get 10000 for each cell types
#MOCA_sub = MOCA[, sample(1:ncol(MOCA),round(ncol(MOCA)/10)) ]
sample.size = 100000
orig.size = ncol(MOCA)
set.seed(111)
sampled.cells<-colnames(MOCA)[sample(x = orig.size, size = sample.size, replace = F)]
MOCA_sub <- subset(MOCA, cells = sampled.cells)
saveRDS(MOCA_sub, file="MOCA.sub.RDS")

MOCA_sub<-SCTransform(MOCA_sub)
#MOCA.epithelial<- FindVariableFeatures(MOCA.epithelial, selection.method = "vst", nfeatures = 2000)
MOCA_sub<-ScaleData(MOCA_sub, verbose = FALSE)
MOCA_sub<-RunPCA(MOCA_sub,verbose = FALSE)
#MOCA.epithelial<- RunPCA(MOCA.epithelial, features = VariableFeatures(object = MOCA.epithelial))
MOCA_sub <- RunUMAP(MOCA_sub, reduction = "pca", dims = 1:40, verbose = FALSE)

t.anchors <- FindTransferAnchors(reference = MOCA_sub, query = tmp,dims = 1:40, reference.reduction = "pca")
predictions <- TransferData(anchorset = t.anchors, refdata = MOCA_sub$Main_cell_type,
                            dims = 1:40)
tmp <- AddMetaData(tmp, metadata = predictions)




#final assignment
#seurat.obj<-readRDS(paste0(pro.name,'/pdfs/mlung.filt.SCTr.RDS'))
DimPlot(seurat.obj, reduction = "tsne", label = T, group.by="seurat_clusters", pt.size = 3, label.size = 4)+ggtitle("m8967KPg")
DimPlot(seurat.obj, reduction = "tsne", label = T, group.by="singler", pt.size = 3, label.size = 4)+ggtitle("m8967KPg")
DimPlot(seurat.obj, reduction = "tsne", label = T, group.by="predicted.ann_level_3", pt.size = 3, label.size = 4)+ggtitle("m8967KPg")


new.cluster.ids=c("0:Macrophages (MF)","1:Epithelial (AT1)","2:Epithelial (AT1)","3:Neutrophils","4:Macrophage (MF.103-11B+24-)",
                  "5:Macrophage (MF.103−11B+24-)","6:T","7:Endothelial",
                  "8:Epithelial (AT2)","9:Endothelial","10:DC (CD103+11B-24+)","11:Epithelial(AT1)","12:B",
                  "13:Endothelial","14:Epithelial (AT1)","15:Rare","16:Macrophages (MFIO5.II+480INT)","17:Macrophages (MFAR-)",
                  "18:Macrophages (MFAR-)", "19:Epithelial (AT2)","20:Macrophages (MFAR-)","21:Epithelial (AT1)","22:Neutrophils",
                  "23:Macrophages (MF)","24:T","25:Fibroblasts")
names(new.cluster.ids) <- levels(seurat.obj)
seurat.obj<- RenameIdents(seurat.obj, new.cluster.ids)
seurat.obj@meta.data$celltype.cluster<-seurat.obj@active.ident

pdf(file = paste0(pro.name,'/pdfs/f16.final.dimplot.pdf'), height = 12, width = 12 )
DimPlot(seurat.obj, reduction = "umap", label = TRUE, pt.size = 3, group.by="celltype.cluster")+ggtitle(pro.name) #8x12 pdf
dev.off()

seurat.obj<-SetIdent(seurat.obj,value = seurat.obj@meta.data$seurat_clusters)


new.cluster.ids=c("Macrophages (MF)","Epithelial (AT1)","Epithelial (AT1)","Neutrophils","Macrophage (MF.103-11B+24-)",
                  "Macrophage (MF.103-11B+24-)","T","Endothelial",
                  "Epithelial (AT2)","Endothelial","DC (CD103+11B-24+)","Epithelial (AT1)","B",
                  "Endothelial","Epithelial (AT1)","Rare","Macrophages (MFIO5.II+480INT)","Macrophages (MFAR-)",
                  "Macrophages (MFAR-)", "Epithelial (AT2)","Macrophages (MFAR-)","Epithelial (AT1)","Neutrophils",
                  "Macrophages (MF)","T","Fibroblasts")

names(new.cluster.ids) <- levels(seurat.obj)
seurat.obj<- RenameIdents(seurat.obj, new.cluster.ids)
seurat.obj@meta.data$celltype<-seurat.obj@active.ident
#count cells in each cell type
tmp<-table(Idents(seurat.obj))
tmp<-as.data.frame(tmp)


pdf(file = paste0(pro.name,'/pdfs/f17.final.dimplot.celltype.pdf'), height = 12, width =12 )
DimPlot(seurat.obj, reduction = "umap", label = TRUE, pt.size = 3)+ggtitle(pro.name) #8x12 pdf
dev.off()


pdf(file = paste0(pro.name,'/pdfs/f17.final.dimplot.celltype.tsne.pdf'), height = 10, width =12 )
DimPlot(seurat.obj, reduction = "tsne", label = TRUE, pt.size = 5)+ggtitle(pro.name) #8x12 pdf
dev.off()


write.csv(seurat.obj@meta.data,file=paste0(pro.name,'/pdfs/mlung.celltype.done.csv'))
saveRDS(seurat.obj, file=paste0(pro.name,'/pdfs/mlung.celltype.done.RDS'))

DimPlot(seurat.obj, reduction = "umap", label = T, group.by="celltype", pt.size = 3, label.size = 4)+ggtitle("m8967KPg")

seurat.obj<-readRDS(file=paste0(pro.name,'/pdfs/mlung.celltype.done.RDS'))
###Step 6:### cell cycle
# for Mus musculus
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)

mlung<-seurat.obj
mlung <- CellCycleScoring(mlung, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

####cutoff for cycling and non-cycling
cutter_s<-0.25
cutter_g2m<-0.25
## mean+2SD
temp<-mlung
cutter_s<-mean(temp@meta.data$S.Score)+1.5*sd(temp@meta.data$S.Score)
cutter_g2m<-mean(temp@meta.data$G2M.Score)+1*sd(temp@meta.data$G2M.Score)
###
pdf(file = paste0(pro.name,'/pdfs/f18.cellcycle_1d.pdf'), height = 8, width =8 )
plot(temp@meta.data$S.Score,temp@meta.data$G2M.Score, xlab="S Score", ylab="G2M Score")
abline(h =cutter_g2m, v = cutter_s, col = "red")
dev.off()
#plot(temp@meta.data$S.Score,temp@meta.data$G2M.Score, xlab="S Score", ylab="G2M Score", group.by="cycling")

mm<-temp@meta.data
mm$cycling_1d<-"PD"
## two color

for(i in 1:dim(mm)[1]){
  if((mm$S.Score[i] <= cutter_s) & (mm$G2M.Score[i]<=cutter_g2m)){
    mm$cycling_1d[i]<-"Non-Cycling"
  }
  if((mm$S.Score[i] > cutter_s) & (mm$G2M.Score[i]>cutter_g2m)){
    mm$cycling_1d[i]<-"Cycling"
  }
  if((mm$S.Score[i] >cutter_s) & (mm$G2M.Score[i]<cutter_g2m)){
    mm$cycling_1d[i]<-"Cycling"
  }
  if((mm$S.Score[i] < cutter_s) & (mm$G2M.Score[i]>cutter_g2m)){
    mm$cycling_1d[i]<-"Cycling"
  }
}

temp@meta.data<-mm
pdf(file = paste0(pro.name,'/pdfs/f19.cycling_vs_nocycling_1d.tsne.pdf'), height = 6, width =6 )
DimPlot(temp,reduction="tsne", group.by ="cycling_1d",pt.size = 2)
dev.off()

seurat.obj@meta.data<-mm

### Step 7:## copykat
mlung<-seurat.obj
exp.rawdata <- as.matrix(mlung@assays$RNA@counts)
rownames(exp.rawdata)<-toupper(rownames(exp.rawdata))
copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", cell.line="no", ngene.chr=1, win.size=25, KS.cut=0.15, sam.name="test", distance="euclidean", n.cores=10)


## write to Seurat

tt<-mlung@meta.data
cnv<-data.frame(copykat.test$prediction)
zz<-merge(tt, cnv, by="row.names", all =TRUE)
zz[is.na(zz)]<-"pd"
a<-zz
#a<-a[,-18]
row.names(a)<-zz[,1]
a<-a[,-1]
#a<-a[,-15]
mlung@meta.data<-a

pdf(file = paste0(pro.name,'/pdfs/f20.copykat.pdf'), height = 8, width =8 )
DimPlot(mlung, reduction = "tsne", group.by="copykat.pred", pt.size = 3)
dev.off()
temp<-table(mlung@meta.data$copykat.pred)
temp<-as.data.frame(temp)
write.csv(a, file=paste0(pro.name,'/pdfs/copykat.csv'))
saveRDS(mlung, file=paste0(pro.name,'/pdfs/final.RDS'))

## cell type top n heatmap opti
#devtools::install_github('caleblareau/BuenColors')
#utils::install.packages(pkgs = "ggstatsplot")
library(BuenColors)
library(ggstatsplot)
pro.name<-"mKP8634g"
mlung<-readRDS(file=paste0(pro.name,'/pdfs/final.RDS'))
mlung.markers <- FindAllMarkers(mlung, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#mlung.markers %>%
#  group_by(cluster) %>%
#  top_n(n = 5, wt = avg_log2FC) -> top5


library(patchwork)
cg <- c("CD4","CD3E","IL7R", "KLF2", "CCR7","TCF7",
        "SELL", "CCL4", "CCL5", "PRF1",  "GZMB",
        "GZMK", "FGFBP2", "CX3CR1", "RORC",
        "CXCR5", "FOXP3", "GATA3", "PTGDR2")

marker_selected_1 <- mlung.markers  %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(pct.1 >= 0.5 & pct.2 <= 0.6) %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::slice_max(order_by = avg_log2FC, n = 10)
DoHeatmapPlot(object = mlung, groupBy = 'celltype', features = marker_selected_1$gene)
DoHeatmapPlot(object = mlung, groupBy = 'celltype', features = cg)

#count percentage of cell types

#utils::install.packages(pkgs = "ggstatsplot")
library(ggstatsplot)
source('PropPlot.R')
table(mlung$seurat_clusters, mlung$celltype)
PropPlot(object = mlung, groupBy = 'celltype') + 
  PropPlot(object = mlung, groupBy = 'seurat_clusters') 
PropPlot(object = mlung, groupBy = 'celltype') + 
  PropPlot(object = mlung, groupBy = 'seurat_clusters') +PropPlot(object = mlung, groupBy = 'cycling_2d')+PropPlot(object = mlung, groupBy = 'copykat.pred')

## CD8T

DimPlot(mlung, reduction = "tsne", label = T, pt.size = 4, label.size = 3)

t<-subset(mlung, ident="CD8+T")

CD8T.matrix<-as.matrix(t@assays$RNA@counts)

CD8T<-CreateSeuratObject(counts = CD8T.matrix, project = "CD8T")
CD8T<-SCTransform(object = CD8T)
CD8T<-RunPCA(object = CD8T)
CD8T<-RunUMAP(object = CD8T, dims = 1:10)
CD8T<-RunTSNE(object = CD8T,dims = 1:10)
CD8T<-FindNeighbors(object = CD8T, dims = 1:10)
CD8T<-FindClusters(object = CD8T, resolution = 1)

DimPlot(CD8T, reduction = "tsne", label = T, pt.size = 8, label.size = 3)

## check markers

#Naive CD8+ T cell
VlnPlot(object = CD8T, features =c("Ccr7","Cd28","Hmgb2","Tcf7"))
#  Cytotoxic CD8+ T cell
# CD8+, CCR7-	IFN-γ+, Perforin+, Granzyme+	
VlnPlot(CD8T, features = c("Gzmb","Hmgb2","Klrg1","Ifng", "Ccr7"))
# Exhausted 
VlnPlot(CD8T, features = c("Lag3","Hmgb2","Pdcd1","Havcr2","Tigit"))
# Figure 7: tumor vs normal of alpha sub-populations

# Memeory T
VlnPlot(CD8T, features = c("Sell","Ccr7","Tigit","Pdcd1","Lag3","Tcf7"))

#Activated T
VlnPlot(CD8T, features = c("Tnf","Ifng","Fos","Jun"))

# effector
VlnPlot(CD8T, features = c("Prf1","Gzmb","Gzmm","Klrg1","Fgfbp2","Klrd1","Cd44"))

CD8T<-SetIdent(CD8T, value = CD8T@meta.data$seurat_clusters)
new.cluster.ids <- c("naive","naive","exhausted","effector")
names(new.cluster.ids) <- levels(CD8T)
CD8T <- RenameIdents(CD8T, new.cluster.ids)
CD8T@meta.data$celltype<-CD8T@active.ident

DimPlot(CD8T, reduction = "tsne", label = T, pt.size = 8, label.size = 3, group.by = "seurat_clusters")
saveRDS(CD8T, file = "CD8T.RDS")

#  immune checkpoint

#stimulatory
DotPlot(mlung, features = c("Cd27","Cd28", "Cd4", "Cd122", "Cd137", "Ox40", "Gitr", "Icos"))
#Inhibitory checkpoint molecules
DotPlot(mlung, features = c("A2ar","B7-H3","B7-H4","Btla","Ctla-4","Ido","Kir","Lag3","Nox2","PD-1","Tim-3","Vista","Siglec7"))
                            
                            
#https://www.bio-rad-antibodies.com/immune-checkpoint-antibody.html?JSESSIONID_STERLING=A9B71F9D8D0136914B577C006B513A34.ecommerce1&evCntryLang=US-en&cntry=US&thirdPartyCookieEnabled=true

DotPlot(mlung, features = c("Cd137","Cd137l","Cd80","Cd86","Cd27","Cd70","Cd28","Cd40","Cd154","Cd122","Cd270","Havem","Cd278","Icos","Gitr","CD357","Gitrl","Cd134","Ox40","Cd252"))
