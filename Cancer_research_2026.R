###################################
### Fig 4.                 ########
###################################

###############################################################################
library(SAVER)
SCRIPTS<-"/Users/xyy/Dropbox (Personal)/single_cell/code/"
#source(paste(SCRIPTS,"seuratMergeWrapper.R",sep = ""))
source(paste(SCRIPTS, "seuratWrapper.R", sep = ""))
library(phateR) 
library(dplyr)   # define funciton %>%
setwd("/Users/xyy/Dropbox (Personal)/single_cell/result")

load('wk4_saver.rda')
load('wk2_saver.rda')
load('ctr_saver.rda')

SCRIPTS<-"/Users/xyy/Dropbox (Personal)/single_cell/code/"
source(paste(SCRIPTS, "seuratWrapper.R", sep = ""))
library(phateR) 
library(dplyr)   # define funciton %>%
setwd("/Users/xyy/Dropbox (Personal)/single_cell/result")

ctr = CreateSeuratObject(ctr.saver)
ctr = NormalizeData(ctr)
ctr = ScaleData(ctr)

wk2 = CreateSeuratObject(wk2.saver)
wk2 = NormalizeData(wk2)
wk2 = ScaleData(wk2)

wk4 = CreateSeuratObject(wk4.saver)
wk4 = NormalizeData(wk4)
wk4 = ScaleData(wk4)

rm(ctr.saver, wk2.saver, wk4.saver)

ctr[["percent.mt"]] <- PercentageFeatureSet(ctr, pattern = "^mt-")
wk2[["percent.mt"]] <- PercentageFeatureSet(wk2, pattern = "^mt-")
wk4[["percent.mt"]] <- PercentageFeatureSet(wk4, pattern = "^mt-")


ctr$stim = 'ctr'
wk2$stim = 'wk2'
wk4$stim = 'wk4'

####plot quality control
VlnPlot(ctr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(wk2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(wk4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

##########################################################################
####  Integration 
##########################################################################
anchor_dim = 30

set.seed(6)
sox2.anchors <- FindIntegrationAnchors(object.list = list(ctr, wk2, wk4), dims = 1:anchor_dim) # may change 30 to 50
sox2.combined <- IntegrateData(anchorset = sox2.anchors, dims = 1:anchor_dim)
DefaultAssay(sox2.combined) <- "integrated"

#select number of PCs
set.seed(6)
sox2.combined <- ScaleData(sox2.combined, verbose = FALSE)
sox2.combined <- RunPCA(sox2.combined, npcs = 100, verbose = FALSE)

DimHeatmap(sox2.combined, dims = 1:15, cells = 8500, balanced = TRUE)
DimHeatmap(sox2.combined, dims = 16:35, cells = 8500, balanced = TRUE)

sox2.combined <- JackStraw(sox2.combined, num.replicate = 100, dims = 30)
sox2.combined <- ScoreJackStraw(sox2.combined, dims = 1:29)
JackStrawPlot(sox2.combined, dims = 1:29)

# Run the standard workflow for visualization and clustering
pca_dim = 30
sox2.combined <- ScaleData(sox2.combined, verbose = FALSE)
sox2.combined <- RunPCA(sox2.combined, npcs = pca_dim, verbose = FALSE)
#UMap and Clustering
sox2.combined <- FindNeighbors(sox2.combined, reduction = "pca", dims = 1:pca_dim) # dim = 50 or 30
sox2.combined <- FindClusters(sox2.combined, resolution = 0.5)
# Visualizationap", split.by = "stim")
sox2.combined <- RunUMAP(sox2.combined, reduction = "pca", dims = 1:30,  n.neighbors = 50)

p1 <- DimPlot(sox2.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(sox2.combined, reduction = "umap", label = TRUE)

CombinePlots(plots = list(p1, p2))

genes = c( 'Trdc', 'Cd2', 'Trac', 'Cd4', 'Cd8a', 'Cd8b1', 'Foxp3', 'Trbc1',
		'Trbc2', 'Gzmb', 'Eomes', 'Icos', 'Cd3g', 'Cd3e', 'Ifng', 'Ncr1', 'Cd19',
		'Cd79a', 'Cd79b') 

pdf('Integrated_Saver_all_cell_marker.pdf', width = 20, height = 16)
FeaturePlot(sox2.combined, features =genes)
dev.off()


pdf('Integrated_Saver_Umap_all_cell.pdf', width = 20, height = 10)
CombinePlots(plots = list(p1, p2))
dev.off()

pdf('Integrated_Umap_Saver_sep_all_cell.pdf', width = 30, height = 10)
DimPlot(sox2.combined, reduction = "umap", split.by = "stim")
dev.off()

marker_ctr = FindAllMarkers(sox2.combined,  only.pos = TRUE)
top_ctr = marker_ctr %>% group_by(cluster) %>% top_n(n=30, wt=avg_logFC)

write.csv(cbind(top_ctr$gene, top_ctr$cluster), file = 'genes__list_integrated_all_cell.csv')

pdf('heatmap_integrated_top100_all_cell.pdf', height = 20, width = 20)
DoHeatmap(sox2.combined, features = top_ctr$gene)
dev.off()

#####################SingleR Annotation
library(SingleR)
dat = as.SingleCellExperiment(sox2.combined)
ImmGe = ImmGenData()

counts <- GetAssayData(sox2.combined)
sox2_clusters = sox2.combined$seurat_clusters

sox2_singler_cluster <- SingleR(test = counts, ref = ImmGe, labels=ImmGe$label.main,  method = "cluster", clusters = sox2_clusters, de.method="wilcox")
#sox2_sub_singler_cluster_cell <- SingleR(test = counts, ref = ImmGe, labels=ImmGe$label.main,  clusters = sox2_sub_clusters, de.method="wilcox")


A = GetAssayData(sox2.combined)
id = match(genes, rownames(A))
write.csv(A[id, ], file = 'all_cell_marker_expression.csv')
write.csv(sox2.combined$seurat_clusters, 'cluster.csv')


#########################
#Fig 4I
##########
library("circlize")
load('sox2_combined_07_14_2020.Rda')
mye_cluster = c(1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 16, 17, 18, 19, 21, 22, 23, 26)


dat = cbind(sox2.combined$stim, sox2.combined$seurat_clusters,
		sox2.combined$seurat_clusters, sox2.combined$seurat_clusters) 
dat[, 4] = as.character(as.numeric(dat[, 4]))

id = dat[, 3] %in% mye_cluster
dat[id, 3] = 'Mye'
dat[!id, 3] = 'Lym'

sum((dat[, 1] %in% 'ctr') + (dat[, 3] %in% 'Mye') == 2) 
sum((dat[, 1] %in% 'ctr') + (dat[, 3] %in% 'Lym') == 2) 
sum((dat[, 1] %in% 'wk2') + (dat[, 3] %in% 'Mye') == 2) 
sum((dat[, 1] %in% 'wk2') + (dat[, 3] %in% 'Lym') == 2) 
sum((dat[, 1] %in% 'wk4') + (dat[, 3] %in% 'Mye') == 2) 
sum((dat[, 1] %in% 'wk4') + (dat[, 3] %in% 'Lym') == 2) 

flow_dat = data.frame(orig_reg = rep(c('Ctr', 'Wk2', 'Wk4'), each = 2), dest_reg = rep(c('Myeloid', 'Lymphoid'), 3),
		flow = c(370, 694, 1997, 230, 5151, 161))


cluster_id = as.character(1:27)

flow_dat2 = data.frame(orig_reg = rep(c('Ctr', 'Wk2', 'Wk4'), each = 27), dest_reg = rep(cluster_id, 3),
				flow =0)

for(i in 1:27){
	flow_dat2[i , 3]     = sum((dat[, 1] %in% 'ctr') + (dat[, 4] %in% cluster_id[i]) == 2)
	flow_dat2[27 + i, 3] = sum((dat[, 1] %in% 'wk2') + (dat[, 4] %in% cluster_id[i]) == 2)
	flow_dat2[54 + i, 3] = sum((dat[, 1] %in% 'wk4') + (dat[, 4] %in% cluster_id[i]) == 2)
}


manualcolors<-c('black', 'orange', 'cornflowerblue', 
		'magenta', 'darkolivegreen4', 'indianred1', 'tan4', 
		'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3',
		"tan1",'darkgray', 'wheat4', '#DDAD4B', 'chartreuse', 
		'seagreen1', 'moccasin', 'mediumvioletred', 'seagreen','cadetblue1',
		"darkolivegreen1" ,"tan2" ,   "tomato3" , "#7CE3D8","gainsboro")

flow_dat1 = data.frame(region = c('Ctr', 'Wk2', 'Wk4', "Myeloid", 'Lymphoid'), order = 1:5, col1 = c('red' ,'green','blue', '#FECC2F', '#A463D7')  )


flow_dat3 = data.frame(region = c('Ctr', 'Wk2', 'Wk4',cluster_id), order = 1:30, col1 = c('red' ,'green','blue',  manualcolors)  )

circos.clear()

circos.par(start.degree = 0, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

 
 chordDiagram(x = flow_dat, grid.col = flow_dat1$col, transparency = 0.25,
              order = flow_dat1$region, directional = 1,
              direction.type = c("arrows", "diffHeight"), diffHeight = -0.04,
              annotationTrack = "grid", annotationTrackHeight = c(0.05, 0.1),
              link.arr.type = "big.arrow", link.sort = TRUE, link.largest.ontop = TRUE)

circos.trackPlotRegion(
   track.index = 1, 
   bg.border = NA, 
   panel.fun = function(x, y) {
     xlim = get.cell.meta.data("xlim")
     sector.index = get.cell.meta.data("sector.index")
     reg1 = flow_dat1$region[flow_dat1$region == sector.index]
     circos.text(x = mean(xlim), y = 3.8, labels = reg1, facing = "bending", cex = 1.4)
     circos.axis(h = "top", 
                 major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2] > 10, yes = 500, no = 100)), 
                 minor.ticks = 1, major.tick.percentage = 0.2,
                 labels.niceFacing = FALSE)
   }
 )
 
 dev.copy2pdf(file = "cfplot_sox2_Mye_Lym_07_16_2020.pdf", height=10, width=10)

circos.clear()

circos.par(start.degree = -20, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

chordDiagram(x = flow_dat2, grid.col = flow_dat3$col, transparency = 0.25,
		order = flow_dat3$region, directional = 1,
		direction.type = c("arrows", "diffHeight"), diffHeight = -0.04,
		annotationTrack = "grid", annotationTrackHeight = c(0.05, 0.1),
		link.arr.type = "big.arrow", link.sort = TRUE, link.largest.ontop = TRUE)

circos.trackPlotRegion(
		track.index = 1, 
		bg.border = NA, 
		panel.fun = function(x, y) {
			xlim = get.cell.meta.data("xlim")
			sector.index = get.cell.meta.data("sector.index")
			reg1 = flow_dat3$region[flow_dat3$region == sector.index]
			circos.text(x = mean(xlim), y = 3.2, labels = reg1, facing = "bending", cex = 1.4)
		}
)

dev.copy2pdf(file = "cfplot_sox2_clusters_07_16_2020.pdf", height=10, width=10)



###################################
### Fig 5.                 ########
###################################
library(SAVER)
library(readr)
library(viridis)
source(paste(SCRIPTS, "seuratWrapper.R", sep = ""))
library(phateR) 
library(dplyr)   

load('wk4_saver_05_20_2020.rda')
load('wk2_saver_05_20_2020.rda')
load('ctr_saver_05_20_2020.rda')

SCRIPTS<-"/Users/xyy/Dropbox (Personal)/single_cell/code/"
source(paste(SCRIPTS, "seuratWrapper.R", sep = ""))
library(phateR) 
library(dplyr)   # define funciton %>%
setwd("/Users/xyy/Dropbox (Personal)/single_cell/result")

ctr = CreateSeuratObject(ctr.saver)
ctr = NormalizeData(ctr)
ctr = ScaleData(ctr)

wk2 = CreateSeuratObject(wk2.saver)
wk2 = NormalizeData(wk2)
wk2 = ScaleData(wk2)

wk4 = CreateSeuratObject(wk4.saver)
wk4 = NormalizeData(wk4)
wk4 = ScaleData(wk4)

rm(ctr.saver, wk2.saver, wk4.saver)

ctr[["percent.mt"]] <- PercentageFeatureSet(ctr, pattern = "^mt-")
wk2[["percent.mt"]] <- PercentageFeatureSet(wk2, pattern = "^mt-")
wk4[["percent.mt"]] <- PercentageFeatureSet(wk4, pattern = "^mt-")


ctr$stim = 'ctr'
wk2$stim = 'wk2'
wk4$stim = 'wk4'

####plot quality control
VlnPlot(ctr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(wk2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(wk4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
##########################################################################
####  Integration 

#anchor_dim = 20
anchor_dim = 30
sox2.anchors <- FindIntegrationAnchors(object.list = list(ctr, wk2, wk4), dims = 1:anchor_dim) # may change 30 to 50
sox2.combined <- IntegrateData(anchorset = sox2.anchors, dims = 1:anchor_dim)
DefaultAssay(sox2.combined) <- "integrated"

sox2.combined
genes = c( 'Trdc', 'Cd2', 'Trac', 'Cd4', 'Cd8a', 'Cd8b1', 'Foxp3', 'Trbc1',
		'Trbc2', 'Gzmb', 'Eomes', 'Icos', 'Cd3g', 'Cd3e', 'Ifng', 'Ncr1', 'Cd19',
		'Cd79a', 'Cd79b') 

load('mye_cell_id.rda') 
mye = SubsetData(sox2.combined, cells = mye_cellname )
mye <- ScaleData(mye, verbose = FALSE)
mye <- FindVariableFeatures(mye, selection.method = "vst", nfeatures = length(genes))
mye <- RunPCA(mye, npcs = 19, verbose = FALSE)

#mye_dat = GetAssayData(mye, assay = 'RNA', slot = 'data')

mye_dat <- GetAssayData(mye)

mye_dat = t(mye_dat)  


day = matrix('0', length(mye$stim), 1)
id = mye$stim == 'wk2'
day[id, 1] = 'wk2'
id = mye$stim == 'wk4'
day[id, 1] = 'wk4'

mye_PHATE = phate(mye_dat)

mye_PHATE <- phate(mye_dat, knn = 10, decay = 40, t = 10, init = mye_PHATE)


gene = 'Cd177'

ggplot(mye_PHATE) +
		geom_point(aes(PHATE1, PHATE2, color=mye_dat[, gene])) +
		labs(color= gene) +
		scale_color_viridis( option = "B")

ggsave(paste('sox2_', gene, '_phate.png', sep = ''), width=5, height=5)


ggplot(mye_PHATE) +
		geom_point(aes(PHATE1, PHATE2, colour = factor(day))) +
		labs(color= "wks") +
		scale_color_viridis( option = "B")

ggsave('sox2_wk_phate_06_16_2020.png', width=5, height=5)


pdf('Phate_wk.pdf', width = 30, height = 10)
par(mfrow = c(1, 3))
id = day == 0
plot(mye_PHATE$embedding[id, ], main = 'Control')
id = day == 2
plot(mye_PHATE$embedding[id, ],  main = 'Wk2')
id = day == 4
plot(mye_PHATE$embedding[id, ],  main = 'Wk4')
dev.off()

load('my_PHATE_06_03_2020.rda')
gene_PHATE = c( 'Cxcr5',  'Cd69', 'Aim2',  'Irf5', 'Irf1', 'Irf3', 'Lgals9', 'Ly6e', 'Nos2',
   'Il6', 'Tnf',  'Isg15', 'Gsdmd',  'Cd8a', 'Cd40', 'Cd80', 'Cd86',  'Il10', 'Tgfb1',
    'Tmem173', 'Batf2', 'Ly6g', 'Cd274', 'Cxcl9', 'Arg1', 'Mrc1')
 
id = match(colnames(ctr_dat), colnames(dat))

plot(ctr_dat[ 'Gata2', ], dat[ 'Gata2', id])


#Fig 5 E-J


gene_PHATE = c('Il1a',  'Slc2a1', 'Stat1' 'Hif1a')

i = 0

i = i + 1
gene = gene_PHATE[i]

ggplot(mye_PHATE) +
		geom_point(aes(PHATE1, PHATE2, color=mye_dat[, gene])) +
		labs(color= gene) +
		scale_color_viridis( option = "B")

ggsave(paste('sox2_', gene, '_phate_07_29_2020.pdf', sep = ''), width=10, height=10)





###################################
### Bulk 658WG             ########
###################################

# building index
salmon index -t Mus.GRCh38.cdna.all.fa.gz -i ms_index

# quantify the  samples

salmon quant -i ms_index -l A \
    -1 _W658G/658WG9/658WG9_1.fastq.gz \
    -2 658_WG/658WG9/658WG9_2.fastq.gz \
    -p 8 -o WG_sample_9_quants_nogccorrection /
library("tximport")
library("readr")
library("tximportData")
library(dplyr)
library("ggplot2")

library(vsn)
dir <- system.file("extdata", package="tximportData")
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
samples$condition <- factor(rep(c("A","B"),each=3))
rownames(samples) <- samples$run
samples[,c("pop","center","run","condition")]

## ----txiFiles------------------------------------------------------------
files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
names(files) <- samples$run
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
## ----tximport, results="hide"--------------------------------------------
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
## ----txi2dds, results="hide"---------------------------------------------
library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi,  colData = samples,
        design = ~ condition)

log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)

dds <- estimateSizeFactors(dds)

df = bind_rows(
          as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
                         mutate(transformation = "log2(x + 1)"),
          as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
          as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog")
  )

colnames(df)[1:2] <- c("x", "y")  
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
          coord_fixed() + facet_grid( . ~ transformation)  

###################################
### Bulk 5211HT             #######
###################################     
---
title: "5211HT"
author: "Zhaoheng Li"
date: "3/7/2022"
output: html_document
---
```{r}
library(DESeq2)
library(dplyr)
library(ggplot2)
library(tibble)
library(factoextra)
```


```{r}
data <- read.delim("gene_expected_count.annot.txt")


WT_IL1R = data[-(55422:55426),paste("Sample_5211.HT.",1:10,sep='')]
rownames(WT_IL1R)=data$gene_id[-(55422:55426)]
WT_KO = data[-(55422:55426),paste("Sample_5211.HT.",11:20,sep='')]
rownames(WT_KO)=data$gene_id[-(55422:55426)]
```

```{r IL1R}

coldata = data.frame(condition = factor(c(rep("WT",5),rep("IL1R",5))))
WT_IL1R <- DESeqDataSetFromMatrix(countData = WT_IL1R,
                              colData = coldata,
                              design = ~ condition)
WT_IL1R
WT_IL1R$condition <- relevel(WT_IL1R$condition, ref = "WT")
WT_IL1R <- DESeq(WT_IL1R)
res.WT_IL1R <- results(WT_IL1R)

```


```{r WT_KO}
coldata = data.frame(condition = factor(c(rep("WT",5),rep("KO",5))))
WT_KO <- DESeqDataSetFromMatrix(countData = WT_KO,
                              colData = coldata,
                              design = ~ condition)
WT_KO
WT_KO$condition <- relevel(WT_KO$condition, ref = "WT")
WT_KO <- DESeq(WT_KO)
res.WT_KO <- results(WT_KO)


```


```{r}
res.WT_IL1R=data.frame(res.WT_IL1R)%>%arrange(log2FoldChange)
res.WT_KO=data.frame(res.WT_KO)%>%arrange(log2FoldChange)

res.WT_IL1R=res.WT_IL1R%>%rownames_to_column()%>%left_join(data%>%select(gene_id,external_gene_name,description),
                                                           by = c("rowname"="gene_id"))
res.WT_KO=res.WT_KO%>%rownames_to_column()%>%left_join(data%>%select(gene_id,external_gene_name,description),
                                                       by = c("rowname"="gene_id"))
head(res.WT_IL1R)
```
```{r}
res.WT_IL1R%>%filter(external_gene_name=='Il1r1')
res.WT_KO%>%filter(external_gene_name=='Il1r1')
```


```{r}
write.csv(res.WT_IL1R,"DEres.1-10.csv")
write.csv(res.WT_KO,"DEres.11-20.csv")
```



```{r}
WT_KO <- estimateSizeFactors(WT_KO)
WT_IL1R=estimateSizeFactors(WT_IL1R)
colData(WT_KO)
sizeFactors(WT_KO)
assay(WT_KO,normalized=TRUE)

pdf("hie_cluster.pdf")
d <- dist(t(assay(WT_KO,normalized=TRUE)), method = "euclidean") # distance matrix
fit <- hclust(d, method="complete")
plot(fit) # display dendogram
d <- dist(t(assay(WT_IL1R,normalized=TRUE)), method = "euclidean") # distance matrix
fit <- hclust(d, method="complete")
plot(fit) # display dendogram
dev.off()
```


```{r}
pdf("MOR.Norm_hie_cluster.pdf")
d <- dist(t(t(t(assay(WT_KO))/sizeFactors(WT_KO))), method = "euclidean") # distance matrix
fit <- hclust(d, method="complete")
plot(fit) # display dendogram
d <- dist(t(t(t(assay(WT_IL1R))/sizeFactors(WT_IL1R))), method = "euclidean") # distance matrix
fit <- hclust(d, method="complete")
plot(fit) # display dendogram
dev.off()
```

```{r}

WT_IL1R.MOR.pca = prcomp(t(counts(WT_IL1R,normalized=TRUE)))
fviz_pca_ind(WT_IL1R.MOR.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
d <- dist(WT_IL1R.MOR.pca$x[,1:2], method = "euclidean") # distance matrix
fit <- hclust(d, method="complete")
plot(fit) # display dendogram

WT_KO.MOR.pca = prcomp(t(counts(WT_KO,normalized=TRUE)))
fviz_pca_ind(WT_KO.MOR.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
d <- dist(WT_KO.MOR.pca$x[,1:2], method = "euclidean") # distance matrix
fit <- hclust(d, method="complete")
plot(fit) # display dendogram
```



# tpm

```{r}
tpm.data <- read.delim("/Users/lizhaoheng/Desktop/5211HT/gene_TPM.annot.txt")



tpm.WT_IL1R = tpm.data[-(55422:55426),paste("Sample_5211.HT.",1:10,sep='')]
rownames(tpm.WT_IL1R)=data$gene_id[-(55422:55426)]
tpm.WT_KO = tpm.data[-(55422:55426),paste("Sample_5211.HT.",11:20,sep='')]
rownames(tpm.WT_KO)=data$gene_id[-(55422:55426)]

```

```{r}

tpm.data%>%filter(external_gene_name=='Il1r1')%>%select(paste("Sample_5211.HT.",1:10,sep=''))

df = data.frame(IL1R =t(tpm.data%>%filter(external_gene_name=='Il1r1')%>%select(paste("Sample_5211.HT.",1:10,sep=''))),condition = factor(rep(c("WT", "IL1R"), each = 5)))
pdf("tpm.IL1R.pdf")
df%>%ggplot(aes(x=condition,y=IL1R,color = condition))+
  geom_boxplot()+
  theme_classic()+
  geom_jitter()
dev.off()
```

```{r}
colnames(scl.tpm.WT_IL1R) = paste("HT",1:10,sep='')
round(scl.tpm.WT_IL1R[1:20,],3)
```

```{r}


sds <- apply(tpm.WT_IL1R,1,sd)
scl.tpm.WT_IL1R = tpm.WT_IL1R[-which(sds==0),]
means <- apply(scl.tpm.WT_IL1R,1,mean)
sds <- apply(scl.tpm.WT_IL1R,1,sd)
scl.tpm.WT_IL1R <- t(scale(t(scl.tpm.WT_IL1R),center=means,scale=sds))




sds <- apply(tpm.WT_KO,1,sd)
scl.tpm.WT_KO  = tpm.WT_KO[-which(sds==0),]
means <- apply(scl.tpm.WT_KO,1,mean)
sds <- apply(scl.tpm.WT_KO,1,sd)
scl.tpm.WT_KO <- t(scale(t(scl.tpm.WT_KO),center=means,scale=sds))




pdf("scl_tpm_hieClu.pdf")
d <- dist(t(scl.tpm.WT_IL1R), method = "euclidean") # distance matrix
fit <- hclust(d, method="complete")
plot(fit, 'Pearson correlation') # display dendogram
d <- dist(t(scl.tpm.WT_KO), method = "euclidean") # distance matrix
fit <- hclust(d, method="complete")
plot(fit, 'Pearson correlation') # display dendogram
#dev.off()


c <- cor((scl.tpm.WT_IL1R), method="pearson") 
d <- as.dist(1-c)
fit <- hclust(d, method="complete")
plot(fit, 'Pearson correlation') # display dendogram


c <- cor((scl.tpm.WT_KO), method="pearson") 
d <- as.dist(1-c)
fit <- hclust(d, method="complete")
plot(fit, main = 'Pearson correlation') # display dendogram
dev.off()
```




## pca

### unscaled

```{r}
tpm.WT_IL1R.pca = prcomp(t(tpm.WT_IL1R))
fviz_pca_ind(tpm.WT_IL1R.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
d <- dist(tpm.WT_IL1R.pca$x[,1:2], method = "euclidean") # distance matrix
fit <- hclust(d, method="complete")
plot(fit) # display dendogram

```

```{r}
tpm.WT_KO.pca = prcomp(t(tpm.WT_KO))
fviz_pca_ind(tpm.WT_KO.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
d <- dist(tpm.WT_KO.pca$x[,1:5], method = "euclidean") # distance matrix
fit <- hclust(d, method="complete")
plot(fit) # display dendogram

```

### scaled

```{r}
#pdf("scl.tpm.WT_IL1R.pca.pdf")
tpm.WT_IL1R.pca = prcomp(t(scl.tpm.WT_IL1R))
fviz_pca_ind(tpm.WT_IL1R.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
d <- dist(tpm.WT_IL1R.pca$x[,1:5], method = "euclidean") # distance matrix
fit <- hclust(d, method="complete")
plot(fit) # display dendogram
#dev.off()
```

```{r}
#pdf("scl.tpm.WT_KO.pca.pdf")
tpm.WT_KO.pca = prcomp(t(scl.tpm.WT_KO))
fviz_pca_ind(tpm.WT_KO.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
d <- dist(tpm.WT_KO.pca$x[,1:5], method = "euclidean") # distance matrix
fit <- hclust(d, method="complete")
plot(fit) # display dendogram
#dev.off()
```
```{r}
library(FactoMineR)
# Compute PCA with ncp = 3
res.pca <- PCA(t(scl.tpm.WT_KO), ncp = 9)
# Compute hierarchical clustering on principal components
res.hcpc <- HCPC(res.pca, graph = FALSE)
```

```{r}
fviz_dend(res.hcpc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
          )
fviz_cluster(res.hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
             )
```

