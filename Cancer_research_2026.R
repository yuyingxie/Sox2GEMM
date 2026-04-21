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
        