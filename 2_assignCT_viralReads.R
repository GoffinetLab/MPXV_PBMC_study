library(Seurat)
library(tidyverse)
library(readxl)
#---------------------------------------

pbmc.int <- readRDS("C:/Users/postmusd/Documents/MPXV_PBMC/hg38_10x/pbmc.1.rds")

umap.path <- "C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/umap"
fplot.path <- "C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/featplot"
hmap.path <- "C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/hmap"
pbmc.int$integrated_snn_res.0.5
## viral figures folder

#### ASSIGN CELLTYPES ####

# create umap with clusters

Idents(pbmc.int) <- "integrated_snn_res.0.5"

umap_clust0.5 <- DimPlot(pbmc.int, label = T) +
  theme(legend.position = "none")

ggsave(filename = paste0(umap.path, "/umap_clusters_res0.5.png"),
       plot = umap_clust0.5,
       width=6.4, height = 4.01)



## create heatmap with marker genes

markers <- c("CD3D","CD8A", "SELL", "CREM", "S100A4", "IL2RA", "FOXP3", "IKZF2", 
             "GNLY", "NKG7", "KLRB1", "TOP2A", 
             "CD79A", "MS4A7", "MS4A1", "TCL1A","IGHA1", "IGHG1", 
             "FCGR3A", "FCER1A", "CST3", "CCL2", "LYZ", "CD14","LILRA4", "CD1C", "CD34", 
             "PPBP", 
             "NME1", "GZMA", "GZMB", "GZMK")

DefaultAssay(pbmc.int) <- "RNA"


## get avg expression of markers
avg <- AverageExpression(pbmc.int, features = markers, return.seurat = T, assays = "RNA") 


h <- DoHeatmap(avg, assay="RNA", size=5, features=markers, label=T, group.bar=F,
               draw.lines = F, hjust = 0)+ 
  scale_fill_gradientn(colors = c("blue", "white","red"), na.value = "white") +
  #NoLegend()+
  theme_bw()#+

h


ggsave(filename = paste0(hmap.path, "/hmap_CTmarks_byclust0.5.png"),
       plot = h,
       width=6.4, height = 7)

## read CT annotation list

ct_assign <- readxl::read_excel("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/CT_cluster_assignxlsx.xlsx")
ct_assign$cluster <- as.character(ct_assign$cluster)
# create named vector
ct_assign_ls <- deframe(ct_assign)
ct_assign_ls
# assign IDs based on named vector
# cts assigned based on 0.5 clustering
pbmc.int$celltype <- pbmc.int$integrated_snn_res.0.5
pbmc.int$celltype <- plyr::revalue(pbmc.int$celltype, replace = ct_assign_ls)
head(pbmc.int@meta.data)
Idents(pbmc.int) <- "celltype"

umap_celltype <- DimPlot(pbmc.int, label = T) +
  theme(legend.position = "none")

pbmc.int$treat_ct <- paste0(pbmc.int$treatment, ".", pbmc.int$celltype)

ggsave(filename = paste0(umap.path, "/umap_celltype.png"), 
       plot = umap_celltype,
       width=6.4, height = 4.01)

#### ASSIGN VIRAL READS ####

## read in viral read counts by barcode

## lib 002 and 004 # infected conditions
mpx_002 <- read.table(file = "C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/viral_reads/run1_combinedref/002_viral_BC.txt")
mpx_002$V1 <- gsub("-1", "-2", mpx_002$V1)
mpx_002$lib <- 2
mpx_004 <- read.table(file = "C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/viral_reads/run1_combinedref/004_viral_BC.txt")
mpx_004$V1 <- gsub("-1", "-4", mpx_004$V1)
mpx_004$lib <- 4
# combine
mpx_all <- rbind(mpx_002, mpx_004)

# count reads by barcode
mpx_rcount <- mpx_all %>%
  group_by(V1,V2)%>% ## group by varcode and umi
  ## deduplicate so that only one bc + umi combo remains
  slice(n=1) %>%
  ## count number of reads by cell(barcode)
  group_by(V1) %>%
  summarise(mpx_reads_raw=n()) %>%
  rename(barcode=V1)
mpx_rcount


# mark cells as mpox pos or neg 
pbmc.int$mpx <- ifelse(pbmc.int$barcode %in% mpx_all$V1, "pos", "neg")

table(pbmc.int$celltype, pbmc.int$mpx)
mpx_all

Idents(pbmc.int) <- "mpx"

## umap
umap_mpox_pos.neg <- DimPlot(pbmc.int, order = "pos", cols = c("grey", "red"))

ggsave(filename = paste0(umap.path, "/umap_mpox_posneg.png"), 
       plot = umap_mpox_pos.neg,
       width=6.4, height = 4.01)

# # add counts to object metadata

meta <- pbmc.int@meta.data
meta <- meta %>%
  plyr::join(., mpx_rcount) %>%
  replace_na(list(mpx_reads_raw = 0)) ## change NAs to 0 reads
# add mpx_reads_raw column to object

pbmc.int$mpx_reads_raw <- meta$mpx_reads_raw

# normalise mpox reads
#   div by total cell reads
#   multiply by 10 000 (scaling factor)
#   take log e

pbmc.int$mpx_reads_norm <- pbmc.int$mpx_reads_raw / pbmc.int$nCount_RNA
pbmc.int$mpx_reads_norm <- pbmc.int$mpx_reads_norm * 10000
pbmc.int$mpx_reads_norm <- log1p(pbmc.int$mpx_reads_norm)

head(pbmc.int@meta.data)
fplot_mpox_raw <- FeaturePlot(pbmc.int, features = c("mpx_reads_raw") ) +
  ggtitle("")

## export Featplots with MPOX expression levels


fplot_mpox_norm <- FeaturePlot(pbmc.int, features = c("mpx_reads_norm")) +
  ggtitle("")

ggsave(filename = paste0(fplot.path, "/featplot_mpox_rawcount.png"), 
       plot = fplot_mpox_raw,
       width=6.4, height = 4.01)
ggsave(filename = paste0(fplot.path, "/featplot_mpox_normcount.png"), 
       plot = fplot_mpox_norm,
       width=6.4, height = 4.01)

fplot_mpox_norm_split <- FeaturePlot(pbmc.int, features = c("mpx_reads_norm"), split.by = "treatment") 

ggsave(filename = paste0(fplot.path, "/featplot_mpox_normcount_bytreat.png"), 
       plot = fplot_mpox_raw,
       width=20, height = 4.01)

#### save annotated object ####

saveRDS(object = pbmc.int, file = "C:/Users/postmusd/Documents/MPXV_PBMC/hg38_10x/pbmc.2.ct_mpox.rds")
