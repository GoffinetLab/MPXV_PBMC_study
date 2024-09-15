library(Seurat)
library(tidyverse)
library(readxl)

#BiocManager::install("Seurat")
#---------------------------------------


pbmc.int = readRDS(file = "C:/Users/postmusd/Documents/MPXV_PBMC/hg38_10x/pbmc.3.ct_mpox_mpoxassay.rds")
head(pbmc.int@meta.data)
unique(pbmc.int$infection)
DefaultAssay(pbmc.int) = "MPXV"

mpx_genes = rownames(pbmc.int@assays$MPXV@data)

mpx_genes

test = order(mpx_genes)
mpx_genes2 = mpx_genes[test]
mpx_genes2
## get avg expression of mpxgenes
## subset to only infected condi
Idents(pbmc.int) = "infection"
pbmc.sub = subset(pbmc.int, idents = c("MPXV"))

Idents(pbmc.sub) = "treat_ct"

avg <- AverageExpression(pbmc.sub, features = mpx_genes, return.seurat = T, assays = "MPXV") 
avg2 <- AverageExpression(pbmc.sub, features = mpx_genes, return.seurat = F, assays = "MPXV") 

avg2_df = as.data.frame(avg2$MPXV)
avg2_df = cbind(" "=rownames(avg2_df), avg2_df)
avg2_df

writexl::write_xlsx(x = avg2_df,path = "C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/viral/MPXV_gene_byCellType_Timepoint.xlsx",)

## reorder mpx genes in spatial order
#avg$MPXV = avg$MPXV[mpx_genes2,]
#avg@assays$MPXV@scale.data = avg@assays$MPXV@scale.data[mpx_genes2,]
#avg@assays$MPXV@scale.data

avg@assays$MPXV
h <- DoHeatmap(avg, assay="MPXV", size=5, features=mpx_genes, label=T, group.bar=F,
               draw.lines = F, hjust = 0)+ 
  scale_fill_gradientn(colors = c("blue", "white","red"), na.value = "white") +
NoLegend()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))

## reorder levels of genes within ggplot2 object
h$data$Feature = factor(h$data$Feature, levels = rev(mpx_genes2))
h


ggsave(plot =h,
       filename = "C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/viral/hmap_allMPXgenes.pdf",
       width = 20, height = 25)


fetch = FetchData(pbmc.sub, vars = c(mpx_genes, "treat_ct"), slot = "data")
fetch2 = pivot_longer(fetch, cols = )
unique(fetch$mpxv.OPG047)
?pivot_longer()
