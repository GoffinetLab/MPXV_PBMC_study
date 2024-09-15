library(Seurat)
library(tidyverse)
library(readxl)
library(writexl)
#---------------------------------------

## read in ct and mpox read annotated object

pbmc.int <- readRDS("C:/Users/postmusd/Documents/MPXV_PBMC/hg38_10x/pbmc.2.ct_mpox.rds")
pbmc.int$mpx_n10 <- ifelse(pbmc.int$mpx_reads_raw > 10, "pos", "neg")

## perform summary statistics

meta <- pbmc.int@meta.data

## % MPOX positive by CT by treatment #####

## all mpox positive cells

df1 <- meta %>%
  group_by(treatment, celltype) %>%
  summarise(pos_cells = sum(mpx == "pos"), 
            tot_cells = n(),
            perc_mpox_pos = pos_cells / tot_cells) 

bplot_mpx <- ggplot(data=df1, aes(celltype, perc_mpox_pos, fill = treatment))+
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  ggtitle("MPXV Positive Fraction")+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))


ggsave(plot = bplot_mpx,
       filename = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/summary_stats", 
                         "/barplot_mpxpositive.png"),
       width = 8, height = 3)


## mpox if only consider cells with v reads > 10

df2 <- meta %>%
  group_by(treatment, celltype) %>%
  summarise(pos_cells = sum(mpx_n10 == "pos"), 
            tot_cells = n(),
            perc_mpox_pos = pos_cells / tot_cells) 



bplot_mpx_n10 <- ggplot(data=df2, aes(celltype, perc_mpox_pos, fill = treatment))+
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  ggtitle("MPXV Positive Fraction (min 10 reads / cell)")+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))
bplot_mpx_n10

ggsave(plot = bplot_mpx_n10,
       filename = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/summary_stats", 
                         "/barplot_mpxpositive_min10reads.png"),
       width = 8, height = 3)


## CELLTYPE representation between treatments #####

## 

# get num of all cells of a celltype
ct <- meta %>%
  group_by(celltype) %>%
  summarise(all_cells = n())
# get num of cells of a celltype within treatment
ct2 <- meta %>%
  group_by(treatment, celltype) %>%
  summarise(n_cell = n()) %>%
  plyr::join(., ct2) %>%
  mutate(perc_ct = n_cell / all_cells)

ggplot(data = ct2, ())

tail(ct2, n = 20)

ggplot(ct2, aes(fill=treatment, y=perc_ct, x=celltype)) + 
  geom_bar(position="fill", stat="identity")

####################################

## get % of celltypes of all cells in lib

# get num of all cells in a treatment
ct <- meta %>%
  group_by(treatment) %>%
  summarise(all_cells = n())
# get num of cells of a celltype within treatment
ct2 <- meta %>%
  group_by(treatment, celltype) %>%
  summarise(n_cell = n()) %>%
  plyr::join(., ct2) %>%
  mutate(perc_ct = n_cell / all_cells)
ct2
ct_frac_allcells <- ggplot(ct2, aes(fill=treatment, y=perc_ct, x=celltype)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  ggtitle("Fraction of celltypes as % of all cells in treatment")+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))

ggsave(plot = ct_frac_allcells,
       filename = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/summary_stats", 
                         "/barplot_fractionct_ofallcellsintreat.png"),
       width = 8, height = 3)


### DEG test
pbmc.int$treat_ct <- paste0(pbmc.int$treatment, ".", pbmc.int$celltype)
pbmc.int$infect_ct <- paste0(pbmc.int$infection, ".", pbmc.int$celltype)
Idents(pbmc.int) <- "infect_ct"

unique(pbmc.int$infect_ct)
d1 <- FindMarkers(pbmc.int, ident.1 = "MPXV.Monocyte", ident.2 = "Mock.Monocyte")
d1$gene <- rownames(d1)
d1
writexl::write_xlsx(d1, path = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/summary_stats", 
                                      "/DEG_MPXV.Monos_v_Mock.Monos.xlsx"))

d2 <- FindMarkers(pbmc.int, ident.1 = "MPXV.CD4_T_Reg", ident.2 = "Mock.CD4_T_Reg")
d2$gene <- rownames(d2)
d2
writexl::write_xlsx(d2, path = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/summary_stats", 
                                      "/DEG_MPXV.Cd4Treg_v_Mock.Cd4Treg.xlsx"))

d3 <- FindMarkers(pbmc.int, ident.1 = "MPXV.NK_cycling", ident.2 = "Mock.NK_cycling")
d3$gene <- rownames(d3)
d3
writexl::write_xlsx(d2, path = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/summary_stats", 
                                      "/DEG_MPXV.NKcycling_v_Mock.NKcycling.xlsx"))

head(pbmc.int@meta.data)


Idents(pbmc.int) <- "treat_ct"
d4 <- FindMarkers(pbmc.int, ident.1 = "MPXV_3dpi.Monocyte", ident.2 = "Mock_3dpi.Monocyte")
d4$gene <- rownames(d4)
writexl::write_xlsx(d4, path = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/summary_stats", 
                                      "/DEG.3dpi_Mono_MPXV_v_Mock.xlsx"))

d5 <- FindMarkers(pbmc.int, ident.1 = "MPXV_5dpi.Monocyte", ident.2 = "Mock_5dpi.Monocyte")
d5$gene <- rownames(d5)
writexl::write_xlsx(d5, path = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/summary_stats", 
                                      "/DEG.5dpi_Mono_MPXV_v_Mock.xlsx"))

markers <- c("STAT1")
avg <- AverageExpression(pbmc.int, features = markers, return.seurat = T, assays = "RNA") 


h <- DoHeatmap(avg, assay="RNA", size=5, features=markers, label=T, group.bar=F,
               draw.lines = F, hjust = 0)+ 
  scale_fill_gradientn(colors = c("blue", "white","red"), na.value = "white") +
  #NoLegend()+
  theme_bw()#+

h

DoHeatmap(pbmc.int, features = c("STAT1"))

d5

d2$gene <- rownames(d2)
d2
