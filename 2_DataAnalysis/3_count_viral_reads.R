library(Seurat)
library(tidyverse)
library(readxl)
#---------------------------------------


#### COUNTVIRAL READS ####

## read GTF

## read in GTF file for gene feature positions
gtf2 = read.table("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/MPXV_NC_063383.1_annotated.gtf")
pbmc.int = readRDS( file = "C:/Users/postmusd/Documents/MPXV_PBMC/hg38_10x/pbmc.2.ct_mpox.rds")


## read in viral read counts by barcode

## lib 002 and 004 # infected conditions
mpx_002 <- read.table(file = "C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/viral_reads/run1_combinedref/002_viral_posMap_BC.txt")
mpx_002$V1 <- gsub("-1", "-2", mpx_002$V1)

mpx_004 <- read.table(file = "C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/viral_reads/run1_combinedref/004_viral_posMap_BC.txt")
mpx_004$V1 <- gsub("-1", "-4", mpx_004$V1)

# combine
mpx_read <- rbind(mpx_002, mpx_004)
head(mpx_read)
## assign names to columns
##############################################
################################################ fix if needed

names(mpx_read) = c("bc", "umi", "tag", "cig", "start_pos")
#names(mpx_read) = c("bc", "umi", "tag", "start_pos")
head(mpx_read)
mpx_read$end_pos = mpx_read$start_pos + 89 # extend read to endpoint by adding 89bp from leftmost map position
mpx_read$gene_start = "unmapped" ## placeholder
mpx_read$gene_end = "unmapped"

mpx_read



#############################################

#### assign gene based on read position ####
print("Counting reads!")


## define function to assign gene feature overlap with left and
## right ends of read based on positions in GTF
## input = will be replaced by read start position

getgene = function(input){
  
  # check which, if any, genes the left end of the read overlaps with
  gene_start = gtf2[with(gtf2, V4 <= input & V5 >= input),]$V16
  # if left end is in region of overlapping genes, flag as in gene overlap
  if(length(gene_start) >1 ){
    gene_start = "overlap"
  }
  ## add read end point
  pos_end = input + 89
  # check which, if any, genes the right end of the read overlaps with
  gene_end = gtf2[with(gtf2, V4 <= pos_end & V5 >= pos_end),]$V16
  # if left end is in region of overlapping genes, flag as in gene overlap
  if(length(gene_end) >1 ){
    gene_end = "overlap"
  }
  #create a result with left and right gene overlaps marked by comma
  # one or more can be empty ""
  return(paste0(gene_start, ",", gene_end))

  
}


## map getgene function to every read start position in mpx read dataframe
mpx_read$out = map_chr(.x = mpx_read$start_pos,
                       .f = getgene, .progress = T)


check = mpx_read[12369,]
check

#### create matrix ####

## filter mapped reads and refine assignment to specific genes

mpx_read2 = mpx_read %>%
  filter(!(grepl(pattern = "overlap", x = out))) %>% ## filter out any reads mapping to overlapping gene positions
  separate(out, sep = ",", into = c("gene_start", "gene_end")) %>%
  filter(!(gene_start =="" & gene_end == "")) %>% ## remove unmapped reads ## both overlap sides are empty
  mutate(gene_assign = ifelse(gene_start == "" & gene_end != "", gene_end, ## overlap feature on right only
                              ifelse(gene_start != "" & gene_end == "", gene_start, ## overlap feature on left only
                                     ifelse(gene_start != "" & gene_end != "" & gene_start == gene_end, gene_start, "multimap")))) %>%
  ## remove any reads where read start and read end overlap different genes
  filter(gene_assign != "multimap") %>% 
  ## count reads by bc+umi+ gene they've been assigned to
  group_by(bc, umi, gene_assign) %>%
  tally() %>%
  group_by(bc, umi) %>%
  ## if same read maps to two/more different gene features, keep only 
  ## gene feature with highest number of mappings
  filter(n == max(n)) %>%
  ## filter out any ties -- read with equal number of mappings to two/more different features 
  ## cannot be confidently assigned to a gene and must be removed
  group_by(bc, umi) %>%
  filter(!(n() > 1))

## convert to sparse matrix and export

mpx_mat = mpx_read2 %>%
  #filter(gene != "unmapped") %>%
  group_by(bc, gene_assign) %>%
  tally() %>%
  pivot_wider(names_from = bc, values_from = n ) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(var = "gene_assign") %>%
  #as.matrix() %>%
  as.sparse()

saveRDS(object = mpx_mat, file = "C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/viral_reads/run1_combinedref/mpx_readmtx_allhits.rds")


## subset to include only cells in pbmc object
mpx_mat_sub = mpx_mat[,colnames(mpx_mat) %in% colnames(pbmc.int)]
saveRDS(object = mpx_mat_sub, file = "C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/viral_reads/run1_combinedref/mpx_readmtx_cellsinobject.rds")

#### Create Seurat Assay Object for Mpx reads to add to main object #####

## all cells in pbmc object that do not have assigned viral reads
all_neg_cells = colnames(pbmc.int)[!(colnames(pbmc.int) %in% colnames(mpx_mat_sub))]
## all cells in pbmc object
all_cells = colnames(pbmc.int)
## all cells with assigned mpox reads
all_pos_cells = colnames(mpx_mat_sub)

#temp_df = data.frame(row.names = row.names(all_pos_cells), )

## create empty df with barcodes of all other cells (mpx-) to append to
temp_mat = matrix(0,
       nrow = nrow(mpx_mat_sub),
       ncol = length(all_neg_cells))

## assign mpx gene names as rownames
row.names(temp_mat) = row.names(mpx_mat_sub)
## assign cell bcs as colnames
colnames(temp_mat) = all_neg_cells
## convert to sparse
temp_mat = as.sparse(temp_mat)
## combine matrices and create Seurat assay object
all_cell_mat = cbind(temp_mat, mpx_mat_sub)
mpx.assay = CreateAssayObject(all_cell_mat)
## add assay
pbmc.int[["MPXV"]] = mpx.assay


#### Preprocess MPXV assay ####
DefaultAssay(pbmc.int) = "MPXV"

pbmc.int = NormalizeData(pbmc.int)
pbmc.int <- FindVariableFeatures(pbmc.int)
pbmc.int <- ScaleData(pbmc.int)

## save object with viral assay 

saveRDS(object = pbmc.int, file = "C:/Users/postmusd/Documents/MPXV_PBMC/hg38_10x/pbmc.3.ct_mpox_mpoxassay.rds")

pbmc.int

head(pbmc.int@meta.data)
Idents(pbmc.int) = "treatment"
unique(pbmc.int$treatment)
FindMarkers(pbmc.int, ident.1 = "Mock_5dpi", ident.2 = "MPXV_5dpi", min.pct = 0)
?FindMarkers()
FeaturePlot()
