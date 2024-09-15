
.libPaths("C:/Program Files/R/R-4.1.0/library")
library(S4Vectors)

library(org.Hs.eg.db)
library(Seurat)
library(writexl)
library(tidyverse)
library(readxl)
library(clusterProfiler, lib.loc = "C:/Program Files/R/R-4.1.0/library")
library(tidygraph)
library(ReactomePA)


#---------------------------------------
#### DEFINE FUNCTIONS ####

## GO Enrichment
run_go = function(input, name, dpi=NULL){
  
  go_res_df = data.frame()
  
  for (ct in unique(pbmc.int$celltype)){
    
    go_input <- input[[ct]] %>%
      #assign ENTREZID ID from symbol
      mutate(ensembl_id = mapIds(org.Hs.eg.db, keys=.$gene, column="ENTREZID", keytype="SYMBOL", multiVals="first")) %>%
      select(c(ensembl_id, avg_log2FC)) %>%
      na.omit () %>% # remove unmapped
      arrange(desc(avg_log2FC)) %>%
      deframe()
    ## run GO enrichment
    goenr <- gseGO(geneList     = go_input,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "BP",
                   minGSSize    = 20,
                   maxGSSize    = 500,
                   pvalueCutoff = 1,
                   verbose      = T)
    
    ## get result
    goenr_tidy <-goenr@result 
    
    if (nrow(goenr_tidy) != 0){
      ## clean and determine GeneRatio
      gene_count <- goenr_tidy %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
      merged_res <- left_join(goenr_tidy, gene_count, by = "ID") %>% 
        mutate(GeneRatio = count/setSize) %>%
        mutate(zscore = scale(NES)) %>%
        mutate(group = ct)
      
      if (!is.null(dpi)){
        merged_res$dpi = dpi
      }
      
      ## add goRes to dataframe
      go_res_df = rbind(go_res_df, merged_res)
      
    }
    
    
  }
  
  ## export to global environment
  assign(x = paste0("gse_GOBP.", name), value = go_res_df, envir = .GlobalEnv)
  saveRDS(object = go_res_df, file = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/geneset_enrich/","gse_GOBP.",
                                            name, ".rds"))
}

## function to run for single DEG res
run_go_singledf = function(input, name, contrast){
  
  go_res_df = data.frame()
    
    go_input <- input %>%
      #assign ENTREZID ID from symbol
      mutate(ensembl_id = mapIds(org.Hs.eg.db, keys=.$gene, column="ENTREZID", keytype="SYMBOL", multiVals="first")) %>%
      select(c(ensembl_id, avg_log2FC)) %>%
      na.omit () %>% # remove unmapped
      arrange(desc(avg_log2FC)) %>%
      deframe()
    ## run GO enrichment
    goenr <- gseGO(geneList     = go_input,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "BP",
                   minGSSize    = 20,
                   maxGSSize    = 500,
                   pvalueCutoff = 1,
                   verbose      = T)
    
    ## get result
    goenr_tidy <-goenr@result 
    
      gene_count <- goenr_tidy %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
      merged_res <- left_join(goenr_tidy, gene_count, by = "ID") %>% 
        mutate(GeneRatio = count/setSize) %>%
        mutate(zscore = scale(NES)) %>%
        mutate(group = contrast)
      
      
      ## add goRes to dataframe
      go_res_df = rbind(go_res_df, merged_res)
      
  
  ## export to global environment
  assign(x = paste0("gse_GOBP.", name), value = go_res_df, envir = .GlobalEnv)
  saveRDS(object = go_res_df, file = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/geneset_enrich/","gse_GOBP.",
                                            name, ".rds"))
}


## REACTOME Enrichment

run_reactome = function(input, name, dpi=NULL){
  
  go_res_df = data.frame()
  
  for (ct in unique(pbmc.int$celltype)){
    
    go_input <- input[[ct]] %>%
      #assign ENTREZID ID from symbol
      mutate(ensembl_id = mapIds(org.Hs.eg.db, keys=.$gene, column="ENTREZID", keytype="SYMBOL", multiVals="first")) %>%
      select(c(ensembl_id, avg_log2FC)) %>%
      #na.omit () %>% # remove unmapped
      arrange(desc(avg_log2FC)) %>%
      deframe()
    
    
    if(length(go_input) > 5){ ## avoid problem of small no. of  ENTREZIDs not being in list
      ## run GO enrichment
      print(ct)
      reactenr <- gsePathway(geneList= go_input,
                             maxGSSize    = 500,
                             pvalueCutoff = 1,
                             pAdjustMethod = "BH",
                             verbose      = T)
      
      ## get result
      reactenr_tidy <-reactenr@result 
      
      if (nrow(reactenr_tidy) != 0){
        ## clean and determine GeneRatio
        gene_count <- reactenr_tidy %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
        merged_res <- left_join(reactenr_tidy, gene_count, by = "ID") %>% 
          mutate(GeneRatio = count/setSize) %>%
          mutate(zscore = scale(NES)) %>%
          mutate(group = ct)
        
        if (!is.null(dpi)){
          merged_res$dpi = dpi
        }
        
        ## add goRes to dataframe
        go_res_df = rbind(go_res_df, merged_res)
        
      }
      
    }
    
  }
  
  ## export to global environment
  assign(x = paste0("gse_REACT.", name), value = go_res_df, envir = .GlobalEnv)
  saveRDS(object = go_res_df, file = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/geneset_enrich/","gse_REACTOME.",
                                            name, ".rds"))
}

## function to run for single DEG res
run_reactome_singledf = function(input, name, contrast){
  
  go_res_df = data.frame()
    
    go_input <- input %>%
      #assign ENTREZID ID from symbol
      mutate(ensembl_id = mapIds(org.Hs.eg.db, keys=.$gene, column="ENTREZID", keytype="SYMBOL", multiVals="first")) %>%
      select(c(ensembl_id, avg_log2FC)) %>%
      #na.omit () %>% # remove unmapped
      arrange(desc(avg_log2FC)) %>%
      deframe()
    
    
    if(length(go_input) > 5){ ## avoid problem of small no. of  ENTREZIDs not being in list
      ## run GO enrichment
      reactenr <- gsePathway(geneList= go_input,
                             maxGSSize    = 500,
                             pvalueCutoff = 1,
                             pAdjustMethod = "BH",
                             verbose      = T)
      
      ## get result
      reactenr_tidy <-reactenr@result 
      
      if (nrow(reactenr_tidy) != 0){
        ## clean and determine GeneRatio
        gene_count <- reactenr_tidy %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
        merged_res <- left_join(reactenr_tidy, gene_count, by = "ID") %>% 
          mutate(GeneRatio = count/setSize) %>%
          mutate(zscore = scale(NES)) %>%
          mutate(group = contrast)
        
        ## add goRes to dataframe
        go_res_df = rbind(go_res_df, merged_res)
        
      }
    }

  
  ## export to global environment
  assign(x = paste0("gse_REACT.", name), value = go_res_df, envir = .GlobalEnv)
  saveRDS(object = go_res_df, file = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/geneset_enrich/","gse_REACTOME.",
                                            name, ".rds"))
}

## graph results of GSE
pway_comboDot <- function(input,
                          num_topPways = 30,
                          minimum_NES = 0.1,
                          max_padj = 0.01,
                          width=20,
                          height=10.5,
                          gs_path = "C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/geneset_enrich/plots/"){
  
  ## get name to save res to
  name_up = paste0(gsr, "_PwayDotPlot_UP")
  name_down = paste0(gsr, "_PwayDotPlot_DOWN")
  
  ### PWAYS UP # select top n pathways by NES
  pways_up <- input %>%
    dplyr::filter(p.adjust < max_padj,
                  NES > minimum_NES) %>%
    group_by(group) %>%
    arrange(desc(NES)) %>%
    slice_head(n = num_topPways) %>%
    pull(var= Description)
  
  pways_up <- unique(pways_up) ## get names of pathways to plot
  
  
  ########################################################
  if (length(pways_up) >0){
    ## cluster dataframe for better graphical representation
    mat <- input %>%
      subset(Description %in% pways_up) %>%
      select(Description, NES, group) %>%
      pivot_wider(names_from = group, values_from = NES) %>%
      column_to_rownames(var = "Description")
    
    ## scale and cluster
    mat <- scale(mat)
    dist_mat <- dist(mat, method = 'euclidean')
    dist_mat<- as.dist(dist_mat)
    ## some NAs created
    dist_mat[is.na(dist_mat)] <- 0
    dist_mat[is.nan(dist_mat)] <- 0
    sum(is.infinite(dist_mat))  # THIS SHOULD BE 0
    hclust_avg <- hclust(dist_mat, method = 'average')
    ## get plotting order
    label_order <- rownames(mat)[hclust_avg$order] 
    
    ## plot pways UP
    plot_up_df <- input %>%
      #filter(contrast %in% contrasts_use, stat == "Significant") %>%
      #mutate(contrast = factor(contrast, levels = contrasts_use)) %>%
      subset(Description %in% pways_up) %>%
      mutate(Description = factor(Description, levels =label_order ))#%>%
    #arrange(match(rownames(.), label_order))# %>%
    
    plot_up <- ggplot(plot_up_df, aes(group, Description)) +
      geom_point(aes(size=GeneRatio, fill= NES, stroke=1), shape=21) +
      scale_fill_gradient(low='#0000FF',high='#FFFF00')+
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
      scale_color_manual(values=c("red", "black"), name="Significance") +
      ggtitle(paste0(gsr, "_Pathways UP"))
    
    ## save as PDF and svg
    #up
    ggsave(filename = paste0(gs_path, name_up, ".pdf"),
           plot = plot_up,
           width=width, height=height)
    
  }
  
  
  #------------------------------------------------------------------------
  
  ### PWAYS UP # select top n pathways by NES
  pways_down <- input %>%
    dplyr::filter(p.adjust < max_padj,
                  NES < -minimum_NES) %>%
    group_by(group) %>%
    arrange(NES) %>%
    slice_head(n = num_topPways) %>%
    pull(var= Description)
  
  pways_down <- unique(pways_down) ## get names of pathways to plot
  
  
  ########################################################
  if (length(pways_down) >0){
    ## cluster dataframe for better graphical representation
    mat <- input %>%
      subset(Description %in% pways_down) %>%
      select(Description, NES, group) %>%
      pivot_wider(names_from = group, values_from = NES) %>%
      column_to_rownames(var = "Description")
    
    ## scale and cluster
    mat <- scale(mat)
    dist_mat <- dist(mat, method = 'euclidean')
    dist_mat<- as.dist(dist_mat)
    ## some NAs created
    dist_mat[is.na(dist_mat)] <- 0
    dist_mat[is.nan(dist_mat)] <- 0
    sum(is.infinite(dist_mat))  # THIS SHOULD BE 0
    hclust_avg <- hclust(dist_mat, method = 'average')
    ## get plotting order
    label_order <- rownames(mat)[hclust_avg$order] 
    
    ## plot pways UP
    plot_up_df <- input %>%
      #filter(contrast %in% contrasts_use, stat == "Significant") %>%
      #mutate(contrast = factor(contrast, levels = contrasts_use)) %>%
      subset(Description %in% pways_down) %>%
      mutate(Description = factor(Description, levels =label_order ))#%>%
    #arrange(match(rownames(.), label_order))# %>%
    
    plot_down <- ggplot(plot_up_df, aes(group, Description)) +
      geom_point(aes(size=GeneRatio, fill= NES, stroke=1), shape=21) +
      scale_fill_gradient(low='#0000FF',high='#FFFF00')+
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
      scale_color_manual(values=c("red", "black"), name="Significance") +
      ggtitle(paste0(gsr, "_Pathways DOWN"))
    
    ## save as PDF and svg
    #down
    ggsave(filename = paste0(gs_path, name_down, ".pdf"),
           plot = plot_down,
           width=width, height=height)
    
  }
  
  
}



go_res_df

pbmc.int = readRDS( file = "C:/Users/postmusd/Documents/MPXV_PBMC/hg38_10x/pbmc.2.ct_mpox.rds")

head(pbmc.int@meta.data)
unique(pbmc.int$infection)


####----------------------------------------------------------------------------
### create DEGs ####

celltypes = unique(pbmc.int@meta.data$celltype)
DefaultAssay(pbmc.int) = "RNA"

#### by celltype in MPXV vs Mock ####
head(pbmc.int@meta.data)

pbmc.int$infect_ct = paste0(pbmc.int$infection, ".", pbmc.int$celltype)
Idents(pbmc.int) = "infect_ct"

deg_Mpx.Mock_ct = list()

for (ct in celltypes){
  id1 = paste0("MPXV.", ct)
  id2 = paste0("Mock.", ct)
  
  
  degs = FindMarkers(pbmc.int, ident.1 = id1, ident.2 = id2 )
  degs$gene <- rownames(degs)
  degs$celltype <- ct
  deg_Mpx.Mock_ct[[ct]] <- degs
  
}

saveRDS(object = deg_Mpx.Mock_ct, file = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/geneset_enrich/","gse_REACTOME.",
                                          "deg_Mpx.Mock_ct", 
                                          ".rds"))


#### in 3dpi treatment ####

deg_Mpx.Mock_3dpi_ct = list()
Idents(pbmc.int) = "treat_ct"

for (ct in celltypes){
  id1 = paste0("MPXV_3dpi.", ct)
  id2 = paste0("Mock_3dpi.", ct)
  
  
  degs = FindMarkers(pbmc.int, ident.1 = id1, ident.2 = id2 )
  degs$gene <- rownames(degs)
  degs$celltype <- ct
  deg_Mpx.Mock_3dpi_ct[[ct]] <- degs
  
}

saveRDS(object = deg_Mpx.Mock_3dpi_ct, file = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/degs/",
                                                "deg_Mpx.Mock_3dpi_ct", 
                                                ".rds"))

## 5dpi 

deg_Mpx.Mock_5dpi_ct = list()
Idents(pbmc.int) = "treat_ct"

for (ct in celltypes){
  id1 = paste0("MPXV_5dpi.", ct)
  id2 = paste0("Mock_5dpi.", ct)
  
  
  degs = FindMarkers(pbmc.int, ident.1 = id1, ident.2 = id2 )
  degs$gene <- rownames(degs)
  degs$celltype <- ct
  deg_Mpx.Mock_5dpi_ct[[ct]] <- degs
  
}

saveRDS(object = deg_Mpx.Mock_5dpi_ct, file = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/degs/",
                                                     "deg_Mpx.Mock_5dpi_ct", 
                                                     ".rds"))


## DEGs Mpox positive v negative cells in select subsets
head(pbmc.int@meta.data)
pbmc.int$mpx_rna = ifelse(pbmc.int$mpx_reads_raw > 0, "pos", "neg")
pbmc.int$ct_mpx = paste0(pbmc.int$celltype,".", pbmc.int$mpx_rna)
pbmc.int$ct_mpx_tp = paste0(pbmc.int$ct_mpx,".", pbmc.int$timepoint)

table(pbmc.int$treat_ct, pbmc.int$mpx_rna)

Idents(pbmc.int) = "ct_mpx_tp"

# NK
deg_NK_PosNeg_5dpi = FindMarkers(pbmc.int, ident.1 = "NK.pos.5dpi", 
                   ident.2 = "NK.neg.5dpi")
deg_NK_PosNeg_5dpi$gene = rownames(deg_NK_PosNeg_5dpi)

# NK cycling
deg_NKcyc_PosNeg_5dpi = FindMarkers(pbmc.int, ident.1 = "NK_cycling.pos.5dpi", 
                                 ident.2 = "NK_cycling.neg.5dpi")
deg_NKcyc_PosNeg_5dpi$gene = rownames(deg_NKcyc_PosNeg_5dpi)

# Monos
deg_mono_PosNeg_5dpi = FindMarkers(pbmc.int, ident.1 = "Monocyte.pos.5dpi", 
                   ident.2 = "Monocyte.neg.5dpi")
deg_mono_PosNeg_5dpi$gene = rownames(deg_mono_PosNeg_5dpi)

# CD4 Treg
deg_cd4treg_PosNeg_5dpi = FindMarkers(pbmc.int, ident.1 = "CD4_T_Reg.pos.5dpi", 
                                   ident.2 = "CD4_T_Reg.neg.5dpi")
deg_cd4treg_PosNeg_5dpi$gene = rownames(deg_cd4treg_PosNeg_5dpi)

#### EXPORT DEGS #####

dg_path = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/degs/",
                 obj, ".xlsx") 
rm(deg_path)
rm(deg_path_rds)
dg_list = ls()[grepl(pattern = "deg",x = ls() )]
for (obj in dg_list){
  dg_path = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/degs/",
                   obj, ".xlsx") 
  dg_path_rds = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/degs/",
                       obj, ".rds") 
  writexl::write_xlsx(x = eval(parse(text = obj)),  ## eval + parse convert string to var name
                      path = dg_path)
  saveRDS(object = eval(parse(text = obj)),
          file = dg_path_rds)
}



### PathwayEnrichment DEGs #####

## GLOBAL GENE CHANGES

## run GO::BP / reactome on overall gene changes
run_go(input = deg_Mpx.Mock_ct, name = "MPXV_v_Mock_all_byct", dpi = NULL)
run_reactome(input = deg_Mpx.Mock_ct, name = "MPXV_v_Mock_all_byct", dpi = NULL)

## 3dpi
run_go(input = deg_Mpx.Mock_3dpi_ct, name = "MPXV_v_Mock_3dpi_byct", dpi =3)
run_reactome(input = deg_Mpx.Mock_3dpi_ct, name = "MPXV_v_Mock_3dpi_byct", dpi = 3)

## 5dpi
run_go(input = deg_Mpx.Mock_5dpi_ct, name = "MPXV_v_Mock_5dpi_byct", dpi = 5)
run_reactome(input = deg_Mpx.Mock_5dpi_ct, name = "MPXV_v_Mock_5dpi_byct", dpi = 5)

## POS v NEG GENE CHANGES

## Monos
run_go_singledf(input = deg_mono_PosNeg_5dpi, name = "Monos5dpi_PosVNeg", contrast = "Monos5dpi_PosVNeg")
run_reactome_singledf(input = deg_mono_PosNeg_5dpi, name = "Monos5dpi_PosVNeg", contrast = "Monos5dpi_PosVNeg")

## NK
run_go_singledf(input = deg_NK_PosNeg_5dpi, name = "NK5dpi_PosVNeg", contrast = "NK5dpi_PosVNeg")
run_reactome_singledf(input = deg_NK_PosNeg_5dpi, name = "NK5dpi_PosVNeg", contrast = "NK5dpi_PosVNeg")

## NK cycling
run_go_singledf(input = deg_NKcyc_PosNeg_5dpi, name = "NKcyc5dpi_PosVNeg", contrast = "NKcyc5dpi_PosVNeg")
run_reactome_singledf(input = deg_NKcyc_PosNeg_5dpi, name = "NKcyc5dpi_PosVNeg", contrast = "NKcyc5dpi_PosVNeg")

## CD4 Treg
run_go_singledf(input = deg_cd4treg_PosNeg_5dpi, name = "cd4treg5dpi_PosVNeg", contrast = "cd4treg5dpi_PosVNeg")
run_reactome_singledf(input = deg_cd4treg_PosNeg_5dpi, name = "cd4treg5dpi_PosVNeg", contrast = "cd4treg5dpi_PosVNeg")


#### GRAPHING RESULTS #####

dg_path = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/degs/",
                  obj, ".xlsx") 
rm(deg_path)
rm(deg_path_rds)
dg_list = ls()[grepl(pattern = "deg",x = ls() )]
for (obj in dg_list){
  dg_path = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/degs/",
                    obj, ".xlsx") 
  dg_path_rds = paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/analysis/hg38_10x/degs/",
                    obj, ".rds") 
  writexl::write_xlsx(x = eval(parse(text = obj)),  ## eval + parse convert string to var name
                      path = dg_path)
  saveRDS(object = eval(parse(text = obj)),
          file = dg_path_rds)
}

dg_list
deg_list = ls()[grepl(pattern = "gse",x = ls() )]




## GSE

head(gse_GOBP.cd4treg5dpi_PosVNeg[,8])


gs_list = ls()[grepl(pattern = "gse",x = ls() )]

for (gsr in gs_list){
  pway_comboDot(eval(parse(text = gsr)))
}









