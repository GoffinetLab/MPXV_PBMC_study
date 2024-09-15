library(Seurat)
library(tidyverse)
library(readxl)
#---------------------------------------


#### COUNTVIRAL READS ####

## read GTF

gtf2 = read.table("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/MPXV_NC_063383.1_annotated.gtf")


## read in viral read counts by barcode

## lib 002 and 004 # infected conditions
mpx_002 <- read.table(file = "C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/viral_reads/run1_combinedref/002_viral_posMap_BC.txt")
mpx_002$V1 <- gsub("-1", "-2", mpx_002$V1)
#mpx_002$lib <- 2
mpx_004 <- read.table(file = "C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/viral_reads/run1_combinedref/004_viral_posMap_BC.txt")
mpx_004$V1 <- gsub("-1", "-4", mpx_004$V1)
#mpx_004$lib <- 4
# combine
mpx_read <- rbind(mpx_002, mpx_004)


#############################################


mpx_read = read.table("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/viral_reads/run1_combinedref/001_viral_posMap_BC.txt")
mpx_read


#mpx_read = read.table(file = "C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/MPXV_PBMC/viral_reads/001_viral_posMap_BC.txt")
names(mpx_read) = c("bc", "umi", "tag", "cig", "start_pos")
mpx_read$end_pos = mpx_read$start_pos + 89
mpx_read$gene_start = "unmapped" ## placeholder
mpx_read$gene_end = "unmapped"
#mpx_read$letters = sub("([[:alpha:]]*).*", "\\1", mpx_read$cig)
mpx_read
## assign gene based on read position
print("Counting reads!")


for (read in seq_along(1:nrow(mpx_read))){ # for every read by row
  for (entry in seq_along(1:nrow(gtf))){ # run through every GTF entry
    if (mpx_read[read, ]$start_pos >= gtf[entry,]$V4 & mpx_read[read, ]$start_pos <= gtf[entry,]$V5){ # check if read maps to within gene
      mpx_read[read, ]$gene_start = gtf[entry,]$V16 # if maps, add gene name
    } 
    if (mpx_read[read, ]$end_pos >= gtf[entry,]$V4 & mpx_read[read, ]$end_pos <= gtf[entry,]$V5){ # check if read maps to within gene
      mpx_read[read, ]$gene_end = gtf[entry,]$V16 # if maps, add gene name
    } 
  }
}


invisible(lapply(1:5,        # Using lapply function
                 function(i) {
                   lapply(1:5)
                   print(paste("Iteration No.", i, "Created by lapply Function"))
                 }
))

############ using lapply

## define func

check_read = function(read){
  for (entry in seq_along(1:nrow(gtf))){ # run through every GTF entry
    if (mpx_read[read, ]$start_pos >= gtf[entry,]$V4 & mpx_read[read, ]$start_pos <= gtf[entry,]$V5){ # check if read maps to within gene
      mpx_read[read, ]$gene_start = gtf[entry,]$V16 # if maps, add gene name
    } 
    if (mpx_read[read, ]$end_pos >= gtf[entry,]$V4 & mpx_read[read, ]$end_pos <= gtf[entry,]$V5){ # check if read maps to within gene
      mpx_read[read, ]$gene_end = gtf[entry,]$V16 # if maps, add gene name
    } 
  }
}


lapply(seq_along(1:nrow(mpx_read)),
       
       function(read){
         for (entry in seq_along(1:nrow(gtf))){ # run through every GTF entry
           if (mpx_read[read, ]$start_pos >= gtf[entry,]$V4 & mpx_read[read, ]$start_pos <= gtf[entry,]$V5){ # check if read maps to within gene
             mpx_read[read, ]$gene_start = gtf[entry,]$V16 # if maps, add gene name
           } 
           if (mpx_read[read, ]$end_pos >= gtf[entry,]$V4 & mpx_read[read, ]$end_pos <= gtf[entry,]$V5){ # check if read maps to within gene
             mpx_read[read, ]$gene_end = gtf[entry,]$V16 # if maps, add gene name
           } 
         }
         assign("mpx_mapped", mpx_read, envir = .GlobalEnv)
       }
       
)


mpx_read

gtf2 = gtf
gtf2$pos_start = 153661
pos_start = 1533661

which(gtf2$V4 <= pos_start && gtf2$V5 >= pos_start)

list = c(153661,170425)
list = mpx_read$start_pos

getgene = function(input){
  #input = 170425
  #start = 170425
  #pos_end = input
  #pos_start = input
  #input= 139372
  gene_start = gtf2[with(gtf2, V4 <= input & V5 >= input),]$V16
  #length(gene_start)
  if(length(gene_start) >1 ){
    gene_start = "overlap"
  }
  #pos_start = as.numeric(input)
  #gene_start = gtf2[with(gtf2, V4 <= pos_start & V5 >= pos_start),]$V16
  #rm(gene_start)
  ## add read end point
  pos_end = input + 89
  #pos_end = 170425
  #pos_start
  gene_end = gtf2[with(gtf2, V4 <= pos_end & V5 >= pos_end),]$V16
  if(length(gene_end) >1 ){
    gene_start = "overlap"
  }
  #blabla = gtf2[with(gtf2, V4 <= pos_end & V5 >= pos_end),]$V16
  return(paste0(gene_start, ",", gene_end))
  #return(gene_start)
  #return(str(pos_start))
  
}

str(pos_end)

getgene(1533661)

pos_start = 10440
list = mpx_read$start_pos

mpx_read$out = map_chr(.x = list,
                       .f = getgene, .progress = T)


test = mpx_read[5:10,]
mpx_read = rbind(mpx_read, test)
######
getwd()
write.table(mpx_read, file = "mpx_read_multi.csv", 
            quote = F, sep = ";")
mpx_read = read.csv("mpx_read_multi.csv", sep = ";")
mpx_read

####

mpx_read
mpx_read2 = mpx_read %>%
  filter(!(grepl(pattern = "overlap", x = out))) %>%
  separate(out, sep = ",", into = c("gene_start", "gene_end")) %>%
  filter(!(gene_start =="" & gene_end == "")) %>% ## remove unmapped reads
  mutate(gene_assign = ifelse(gene_start == "" & gene_end != "", gene_end, ## overlap feature on right only
                              ifelse(gene_start != "" & gene_end == "", gene_start, ## overlap feature on left only
                                     ifelse(gene_start != "" & gene_end != "" & gene_start == gene_end, gene_start, "multimap")))) %>%
  filter(gene_assign != "multimap") %>%
  group_by(bc, umi, gene_assign) %>%
  tally() %>%
  group_by(bc, umi) %>%
  #filter(n_distinct(n) >1)
  ## get bc+umi+gene combo with highest counts
  filter(n == max(n)) %>%
  ## filter out any ties
  group_by(bc, umi) %>%
  #distinct(bc, umi) 
  filter(!(n() > 1))



mpx_read2 


print(mpx_read2, n= 56)[mpx_read2$bc == "AGAGAGCAGTAAACAC-1",
]
print(mpx_read2, n= 56)[mpx_read2$bc == "GATCAGTGTAATACCC-1",
]

test =print(mpx_read2, n= 56)[mpx_read2$bc == "AGAGAGCAGTAAACAC-1",
]

test2 %>%
  filter(n_distinct())


## both ends overlap feature and overlap is only within one gene

mpx_mat = mpx_read2 %>%
  #filter(gene != "unmapped") %>%
  group_by(bc, gene_assign) %>%
  tally() %>%
  pivot_wider(names_from = bc, values_from = n ) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(var = "gene_assign") %>%
  #as.matrix() %>%
  as.sparse()
mpx_mat

saveRDS(mpx_mat, file = "mpx_readcount_matrix.rds")