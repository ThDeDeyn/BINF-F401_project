library(tidyverse)

setwd("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland")
##
# Merging in a single dataframe
## 

### Reading the DF

RNA_read_counts = read_tsv("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland/Data_ThyroidGland/OG/RNA_read_counts.tsv")
morphological_counts_lunit_dino = read_tsv("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland/Data_ThyroidGland/OG/morphological_counts_lunit_dino.tsv")
clinical_data = read_tsv("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland/Data_ThyroidGland/OG/clinical_data.tsv")


### Pre-processing df
    # Invert 
RNA_read_counts = t(RNA_read_counts)
RNA_read_counts = as.data.frame(RNA_read_counts)

    # Fine tune columns and rownames


dico = RNA_read_counts[1:2,]
colnames(RNA_read_counts) = RNA_read_counts[1,]

RNA_read_counts = RNA_read_counts[-c(1:2),]
RNA_read_counts = RNA_read_counts %>% mutate(SMPLID = rownames(RNA_read_counts))
  
### Merging df

study_dataframe = clinical_data
study_dataframe = study_dataframe %>% left_join(y = morphological_counts_lunit_dino, by = "SMPLID")
study_dataframe = study_dataframe %>% left_join(y = RNA_read_counts, by = "SMPLID")


### Save new df

path_out = "Data_ThyroidGland/OG"
saveRDS(file = paste0(path_out, "/study_dataframe.rds"),object = study_dataframe)

# tmp = readRDS(file = paste0(path_out, "/study_dataframe.rds"))
