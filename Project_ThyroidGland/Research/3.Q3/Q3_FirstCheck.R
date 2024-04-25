library(tidyverse)

setwd("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland")

RNA_read = read_tsv("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland/Data_ThyroidGland/OG/RNA_read_counts.tsv")
morphological_counts = read_tsv("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland/Data_ThyroidGland/OG/morphological_counts_lunit_dino.tsv")

clinical_data = read_tsv("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland/Data_ThyroidGland/OG/clinical_data.tsv")


rownames(RNA_read) = RNA_read$Name 
df = RNA_read %>% select_if(is.numeric)
rownames(df) = RNA_read$Name 

med_avg = apply(df,1, function(var) mad(var, na.rm = TRUE))

hist(med_avg)
hist(log(med_avg))

max(med_avg)

# Suspectly high value --> too 
pos.max = which.max(med_avg)


#View(t(df[pos.max,])[,1])
hist(t(df[pos.max ,])[,1])
hist(log10(t(df[pos.max ,])[,1]))

