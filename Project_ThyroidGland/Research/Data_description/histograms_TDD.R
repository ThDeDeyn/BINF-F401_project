library(tidyverse)

## 
# Loading data 
##

setwd("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland")

clinical_data = read_tsv("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland/Data_ThyroidGland/OG/clinical_data.tsv")
morphological_counts_lunit_dino = read_tsv("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland/Data_ThyroidGland/OG/morphological_counts_lunit_dino.tsv")
RNA_read_counts = read_tsv("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland/Data_ThyroidGland/OG/RNA_read_counts.tsv")

          ## 
          # clinical_data analysis
          ##


## 
# Performing first analysis: histograms with mean & sd
##

clinical_class <- sapply(clinical_data, class)
clinical_numeric =  clinical_data %>% select(where(is.numeric)) %>% colnames()
path_out = "Result/Data_description"


cairo_pdf(paste(path_out, "clinical_histograms.pdf", sep = "/"), onefile = T)

for(var in clinical_numeric){
  print(var)
  
  hist = hist(clinical_data[[var]],
              col = "#DE2F56",        # Setting color of bars
              border = "black",       # Setting color of border
              main = "Histogram of numeric variables in clinical_data",  # Setting main title
              xlab = var,         # Setting x-axis label
              ylab = "Count")     # Setting y-axis label
  
}
dev.off()

