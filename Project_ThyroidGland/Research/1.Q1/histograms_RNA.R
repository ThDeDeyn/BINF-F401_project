library(tidyverse)
library(corrplot)

## 
# Loading data 
##

setwd("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland")

RNA_read_counts = read_tsv("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland/Data_ThyroidGland/OG/RNA_read_counts.tsv")

##
# Setting the parameters used throughout the code
##
p_value = .05
          ## 
          # clinical_data analysis
          ##


## 
# Performing first analysis: histograms with mean & sd
##

RNA_class <- sapply(RNA_read_counts, class)
RNA_numeric =  RNA_read_counts %>% select(where(is.numeric)) %>% colnames()
path_out = "Result/Data_description"


cairo_pdf(paste(path_out, "RNA_histograms.pdf", sep = "/"), onefile = T)

for(var in RNA_numeric){
  print(var)
  
  hist = hist(RNA_numeric[[var]],
              col = "#DE2F56",        # Setting color of bars
              border = "black",       # Setting color of border
              main = "Histogram of numeric variables in clinical_data",  # Setting main title
              xlab = var,         # Setting x-axis label
              ylab = "Count")     # Setting y-axis label
  
}
dev.off()


## 
# Performing second analysis: correlations within clinical_data  
##
