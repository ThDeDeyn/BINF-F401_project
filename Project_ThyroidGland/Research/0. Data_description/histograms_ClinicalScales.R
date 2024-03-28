library(tidyverse)
library(corrplot)

## 
# Loading data 
##

setwd("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland")

clinical_data = read_tsv("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland/Data_ThyroidGland/OG/clinical_data.tsv")

##
# Setting the parameters used throughout the code
##
p_value = .05


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


## 
# Performing second analysis: correlations within clinical_data  
##

cor = cor(clinical_data %>% select(clinical_numeric), use = "pairwise.complete.obs")

#colnames(cor) = rep("", ncol(cor))

Mtest = cor.mtest(clinical_data %>% select(clinical_numeric),
                  conf.level = p_value)
cor[Mtest$p > 1 - p_value] = 0


cairo_pdf(paste(path_out, "clinical_correlations.pdf", sep = "/"), onefile = T)
corrplot <- corrplot::corrplot(corr = cor, method = "number",
                              diag = FALSE, 
                              tl.cex = .75, tl.col = "black",
                              mar = c(0,0,1,0), type = 'upper',
                              col = colorRampPalette(c('#DE2F56','white','#0B65F9'))(100))

dev.off()

