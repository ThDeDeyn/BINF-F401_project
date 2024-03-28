library(tidyverse)
library(corrplot)

## 
# Loading data 
##

setwd("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland")

morphological_counts_lunit_dino = read_tsv("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland/Data_ThyroidGland/OG/morphological_counts_lunit_dino.tsv")

##
# Setting the parameters used throughout the code
##
p_value = .05



## 
# Performing first analysis: histograms with mean & sd
##

morpho_class <- sapply(morphological_counts_lunit_dino, class)
morpho_numeric =  morphological_counts_lunit_dino %>% select(where(is.numeric)) %>% colnames()
path_out = "Result/Data_description"


cairo_pdf(paste(path_out, "morphological_histograms.pdf", sep = "/"), onefile = T)

for(var in morpho_numeric){
  print(var)
  
  hist = hist(morphological_counts_lunit_dino[[var]],
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

cor = cor(morphological_counts_lunit_dino %>% select(morpho_numeric), use = "pairwise.complete.obs")


prefix = nchar("Mophological.cluster.")
rownames(cor) = substr(rownames(cor), prefix + 1, nchar(rownames(cor)))
colnames(cor) = rep("", ncol(cor))

Mtest = cor.mtest(morphological_counts_lunit_dino %>% select(morpho_numeric),
                  conf.level = p_value)
cor[Mtest$p > 1 - p_value] = 0


cairo_pdf(paste(path_out, "morphological_correlations.pdf", sep = "/"), onefile = T)
corrplot <- corrplot::corrplot(corr = cor, method = "color",
                              diag = FALSE, 
                              tl.cex = .75, tl.col = "black",
                              mar = c(0,0,1,0), type = 'upper',
                              col = colorRampPalette(c('#DE2F56','white','#0B65F9'))(100))

dev.off()

