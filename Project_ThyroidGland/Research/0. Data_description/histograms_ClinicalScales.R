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

## 
# Performing third analysis: correlations within clinical_data  
##

highest_cor = cor
# Keep high correlations
highest_cor[abs(cor) < 0.2] = 0

# Keep significant correlations
highest_cor[(Mtest$p) > 0.05] = 0

# Remove self correlations
highest_cor[abs(cor) == 1] = 0

cor_list = data.frame(var.x = character(0), var.y  = character(0),
                      cor = numeric(0), p = numeric(0))
for(var.x in colnames(cor)){
  for(var.y in colnames(cor)){
    if(var.x != var.y){
      new_row = data.frame(var.x = var.x,
                           var.y = var.y,
                           cor = cor(clinical_data[[var.x]], clinical_data[[var.y]],
                                     use = "pairwise.complete.obs"), 
                           p = cor.test(clinical_data[[var.x]], clinical_data[[var.y]],
                                     use = "pairwise.complete.obs")$p.value)
      
      if(new_row$p < 0.05 & abs(new_row$cor) > 0.15){
        cor_list = rbind(cor_list, new_row)
      }                  
    }
  }
}

cairo_pdf(paste(path_out, "clinical_point.pdf", sep = "/"), onefile = T)

for(comp in 1:nrow(cor_list)){
  var.x = cor_list[comp, "var.x"]
  var.y = cor_list[comp, "var.y"]
  
 clinical_data %>% ggplot() + geom_point(aes( x = clinical_data[[var.x]], y = clinical_data[[var.y]]))
}

dev.off()

