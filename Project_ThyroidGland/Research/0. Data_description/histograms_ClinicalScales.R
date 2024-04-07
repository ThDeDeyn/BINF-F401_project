library(tidyverse)
library(corrplot)

## 
# Loading data 
##

setwd("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland")

clinical_data = read_tsv("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland/Data_ThyroidGland/OG/clinical_data.tsv")


#clinical_data$SEX = as.factor(clinical_data$SEX)
clinical_data$DTHHRDY = as.factor(clinical_data$DTHHRDY)
clinical_data$COHORT = as.factor(clinical_data$COHORT)

clinical_data$DTHVNT[clinical_data$DTHVNT == 99] = NA

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
cor[Mtest$p > p_value] = 0


cairo_pdf(paste(path_out, "clinical_correlations.pdf", sep = "/"), onefile = T)
corrplot <- corrplot::corrplot(corr = cor, method = "number",
                              diag = FALSE, 
                              tl.cex = .75, tl.col = "black",
                              mar = c(0,0,1,0), type = 'upper',
                              col = colorRampPalette(c('#DE2F56','white','#0B65F9'))(100))

dev.off()

## 
# Performing third analysis: correlations within clinical_data (numeric VS numeric)  
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

var.y_list = colnames(cor)

for(var.x in colnames(cor)){
  for(var.y in var.y_list){
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
  var.y_list = var.y_list[-which(var.y_list == var.x)]
}


cairo_pdf(paste(path_out, "clinical_point.pdf", sep = "/"), onefile = T)

for(comp in 1:nrow(cor_list)){
  var.x = cor_list[comp, "var.x"]
  var.y = cor_list[comp, "var.y"]
  
 plot = clinical_data %>% ggplot() + geom_point(aes( x = clinical_data[[var.x]], y = clinical_data[[var.y]], 
                                                    alpha = 0.5), 
                                                #position = position_jitter(w = 0, h = 0),  
                                                col = '#DE2F56') + 
                                      xlab(var.x) + ylab(var.y) + 
                                      theme_minimal()
 print(plot)
}

dev.off()

## 
# Performing fourth analysis: boxplots of values separated using most correlated (numeric VS factors)
##

clinical_numeric = clinical_numeric
clinical_factor =  clinical_data %>% select(where(is.factor)) %>% colnames()

    ## 1: COHORT vs numeric

cairo_pdf(paste(path_out, "clinical_boxplots.pdf", sep = "/"), onefile = T)

comp_list = data.frame(var.x = character(0), var.y  = character(0),
                      mean.x = numeric(0), mean.y = numeric(0),
                      p = numeric(0))

for(var.x in clinical_numeric){
  tmp.A = clinical_data %>% filter(clinical_data$COHORT == "Postmortem") %>% select(var.x)
  tmp.B = clinical_data %>% filter(clinical_data$COHORT != "Postmortem") %>% select(var.x)
  
  t.test = t.test(tmp.A, tmp.B)
  
  new_row = data.frame(var.x = var.x,
                       var.y = "COHORT",
                       mean.x = t.test$estimate[1],
                       mean.y = t.test$estimate[2], 
                       p = t.test$p.value)
  
  if(new_row$p < 0.05){
    comp_list = rbind(comp_list, new_row)
    plot = clinical_data %>% ggplot(aes(x = clinical_data[[var.x]], y = clinical_data[["COHORT"]], fill = clinical_data[["COHORT"]])) +
      geom_boxplot(outlier.color = "black") + 
      xlab(var.x) +
      theme_minimal()
    print(plot)
    }  
}

    ## 2: DTHHRDY vs numeric

comp_list = data.frame(var.x = character(0), var.y  = character(0),
                       mean.x = numeric(0), mean.y = numeric(0),
                       p = numeric(0))

for(var.x in clinical_numeric){

  eq = reformulate(var.x, response = "DTHHRDY")
  #eq = paste0(var.x, sep = " ~ ", "DTHHRDY")
  aov = aov(eq, data = clinical_data)
  aov$coefficients
  new_row = data.frame(var.x = var.x,
                       var.y = "COHORT",
                       mean.x = t.test$estimate[1],
                       mean.y = t.test$estimate[2], 
                       p = t.test$p.value)
  
    comp_list = rbind(comp_list, new_row)
    plot = clinical_data %>% ggplot(aes(x = clinical_data[[var.x]], y = clinical_data[["DTHHRDY"]], fill = clinical_data[["DTHHRDY"]])) +
      geom_boxplot(outlier.color = "black") + 
      xlab(var.x) +
      theme_minimal()
    print(plot)
}
dev.off()

## 
# Performing fifth analysis: PCA
##

library(caret)
library(FactoMineR)


# A: scaling numeric df

cairo_pdf(paste(path_out, "clinical_PCA_1.pdf", sep = "/"), onefile = T)

pca_df = clinical_data %>% dplyr::select(clinical_numeric)
pca_df = as.data.frame(scale(pca_df) )

PCA = PCA(pca_df, graph = TRUE)

eigenvalues = PCA$eig
barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], 
      type="b", pch=19, col = "red")

cumul.eigenvalues = cumsum(eigenvalues[, 2])
barplot(cumul.eigenvalues, names.arg=1:nrow(eigenvalues), 
        main = "Cumulative Explained Variance",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
lines(x = 1:nrow(eigenvalues), cumul.eigenvalues, 
      type="b", pch=19, col = "red")

dev.off()

## 
# Performing fifth analysis: MFA (Multiple Factors Analysis)
##
library("factoextra")

cairo_pdf(paste(path_out, "clinical_MFA_1.pdf", sep = "/"), onefile = T)

#lapply(mfa_df, function(var) class(var))
mfa_df = clinical_data %>% select(-c("IMGURL", "SUBJID","SMPTHNTS"), -matches("SMPLID"))
#mfa_df$SEX = as.factor(mfa_df$SEX)

#mfa_df %>% colnames()
mfa_df <- mfa_df[, c(  "SEX",     "AGE",     "HGHT",    "WGHT",    "BMI",     "TRISCHD", "DTHVNT",  "DTHHRDY", "COHORT")]

group.name = c("patients", "death_numeric", "death_factor" )
MFA = MFA(mfa_df, 
          group = c(5, 2,2), 
          type = c("c","c","n"),
          name.group = group.name)


# Study the eigenvalues
eigenvalues = MFA$eig
barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], 
      type="b", pch=19, col = "red")

cumul.eigenvalues = cumsum(eigenvalues[, 2])
barplot(cumul.eigenvalues, names.arg=1:nrow(eigenvalues), 
        main = "Cumulative Explained Variance",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
lines(x = 1:nrow(eigenvalues), cumul.eigenvalues, 
      type="b", pch=19, col = "red")

#Study the contributions
fviz_mfa_var(MFA, "group")

dev.off()
