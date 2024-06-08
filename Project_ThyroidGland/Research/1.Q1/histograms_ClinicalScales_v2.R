library(tidyverse)
library(corrplot)

## 
# Loading data 
##

setwd("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland")
clinical_data = read_tsv("D:/Polytech/MA2 2023-2024/Q2/Genomics/Project_ThyroidGland/Data_ThyroidGland/OG/clinical_data.tsv")


clinical_data$SEX = as.factor(clinical_data$SEX)
clinical_data$DTHHRDY = as.factor(clinical_data$DTHHRDY)
clinical_data$COHORT = as.factor(clinical_data$COHORT)
clinical_data$DTHVNT = as.factor(clinical_data$DTHVNT)


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
      ylab("COHORT") +
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
                       var.y = "DTHHRDY",
                       mean.x = t.test$estimate[1],
                       mean.y = t.test$estimate[2], 
                       p = t.test$p.value)
  
  comp_list = rbind(comp_list, new_row)
  plot = clinical_data %>% ggplot(aes(x = clinical_data[[var.x]], y = clinical_data[["DTHHRDY"]], fill = clinical_data[["DTHHRDY"]])) +
    geom_boxplot(outlier.color = "black") + 
    xlab(var.x) +
    ylab("DTHHRDY") + 
    theme_minimal()
  print(plot)
}

## 3: SEX vs numeric

comp_list = data.frame(var.x = character(0), var.y  = character(0),
                       mean.x = numeric(0), mean.y = numeric(0),
                       p = numeric(0))

for(var.x in clinical_numeric){
  tmp.A = clinical_data %>% filter(clinical_data$SEX == "1") %>% select(var.x)
  tmp.B = clinical_data %>% filter(clinical_data$SEX != "1") %>% select(var.x)
  
  t.test = t.test(tmp.A, tmp.B)
  
  new_row = data.frame(var.x = var.x,
                       var.y = "SEX",
                       mean.x = t.test$estimate[1],
                       mean.y = t.test$estimate[2], 
                       p = t.test$p.value)
  
  if(new_row$p < 0.05){
    print(paste0("Sex and ", var.x))
    comp_list = rbind(comp_list, new_row)
    plot = clinical_data %>% ggplot(aes(x = clinical_data[[var.x]], y = clinical_data[["SEX"]], fill = clinical_data[["SEX"]])) +
      geom_boxplot(outlier.color = "black") + 
      xlab(var.x) +
      ylab("SEX") +
      theme_minimal()
    print(plot)
  }  
}

## 4: DTHVNT vs numeric

comp_list = data.frame(var.x = character(0), var.y  = character(0),
                       mean.x = numeric(0), mean.y = numeric(0),
                       p = numeric(0))

for(var.x in clinical_numeric){
  tmp.A = clinical_data %>% filter(clinical_data$DTHVNT == "0") %>% select(var.x)
  tmp.B = clinical_data %>% filter(clinical_data$DTHVNT != "0") %>% select(var.x)
  
  t.test = t.test(tmp.A, tmp.B)
  
  new_row = data.frame(var.x = var.x,
                       var.y = "DTHVNT",
                       mean.x = t.test$estimate[1],
                       mean.y = t.test$estimate[2], 
                       p = t.test$p.value)
  
  if(new_row$p < 0.05){
    print(paste0("DTHVNT and ", var.x))
    comp_list = rbind(comp_list, new_row)
    plot = clinical_data %>% ggplot(aes(x = clinical_data[[var.x]], y = clinical_data[["DTHVNT"]], fill = clinical_data[["DTHVNT"]])) +
      geom_boxplot(outlier.color = "black") + 
      xlab(var.x) +
      ylab("DTHVNT") +
      theme_minimal()
    print(plot)
  }  
}

dev.off()
