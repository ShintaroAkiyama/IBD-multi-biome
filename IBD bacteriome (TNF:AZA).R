# Analysis of bacteriome in patients on TNF inhibitors
library(dplyr)
library(tidyverse)
library(Maaslin2)
library(ComplexHeatmap)

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/SP_mOTU3(METAFは番号に修正済みIBD症状薬剤解析)")
#IBD
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% mutate("TNF" = if_else(Meta$"Infliximab.or.adalimumab" == 1 |Meta$"Immunomodulators" == 1, 1, 0)) -> Meta
Meta %>% select(Age, Sex, BMI, "TNF") -> JP_meta

rownames_to_column(JP_meta, "ID") -> JP_meta
JP_meta %>% mutate(ID = as.numeric(JP_meta$ID)) ->JP_meta

SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalizaed)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input

left_join(JP_meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("Age", "Sex", "BMI", "TNF", "METAF", "IBD", "UC", "CD"))] 

column_to_rownames(SP_input2, "ID") -> SP_input2
column_to_rownames(JP_meta, "ID") -> JP_meta

keep <- apply(SP_input2, 2, mean) > 1E-4 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
SP_input2  <- SP_input2[, keep]

fit_data = Maaslin2(
  input_data = SP_input2, 
  input_metadata = JP_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "JP_IBD", 
  normalization = "NONE",
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  standardize = FALSE, 
  cores = 4,
  fixed_effects=c("TNF,Age,Sex,BMI"),
  reference = c("TNF,0")) 

fit_data$results -> JP_maaslin
JP_maaslin %>% filter(JP_maaslin$"metadata" == "TNF") -> JP_maaslin_IBD
JP_maaslin_IBD <- JP_maaslin_IBD[, -which (colnames(JP_maaslin_IBD) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=JP_maaslin_IBD, "Coefficient (IBD TNF/AZA)" = "coef") -> JP_maaslin_IBD
rename(.data=JP_maaslin_IBD, "P-value (IBD TNF/AZA)" = "pval") -> JP_maaslin_IBD
rename(.data=JP_maaslin_IBD, "Q-value (IBD TNF/AZA)" = "qval") -> JP_maaslin_IBD

#UC
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% mutate("TNF" = if_else(Meta$"Infliximab.or.adalimumab" == 1 |Meta$"Immunomodulators" == 1, 1, 0)) -> Meta
Meta %>% select(Age, Sex, BMI, UC, Crohn, TNF) -> JP_meta

rownames_to_column(JP_meta, "ID") -> JP_meta
JP_meta %>% mutate(ID = as.numeric(JP_meta$ID)) ->JP_meta
JP_meta %>% filter(JP_meta$UC == 1) -> JP_meta #Select UC

SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalizaed)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input

left_join(JP_meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("Age", "Sex", "BMI", "UC.x", "Crohn", "TNF", "METAF", "IBD", "UC.y", "CD"))] 

column_to_rownames(SP_input2, "ID") -> SP_input2
column_to_rownames(JP_meta, "ID") -> JP_meta

keep <- apply(SP_input2, 2, mean) > 1E-4 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
SP_input2  <- SP_input2[, keep]

fit_data = Maaslin2(
  input_data = SP_input2, 
  input_metadata = JP_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "JP_UC", 
  normalization = "NONE",
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  standardize = FALSE, 
  cores = 4,
  fixed_effects=c("TNF,Age,Sex,BMI"),
  reference = c("TNF,0")) 

fit_data$results -> JP_maaslin
JP_maaslin %>% filter(JP_maaslin$"metadata" == "TNF") -> JP_maaslin_UC
JP_maaslin_UC <- JP_maaslin_UC[, -which (colnames(JP_maaslin_UC) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=JP_maaslin_UC, "Coefficient (UC TNF/AZA)" = "coef") -> JP_maaslin_UC
rename(.data=JP_maaslin_UC, "P-value (UC TNF/AZA)" = "pval") -> JP_maaslin_UC
rename(.data=JP_maaslin_UC, "Q-value (UC TNF/AZA)" = "qval") -> JP_maaslin_UC

#CD
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% mutate("TNF" = if_else(Meta$"Infliximab.or.adalimumab" == 1 |Meta$"Immunomodulators" == 1, 1, 0)) -> Meta
Meta %>% select(Age, Sex, BMI, UC, Crohn, TNF) -> JP_meta

rownames_to_column(JP_meta, "ID") -> JP_meta
JP_meta %>% mutate(ID = as.numeric(JP_meta$ID)) ->JP_meta
JP_meta %>% filter(JP_meta$Crohn == 1) -> JP_meta #Select CD
SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalizaed)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input

left_join(JP_meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("Age", "Sex", "BMI", "UC.x", "Crohn", "TNF", "METAF", "IBD", "UC.y", "CD"))] 

column_to_rownames(SP_input2, "ID") -> SP_input2
column_to_rownames(JP_meta, "ID") -> JP_meta

keep <- apply(SP_input2, 2, mean) > 1E-4 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
SP_input2  <- SP_input2[, keep]

fit_data = Maaslin2(
  input_data = SP_input2, 
  input_metadata = JP_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "JP_UC", 
  normalization = "NONE",
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  standardize = FALSE, 
  cores = 4,
  fixed_effects=c("TNF,Age,Sex,BMI"),
  reference = c("TNF,0")) 

fit_data$results -> JP_maaslin
JP_maaslin %>% filter(JP_maaslin$"metadata" == "TNF") -> JP_maaslin_CD
JP_maaslin_CD <- JP_maaslin_CD[, -which (colnames(JP_maaslin_CD) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=JP_maaslin_CD, "Coefficient (CD TNF/AZA)" = "coef") -> JP_maaslin_CD
rename(.data=JP_maaslin_CD, "P-value (CD TNF/AZA)" = "pval") -> JP_maaslin_CD
rename(.data=JP_maaslin_CD, "Q-value (CD TNF/AZA)" = "qval") -> JP_maaslin_CD

#Data integration
JP_MA <- full_join(full_join(JP_maaslin_IBD, JP_maaslin_UC, by = "feature"), JP_maaslin_CD, by = "feature")
JP_MA %>% mutate(ID = str_sub(JP_MA$feature, start = -6, end = -2)) -> JP_MA

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/mOTU3_original_data")
ID <- read.csv("Feature.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
JP_MA %>% mutate(ID = as.numeric(JP_MA$ID)) -> JP_MA
JP_MA <- full_join(JP_MA, ID, by ="ID")

#Heatmap creation
#Select speceis with abs coef >0.8 in UC or CD
JP_MA %>% filter(abs(JP_MA$"Coefficient (CD TNF/AZA)") > 0.8 | abs(JP_MA$"Coefficient (UC TNF/AZA)") > 0.8) -> SP_input
SP_input %>% arrange(-SP_input$"Coefficient (IBD TNF/AZA)") -> SP_input2 
column_to_rownames(SP_input2, "Feature") -> SP_input2

SP_input2 %>% select("Coefficient (IBD TNF/AZA)", "P-value (IBD TNF/AZA)", "Q-value (IBD TNF/AZA)", 
                     "Coefficient (UC TNF/AZA)", "P-value (UC TNF/AZA)", "Q-value (UC TNF/AZA)", 
                     "Coefficient (CD TNF/AZA)", "P-value (CD TNF/AZA)", "Q-value (CD TNF/AZA)") ->  Fig1 

Fig1 %>% select("Q-value (IBD TNF/AZA)") -> qval_IBD
Fig1 %>% select("Q-value (UC TNF/AZA)", "Q-value (CD TNF/AZA)") -> qval_UC_CD
qval_IBD[is.na(qval_IBD)] <- 1
qval_UC_CD[is.na(qval_UC_CD)] <- 1

Fig1 %>% select("Coefficient (IBD TNF/AZA)") -> Gram_coef_IBD
Fig1 %>% select("Coefficient (UC TNF/AZA)", "Coefficient (CD TNF/AZA)") -> Gram_coef_UC_CD

Fig1 %>% select("P-value (IBD TNF/AZA)") -> pval_IBD
Fig1 %>% select("P-value (UC TNF/AZA)", "P-value (CD TNF/AZA)") -> pval_UC_CD
pval_IBD[is.na(pval_IBD)] <- 1
pval_UC_CD[is.na(pval_UC_CD)] <- 1

anno_width = unit(2, "cm")

lgd_sig = Legend(pch = "*", type = "points", labels = "P < 0.05")

p1=pheatmap(as.matrix(Gram_coef_IBD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(pval_IBD < 0.05,"*", ""), nrow(pval_IBD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-2.5,0,1.5), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))
p2=pheatmap(as.matrix(Gram_coef_UC_CD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(pval_UC_CD < 0.05,"*", ""), nrow(pval_UC_CD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-2.5,0,1.5), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))

draw(p1+p2, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig))