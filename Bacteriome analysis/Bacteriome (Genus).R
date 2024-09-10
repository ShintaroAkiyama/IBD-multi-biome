# Bacteriome (Genus) analysis (Supplementary Figure 9)
library(dplyr)
library(tidyverse)
library(Maaslin2)
library(ComplexHeatmap)

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Genus")
#UCvsControl
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% filter(Meta$UC_Dx == 1| Meta$UC_Dx ==0) -> Meta
Meta %>% select(UC_Dx, Age, sex, BMI) -> JP_meta

SP_input <- read.csv("Genus.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
SP_input %>% filter(SP_input$IBD ==2 |SP_input$UC ==1) -> SP_input
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("ResID", "IBD", "UC", "CD"))] 

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
  fixed_effects=c("UC_Dx,Age,sex,BMI"),
  reference = c("UC_Dx,0")) 

fit_data$results -> JP_maaslin
JP_maaslin %>% filter(JP_maaslin$"metadata" == "UC_Dx") -> JP_maaslin_UC
JP_maaslin_UC <- JP_maaslin_UC[, -which (colnames(JP_maaslin_UC) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=JP_maaslin_UC, "Coefficient (UC JP)" = "coef") -> JP_maaslin_UC
rename(.data=JP_maaslin_UC, "P-value (UC JP)" = "pval") -> JP_maaslin_UC
rename(.data=JP_maaslin_UC, "Q-value (UC JP)" = "qval") -> JP_maaslin_UC

#CDvsControl
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% filter(Meta$CD_Dx == 1| Meta$CD_Dx ==0) -> Meta
Meta %>% select(CD_Dx, Age, sex, BMI) -> JP_meta

SP_input <- read.csv("Genus.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) # abundance data (normalized)
SP_input %>% filter(SP_input$IBD ==2 |SP_input$CD ==1) -> SP_input
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("ResID", "IBD", "UC", "CD"))] 

keep <- apply(SP_input2, 2, mean) > 1E-4 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
SP_input2  <- SP_input2[, keep]

fit_data = Maaslin2(
  input_data = SP_input2, 
  input_metadata = JP_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "JP_CD", 
  normalization = "NONE",
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  standardize = FALSE, 
  cores = 4,
  fixed_effects=c("CD_Dx,Age,sex,BMI"),
  reference = c("CD_Dx,0")) 

fit_data$results -> JP_maaslin
JP_maaslin %>% filter(JP_maaslin$"metadata" == "CD_Dx") -> JP_maaslin_CD
JP_maaslin_CD <- JP_maaslin_CD[, -which (colnames(JP_maaslin_CD) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=JP_maaslin_CD, "Coefficient (CD JP)" = "coef") -> JP_maaslin_CD
rename(.data=JP_maaslin_CD, "P-value (CD JP)" = "pval") -> JP_maaslin_CD
rename(.data=JP_maaslin_CD, "Q-value (CD JP)" = "qval") -> JP_maaslin_CD

#IBDvsControl
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% select(IBD, Age, sex, BMI) -> JP_meta

SP_input <- read.csv("Genus.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) # abundance data (normalized)
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("ResID", "IBD", "UC", "CD"))] 

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
  fixed_effects=c("IBD,Age,sex,BMI"),
  reference = c("IBD,0")) 

fit_data$results -> JP_maaslin
JP_maaslin %>% filter(JP_maaslin$"metadata" == "IBD") -> JP_maaslin_IBD
JP_maaslin_IBD <- JP_maaslin_IBD[, -which (colnames(JP_maaslin_IBD) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=JP_maaslin_IBD, "Coefficient (IBD JP)" = "coef") -> JP_maaslin_IBD
rename(.data=JP_maaslin_IBD, "P-value (IBD JP)" = "pval") -> JP_maaslin_IBD
rename(.data=JP_maaslin_IBD, "Q-value (IBD JP)" = "qval") -> JP_maaslin_IBD

#Data integration
JP_MA <- full_join(full_join(JP_maaslin_IBD, JP_maaslin_UC, by = "feature"), JP_maaslin_CD, by = "feature")

#Heatmap creation
#Select species with abs coef >1 and FDR<0.1 in UC or CD 
JP_MA %>% filter(abs(JP_MA$"Coefficient (CD JP)") > 1 & JP_MA$"Q-value (CD JP)" < 0.1 | abs(JP_MA$"Coefficient (UC JP)") > 1 & JP_MA$"Q-value (UC JP)" < 0.1) -> SP_input
SP_input %>% arrange(-SP_input$"Coefficient (IBD JP)") -> SP_input2 
column_to_rownames(SP_input2, "feature") -> SP_input2

SP_input2 %>% select("Coefficient (IBD JP)", "P-value (IBD JP)", "Q-value (IBD JP)", "Coefficient (UC JP)", "P-value (UC JP)", "Q-value (UC JP)", "Coefficient (CD JP)", "P-value (CD JP)", "Q-value (CD JP)") ->  Fig1 

Fig1 %>% select("Q-value (IBD JP)") -> qval_IBD
Fig1 %>% select("Q-value (UC JP)", "Q-value (CD JP)") -> qval_UC_CD

qval_IBD[is.na(qval_IBD)] <- 1
qval_UC_CD[is.na(qval_UC_CD)] <- 1

Fig1 %>% select("Coefficient (IBD JP)") -> Gram_coef_IBD
Fig1 %>% select("Coefficient (UC JP)", "Coefficient (CD JP)") -> Gram_coef_UC_CD

rename(.data= Gram_coef_IBD, "IBD" = "Coefficient (IBD JP)") -> Gram_coef_IBD
rename(.data= Gram_coef_UC_CD, "UC" = "Coefficient (UC JP)") -> Gram_coef_UC_CD 
rename(.data= Gram_coef_UC_CD, "CD" = "Coefficient (CD JP)") -> Gram_coef_UC_CD

anno_width = unit(2, "cm")

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")

p1=pheatmap(as.matrix(Gram_coef_IBD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(qval_IBD < 0.1,"*", ""), nrow(qval_IBD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-4,0,3), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))
p2=pheatmap(as.matrix(Gram_coef_UC_CD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(qval_UC_CD < 0.1,"*", ""), nrow(qval_UC_CD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-4,0,3), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))

draw(p1+p2, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig))
