# Bacteriome (Pathobiont) analysis 
library(dplyr)
library(tidyverse)
library(Maaslin2)
library(ComplexHeatmap)

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Pathobiont")
#UCvsControl
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% filter(Meta$UC_Dx == 1| Meta$UC_Dx ==0) -> Meta

SP_input <- read.csv("■Pathogen_ncgm.4198.pachogen.rpkm.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(Meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("IBD_HC2", "Age", "sex", "BMI", "IBD_HC2.1", "UC", "Crohn", "UC_Dx", "CD_Dx", "IBD"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2
SP_input2 %>% select("Enterococcus.faecium", "Staphylococcus.aureus",  "Klebsiella.pneumoniae", "Acinetobacter.baumannii", "Pseudomonas.aeruginosa") -> SP_input2 #ESKAPE are selected.

rownames(Meta) <- NULL　
column_to_rownames(Meta, "ID") -> Meta
Meta %>% select(UC_Dx, Age, sex, BMI) -> JP_meta

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

SP_input <- read.csv("■Pathogen_ncgm.4198.pachogen.rpkm.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(Meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("IBD_HC2", "Age", "sex", "BMI", "IBD_HC2.1", "UC", "Crohn", "UC_Dx", "CD_Dx", "IBD"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2
SP_input2 %>% select("Enterococcus.faecium", "Staphylococcus.aureus",  "Klebsiella.pneumoniae", "Acinetobacter.baumannii", "Pseudomonas.aeruginosa") -> SP_input2 #ESKAPE are selected.

rownames(Meta) <- NULL
column_to_rownames(Meta, "ID") -> Meta
Meta %>% select(CD_Dx, Age, sex, BMI) -> JP_meta

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
Meta %>% filter(Meta$IBD == 1| Meta$IBD ==0) -> Meta

SP_input <- read.csv("■Pathogen_ncgm.4198.pachogen.rpkm.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(Meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("IBD_HC2", "Age", "sex", "BMI", "IBD_HC2.1", "UC", "Crohn", "UC_Dx", "CD_Dx", "IBD"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2
SP_input2 %>% select("Enterococcus.faecium", "Staphylococcus.aureus",  "Klebsiella.pneumoniae", "Acinetobacter.baumannii", "Pseudomonas.aeruginosa") -> SP_input2 #ESKAPE are selected.

rownames(Meta) <- NULL　
column_to_rownames(Meta, "ID") -> Meta
Meta %>% select(IBD, Age, sex, BMI) -> JP_meta

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

Profile <- read.csv("Profile.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
rownames_to_column(Profile, "feature") -> Profile
Patho_input <- full_join(full_join(JP_maaslin_IBD, JP_maaslin_UC, by = "feature"), JP_maaslin_CD, by = "feature")
Patho_input <- left_join(Patho_input, Profile, by = "feature")

Patho_input$feature <- str_replace_all(Patho_input $feature, pattern="[.]", replacement=" ") 
column_to_rownames(Patho_input, "feature") -> Patho_input
Patho_input %>% arrange(-Patho_input$"Coefficient (CD JP)") -> Patho_input

Patho_input %>% select("Q-value (IBD JP)") -> qval_IBD
Patho_input %>% select("Q-value (UC JP)", "Q-value (CD JP)") -> qval_UC_CD
qval_IBD[is.na(qval_IBD)] <- 1
qval_UC_CD[is.na(qval_UC_CD)] <- 1

Patho_input  %>% select("Coefficient (IBD JP)") -> Patho_coef_IBD
Patho_input  %>% select("Coefficient (UC JP)" , "Coefficient (CD JP)") -> Patho_coef_UC_CD
Patho_input  %>% select(Gram) -> Gram
Patho_input  %>% select(Phylum) -> Phylum

rename(.data= Patho_coef_IBD, "IBD" = "Coefficient (IBD JP)") -> Patho_coef_IBD
rename(.data= Patho_coef_UC_CD, "UC" = "Coefficient (UC JP)") -> Patho_coef_UC_CD 
rename(.data= Patho_coef_UC_CD, "CD" = "Coefficient (CD JP)") -> Patho_coef_UC_CD

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1 (MaAsLin)")

ra = rowAnnotation("Phylum"=Patho_input$Phylum, "Gram staining"=Patho_input$Gram,  col=list(Phylum= c("Actinobacteria"="#03af7a", "Bacteroidetes" = "#c9ace6", "Firmicutes" = "#ff8082", "Proteobacteria" = "#804000", "Bacteria phylum incertae sedis" = "#ffff80", "Verrucomicrobia"="blue", "Fusobacteria"="#bfe4ff", "Synergistetes"="#ffca80", "Spirochaetes"="red", "Unknown"="#000000"), "Gram staining"= c("positive"="purple", "negative"="#ffcabf", "unknown"="#c8c8cb")), simple_anno_size = unit(3, "mm"), gp = gpar(col = "black"), annotation_name_gp= gpar(fontsize = 8))

p1=pheatmap(as.matrix(Patho_coef_IBD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(qval_IBD<0.1,"*", ""), nrow(qval_IBD)), fontsize_number = 7, border_color = TRUE, col = circlize::colorRamp2(c(-1, 0, 4), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))
p2=pheatmap(as.matrix(Patho_coef_UC_CD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(qval_UC_CD<0.1,"*", ""), nrow(qval_UC_CD)), fontsize_number = 7, border_color = TRUE, col = circlize::colorRamp2(c(-1, 0, 4), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))

draw(ra+p1+p2,  heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig)) 
