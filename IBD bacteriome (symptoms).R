# Analysis of bacteriome in patients with active symptoms
library(dplyr)
library(tidyverse)
library(Maaslin2)
library(ComplexHeatmap)
library(ggpubr)
library(ggstatsplot)

#4D JP cohort
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/SP_mOTU3_JP") 
#UCvsControl
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% filter(Meta$UC_Dx == 1 | Meta$UC_Dx ==0) -> Meta 
Meta %>% select(UC_Dx, Age, sex, BMI) -> JP_meta

SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
SP_input %>% filter(SP_input$IBD ==2 |SP_input$UC ==1) -> SP_input
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("ID", "IBD", "UC", "CD"))] 

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

SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
SP_input %>% filter(SP_input$IBD ==2 |SP_input$CD ==1) -> SP_input
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("ID", "IBD", "UC", "CD"))] 

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

SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("ID", "IBD", "UC", "CD"))] 

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
JP_MA %>% mutate(ID = str_sub(JP_MA$feature, start = -6, end = -2)) -> JP_MA

#For clinical symptoms
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/SP_mOTU3(METAFは番号に修正済みIBD症状薬剤解析)")
#UCvsControl
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)

Meta %>% filter(Meta$UC_Dx == 1 & Meta$"Any.clinical.symptoms" == 0 | Meta$UC_Dx ==0) -> Meta #UC without symptoms and control
rownames_to_column(Meta, "METAF") -> Meta
Meta %>% select(UC_Dx, METAF, Age, sex, BMI) -> JP_meta

SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
left_join(Meta, SP_input, by = "METAF") -> SP_input
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("ID", "IBD_HC2", "Age", "sex", "BMI", "IBD_HC2.1", "UC.x", 
                                                          "Crohn", "UC_Dx", "CD_Dx", "IBD.x", "Immunomodulators", 
                                                          "Infliximab.or.adalimumab", "Any.clinical.symptoms", "IBD.y", "UC.y", "CD"))] 

column_to_rownames(JP_meta, "METAF") -> JP_meta
column_to_rownames(SP_input2, "METAF") -> SP_input2

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
rename(.data=JP_maaslin_UC, "Coefficient (UC w/o symptoms)" = "coef") -> JP_maaslin_UC
rename(.data=JP_maaslin_UC, "P-value (UC w/o symptoms)" = "pval") -> JP_maaslin_UC
rename(.data=JP_maaslin_UC, "Q-value (UC w/o symptoms)" = "qval") -> JP_maaslin_UC

#CDvsControl
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)

Meta %>% filter(Meta$CD_Dx == 1 & Meta$"Any.clinical.symptoms" == 0 | Meta$CD_Dx ==0) -> Meta #CD without TNF and control
rownames_to_column(Meta, "METAF") -> Meta
Meta %>% select(CD_Dx, METAF, Age, sex, BMI) -> JP_meta

SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
left_join(Meta, SP_input, by = "METAF") -> SP_input
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("ID", "IBD_HC2", "Age", "sex", "BMI", "IBD_HC2.1", "UC.x", 
                                                          "Crohn", "UC_Dx", "CD_Dx", "IBD.x", "Immunomodulators", 
                                                          "Infliximab.or.adalimumab", "Any.clinical.symptoms",  "IBD.y", "UC.y", "CD"))] 

column_to_rownames(JP_meta, "METAF") -> JP_meta
column_to_rownames(SP_input2, "METAF") -> SP_input2

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
rename(.data=JP_maaslin_CD, "Coefficient (CD w/o symptoms)" = "coef") -> JP_maaslin_CD
rename(.data=JP_maaslin_CD, "P-value (CD w/o symptoms)" = "pval") -> JP_maaslin_CD
rename(.data=JP_maaslin_CD, "Q-value (CD w/o symptoms)" = "qval") -> JP_maaslin_CD

#IBDvsControl
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)

Meta %>% filter(Meta$IBD == 1 & Meta$"Any.clinical.symptoms" == 0 | Meta$IBD ==0) -> Meta #IBD without TNF and control
rownames_to_column(Meta, "METAF") -> Meta
Meta %>% select(IBD, METAF, Age, sex, BMI) -> JP_meta

SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
left_join(Meta, SP_input, by = "METAF") -> SP_input
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("ID", "IBD_HC2", "Age", "sex", "BMI", "IBD_HC2.1", "UC.x", 
                                                          "Crohn", "UC_Dx", "CD_Dx", "IBD.x", "Immunomodulators", 
                                                          "Infliximab.or.adalimumab", "Any.clinical.symptoms",  "IBD.y", "UC.y", "CD"))] 

column_to_rownames(JP_meta, "METAF") -> JP_meta
column_to_rownames(SP_input2, "METAF") -> SP_input2

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
rename(.data=JP_maaslin_IBD, "Coefficient (IBD w/o symptoms)" = "coef") -> JP_maaslin_IBD
rename(.data=JP_maaslin_IBD, "P-value (IBD w/o symptoms)" = "pval") -> JP_maaslin_IBD
rename(.data=JP_maaslin_IBD, "Q-value (IBD w/o symptoms)" = "qval") -> JP_maaslin_IBD

#Data integration
JP_MA_sym <- full_join(full_join(JP_maaslin_IBD, JP_maaslin_UC, by = "feature"), JP_maaslin_CD, by = "feature")
JP_MA_sym %>% mutate(ID = str_sub(JP_MA_sym$feature, start = -6, end = -2)) -> JP_MA_sym

JP_MA2 <- full_join(JP_MA, JP_MA_sym, by = "ID")

#Spearman scatter plots
sp <- ggplot(JP_MA2, aes(x = JP_MA2$"Coefficient (IBD JP)", y = JP_MA2$"Coefficient (IBD w/o symptoms)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "", subtitle = "IBD", x = "Coefficient value (all IBD)", y = "Coefficient value (IBD w/o symptoms)", tag = "") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_ibd <- sp +  stat_cor(method = "spearman", label.x = -2, label.y = 2, cor.coef.name = c("rho"))

sp <- ggplot(JP_MA2, aes(x = JP_MA2$"Coefficient (UC JP)", y = JP_MA2$"Coefficient (UC w/o symptoms)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "", subtitle = "UC", x = "Coefficient value (all UC)", y = "Coefficient value (UC w/o symptoms)", tag = "") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_uc <- sp +  stat_cor(method = "spearman", label.x = -2, label.y = 2, cor.coef.name = c("rho"))

sp <- ggplot(JP_MA2, aes(x = JP_MA2$"Coefficient (CD JP)", y = JP_MA2$"Coefficient (CD w/o symptoms)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "", subtitle = "CD", x = "Coefficient value (all CD)", y = "Coefficient value (CD w/o symptoms)", tag = "") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_cd <- sp +  stat_cor(method = "spearman", label.x = -2, label.y = 2, cor.coef.name = c("rho"))

combine_plots(
  list(sp_ibd, sp_uc, sp_cd),
  plotgrid.args = list(nrow = 1),
  annotation.args = list(
    title = "Sensitivity analysis (all IBD vs IBD without clinically active symptoms)",
    caption = ""
  )
)

#Heatmap creation
#Select species with abs coef >1 and FDR<0.1 in UC or CD 
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/mOTU3_original_data")
ID <- read.csv("Feature.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
JP_MA2 %>% mutate(ID = as.numeric(JP_MA2$ID)) -> JP_MA2
JP_MA2 <- full_join(JP_MA2, ID, by ="ID")

JP_MA2 %>% filter(abs(JP_MA2$"Coefficient (CD JP)") > 1 & JP_MA2$"Q-value (CD JP)" < 0.1 | abs(JP_MA2$"Coefficient (UC JP)") > 1 & JP_MA2$"Q-value (UC JP)" < 0.1) -> SP_input
SP_input %>% arrange(-SP_input$"Coefficient (IBD JP)") -> SP_input2 

column_to_rownames(SP_input2, "Feature") -> SP_input2_bacteriome

SP_input2_bacteriome %>% select("Coefficient (IBD JP)", "P-value (IBD JP)", "Q-value (IBD JP)", "Coefficient (UC JP)", "P-value (UC JP)", "Q-value (UC JP)", "Coefficient (CD JP)", "P-value (CD JP)", "Q-value (CD JP)", 
                                "Coefficient (IBD w/o symptoms)", "P-value (IBD w/o symptoms)", "Q-value (IBD w/o symptoms)", "Coefficient (UC w/o symptoms)", "P-value (UC w/o symptoms)", "Q-value (UC w/o symptoms)", 
                                "Coefficient (CD w/o symptoms)", "P-value (CD w/o symptoms)", "Q-value (CD w/o symptoms)") ->  Fig1 

Fig1 %>% select("Q-value (IBD JP)", "Q-value (IBD w/o symptoms)") -> qval_IBD
Fig1 %>% select("Q-value (UC JP)",  "Q-value (UC w/o symptoms)") -> qval_UC
Fig1 %>% select("Q-value (CD JP)", "Q-value (CD w/o symptoms)") -> qval_CD

qval_IBD[is.na(qval_IBD)] <- 1
qval_UC[is.na(qval_UC)] <- 1
qval_CD[is.na(qval_CD)] <- 1

Fig1 %>% select("Coefficient (IBD JP)", "Coefficient (IBD w/o symptoms)") -> Gram_coef_IBD
Fig1 %>% select("Coefficient (UC JP)",  "Coefficient (UC w/o symptoms)") -> Gram_coef_UC
Fig1 %>% select("Coefficient (CD JP)", "Coefficient (CD w/o symptoms)") -> Gram_coef_CD

anno_width = unit(2, "cm")

rename(.data= Gram_coef_IBD, "All IBD" = "Coefficient (IBD JP)") -> Gram_coef_IBD
rename(.data= Gram_coef_IBD, "IBD w/o symptoms" = "Coefficient (IBD w/o symptoms)") -> Gram_coef_IBD #Franzosa_2018

rename(.data= Gram_coef_UC, "All UC" = "Coefficient (UC JP)") -> Gram_coef_UC
rename(.data= Gram_coef_UC, "UC w/o symptoms" = "Coefficient (UC w/o symptoms)") -> Gram_coef_UC #Franzosa_2018

rename(.data= Gram_coef_CD, "All CD" = "Coefficient (CD JP)") -> Gram_coef_CD
rename(.data= Gram_coef_CD, "CD w/o symptoms" = "Coefficient (CD w/o symptoms)") -> Gram_coef_CD #Franzosa_2018

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")

p1=pheatmap(as.matrix(Gram_coef_IBD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(qval_IBD < 0.1,"*", ""), nrow(qval_IBD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-2.5,0,1.5), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))
p2=pheatmap(as.matrix(Gram_coef_UC), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(qval_UC < 0.1,"*", ""), nrow(qval_UC)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-2.5,0,1.5), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))
p3=pheatmap(as.matrix(Gram_coef_CD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(qval_CD < 0.1,"*", ""), nrow(qval_CD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-2.5,0,1.5), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))

draw(p1+p2+p3, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig))


