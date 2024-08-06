#Propensity score matching analysis
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(Maaslin2)
library(MatchIt)
library(tableone)
library(cobalt)
library(ggpubr)
library(ggstatsplot)

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/PSM")
set.seed(123)

Meta <- read.csv("240612metadataforPSM.csv", header = TRUE, na.strings = c(NA, '')) 
#column_to_rownames(Meta, "Stool_No.") -> Meta

matched_1 <- matchit(IBDforPSM ~ `Age50over` + `sex` + `BMI25over` + `Alchol` + 
                       `CurrentSmoking` + `Bread`+ `Vegetables` + `Fruit` + `Seafood` + 
                       `Meat` + `Egg` + `Coffee` + `StrongActive` + `ModerateActive` + 
                       `WeekActive` +`DM` +`Dyslipidemia` + `Hypertension` + 
                       `AcuteCoronarySyndrome` + `ChronicHeartFaulure` + `CVD` + 
                       `Depression` + `AllergicDisease` + `CollagenDisease` + 
                       `BSS12` + `BSS345`+ `BSS67`  + `PPI` + `Antiinfectives`,
                     data=Meta,
                     method="nearest",
                     distance = "glm",
                     link = "logit",
                     std.caliper = TRUE,
                     caliper = 0.2,
                     estimand = "ATT",
                     replace=FALSE,
                     ratio=4,
                     seed=1234)

matched_1
summary(matched_1)
#plot(matched_1, type="jitter")
love.plot(matched_1, threshold = 0.1, abs = TRUE, stars="std")
data_ps <- match.data(matched_1) %>% as.data.frame() #Go to the CreateTableOne if table is needed

data_ps %>% dplyr::select("Stool_No", "Age", "sex", "BMI", "IBDforPSM") -> PSM_IBD

#UC for PSM
Meta <- read.csv("240612metadataforPSM.csv", header = TRUE, na.strings = c(NA, '')) 
Meta %>% filter(!is.na(Meta$"UCforPSM")) -> Meta

matched_1 <- matchit(UCforPSM ~ `Age50over` + `sex` + `BMI25over` +  `Alchol` + 
                       `CurrentSmoking` + `Bread`+ `Vegetables` + `Fruit` + `Seafood` + 
                       `Meat` + `Egg` + `Coffee` + `StrongActive` + `ModerateActive` + 
                       `WeekActive` +`DM` +`Dyslipidemia` + `Hypertension` + 
                       `AcuteCoronarySyndrome` + `ChronicHeartFaulure` + `CVD` + 
                       `Depression` + `AllergicDisease` + `CollagenDisease` + 
                       `BSS12`  + `BSS345`+  `BSS67`  + `PPI` + `Antiinfectives`,
                     data=Meta,
                     method="nearest",
                     distance = "glm",
                     link = "logit",
                     std.caliper = TRUE,
                     caliper = 0.2,
                     estimand = "ATT",
                     replace=FALSE,
                     ratio=4,
                     seed=1234)

matched_1
summary(matched_1)
#plot(matched_1, type="jitter")
love.plot(matched_1, threshold = 0.1, abs = TRUE, stars="std")
data_ps <- match.data(matched_1) %>% as.data.frame() #Go to the CreateTableOne if table is needed

data_ps %>% dplyr::select("Stool_No", "Age", "sex", "BMI", "UCforPSM") -> PSM_UC

#CD for PSM
Meta <- read.csv("240612metadataforPSM.csv", header = TRUE, na.strings = c(NA, '')) 
Meta %>% filter(!is.na(Meta$"CDforPSM")) -> Meta
matched_1 <- matchit(CDforPSM ~ `Age50over` + `sex` + `BMI25over` +  `Alchol` + 
                       `CurrentSmoking` + `Bread`+ `Vegetables` + `Fruit` + `Seafood` + 
                       `Meat` + `Egg` + `Coffee` + `StrongActive` + `ModerateActive` + 
                       `WeekActive` +`DM` +`Dyslipidemia` + `Hypertension` + 
                       `AcuteCoronarySyndrome` + `ChronicHeartFaulure` + `CVD` + 
                       `Depression` + `AllergicDisease` + `CollagenDisease` + 
                       `BSS12` + `BSS345`+  `BSS67`  + `PPI` + `Antiinfectives`,
                     data=Meta,
                     method="nearest",
                     distance = "glm",
                     link = "logit",
                     std.caliper = TRUE,
                     caliper = 0.2,
                     estimand = "ATT",
                     replace=FALSE,
                     ratio=4,
                     seed=1234)

matched_1
summary(matched_1)
#plot(matched_1, type="jitter")
love.plot(matched_1, threshold = 0.1, abs = TRUE, stars="std")
data_ps <- match.data(matched_1) %>% as.data.frame() #Go to the CreateTableOne if table is needed

data_ps %>% dplyr::select("Stool_No", "Age", "sex", "BMI", "CDforPSM") -> PSM_CD

#To create Table of PSM cohort
tableone_matched <-CreateTableOne(vars = c("Age50over", "sex", "BMI25over", "Alchol", "CurrentSmoking", "Bread", 
                                           "Vegetables", "Fruit", "Seafood", "Meat", "Egg", "Coffee", "StrongActive", 
                                           "ModerateActive", "WeekActive", "DM", "Dyslipidemia", "Hypertension", 
                                           "AcuteCoronarySyndrome", "ChronicHeartFaulure", "CVD", "Depression", 
                                           "AllergicDisease", "CollagenDisease", "BSS12", "BSS345", "BSS67",  "PPI", 
                                           "Antiinfectives"), 
                                  strata  = "IBDforPSM", data = data_ps, factorVars = c("Age50over", 
                                                                                        "sex", "BMI25over", "Alchol", "CurrentSmoking", "Bread", "Vegetables", 
                                                                                        "Fruit", "Seafood", "Meat", "Egg", "Coffee", "StrongActive", 
                                                                                        "ModerateActive", "WeekActive", "DM", "Dyslipidemia", "Hypertension", 
                                                                                        "AcuteCoronarySyndrome", "ChronicHeartFaulure", "CVD", "Depression", 
                                                                                        "AllergicDisease", "CollagenDisease", "BSS12", "BSS345", "BSS67", "PPI", 
                                                                                        "Antiinfectives"), test = TRUE) #Select IBDforPSM, UCforPSM, or CDforPSM

TestTable<- print(tableone_matched, smd = TRUE) 　#SMD<0.1 is the best condition.
write.csv(TestTable, file = "TestTable.csv", fileEncoding = "CP932")

#MaAsLin2 analysis of PSM cohort
#IBDvsControl
SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
rownames_to_column(SP_input, "Stool_No") -> SP_input
SP_input %>% mutate(Stool_No = as.numeric(SP_input$Stool_No)) -> SP_input

left_join(PSM_IBD, SP_input, by = "Stool_No") -> SP_input2

SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("Age", "sex", "BMI", "IBDforPSM", "METAF", "Selection", "IBDForPSM", "UCForPSM", "CDForPSM"))] 
column_to_rownames(SP_input2, "Stool_No") -> SP_input2

rownames(PSM_IBD) <- NULL
column_to_rownames(PSM_IBD, "Stool_No") -> JP_meta

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
  fixed_effects=c("IBDforPSM"), #No need to include Age Sex BMI as this is the PSM cohort.
  reference = c("IBDforPSM,0")) 

fit_data$results -> JP_maaslin
JP_maaslin %>% filter(JP_maaslin$"metadata" == "IBDforPSM") -> JP_maaslin_IBD
JP_maaslin_IBD <- JP_maaslin_IBD[, -which (colnames(JP_maaslin_IBD) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=JP_maaslin_IBD, "Coefficient (IBD PSM)" = "coef") -> JP_maaslin_IBD
rename(.data=JP_maaslin_IBD, "P-value (IBD PSM)" = "pval") -> JP_maaslin_IBD
rename(.data=JP_maaslin_IBD, "Q-value (IBD PSM)" = "qval") -> JP_maaslin_IBD

#UC
SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
rownames_to_column(SP_input, "Stool_No") -> SP_input
SP_input %>% mutate(Stool_No = as.numeric(SP_input$Stool_No)) -> SP_input

left_join(PSM_UC, SP_input, by = "Stool_No") -> SP_input2

SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("Age", "sex", "BMI", "UCforPSM", "METAF", "Selection", "IBDForPSM", "UCForPSM", "CDForPSM"))] 
column_to_rownames(SP_input2, "Stool_No") -> SP_input2

rownames(PSM_UC) <- NULL
column_to_rownames(PSM_UC, "Stool_No") -> JP_meta

keep <- apply(SP_input2, 2, mean) > 1E-4 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-4、prevalence > 0.1
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
  fixed_effects=c("UCforPSM"), #No need to include Age Sex BMI as this is the PSM cohort.
  reference = c("UCforPSM,0")) 

fit_data$results -> JP_maaslin
JP_maaslin %>% filter(JP_maaslin$"metadata" == "UCforPSM") -> JP_maaslin_UC
JP_maaslin_UC <- JP_maaslin_UC[, -which (colnames(JP_maaslin_UC) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=JP_maaslin_UC, "Coefficient (UC PSM)" = "coef") -> JP_maaslin_UC
rename(.data=JP_maaslin_UC, "P-value (UC PSM)" = "pval") -> JP_maaslin_UC
rename(.data=JP_maaslin_UC, "Q-value (UC PSM)" = "qval") -> JP_maaslin_UC

#CD
SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
rownames_to_column(SP_input, "Stool_No") -> SP_input
SP_input %>% mutate(Stool_No = as.numeric(SP_input$Stool_No)) -> SP_input

left_join(PSM_CD, SP_input, by = "Stool_No") -> SP_input2

SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("Age", "sex", "BMI", "CDforPSM", "METAF", "Selection", "IBDForPSM", "UCForPSM", "CDForPSM"))] 
column_to_rownames(SP_input2, "Stool_No") -> SP_input2

rownames(PSM_CD) <- NULL
column_to_rownames(PSM_CD, "Stool_No") -> JP_meta

keep <- apply(SP_input2, 2, mean) > 1E-4 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-4、prevalence > 0.1
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
  fixed_effects=c("CDforPSM"), #No need to include Age Sex BMI as this is the PSM cohort.
  reference = c("CDforPSM,0")) 

fit_data$results -> JP_maaslin
JP_maaslin %>% filter(JP_maaslin$"metadata" == "CDforPSM") -> JP_maaslin_CD
JP_maaslin_CD <- JP_maaslin_CD[, -which (colnames(JP_maaslin_CD) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=JP_maaslin_CD, "Coefficient (CD PSM)" = "coef") -> JP_maaslin_CD
rename(.data=JP_maaslin_CD, "P-value (CD PSM)" = "pval") -> JP_maaslin_CD
rename(.data=JP_maaslin_CD, "Q-value (CD PSM)" = "qval") -> JP_maaslin_CD

#Data integration
JP_MA <- full_join(full_join(JP_maaslin_IBD, JP_maaslin_UC, by = "feature"), JP_maaslin_CD, by = "feature")
JP_MA %>% mutate(ID = str_sub(JP_MA$feature, start = -6, end = -2)) -> JP_MA

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/mOTU3_original_data")
ID <- read.csv("Feature.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
JP_MA %>% mutate(ID = as.numeric(JP_MA$ID)) -> JP_MA
JP_MA <- full_join(JP_MA, ID, by ="ID")

#Heatmap creation
#Select species with abs coef>1 and FDR<0.1 in UC or CD
JP_MA %>% filter(abs(JP_MA$"Coefficient (CD PSM)") > 1 & JP_MA$"Q-value (CD PSM)" < 0.1 | abs(JP_MA$"Coefficient (UC PSM)") > 1 & JP_MA$"Q-value (UC PSM)" < 0.1) -> PSM_cohort
PSM_cohort %>% arrange(-PSM_cohort$"Coefficient (IBD PSM)") -> PSM_cohort2 
column_to_rownames(PSM_cohort2, "Feature") -> PSM_cohort2

PSM_cohort2 %>% select("Coefficient (IBD PSM)", "P-value (IBD PSM)", "Q-value (IBD PSM)", 
                       "Coefficient (UC PSM)", "P-value (UC PSM)", "Q-value (UC PSM)", 
                       "Coefficient (CD PSM)", "P-value (CD PSM)", "Q-value (CD PSM)") ->  Fig1 

Fig1 %>% select("Q-value (IBD PSM)") -> qval_IBD
Fig1 %>% select("Q-value (UC PSM)", "Q-value (CD PSM)") -> qval_UC_CD
qval_IBD[is.na(qval_IBD)] <- 1
qval_UC_CD[is.na(qval_UC_CD)] <- 1

Fig1 %>% select("Coefficient (IBD PSM)") -> Gram_coef_IBD
Fig1 %>% select("Coefficient (UC PSM)", "Coefficient (CD PSM)") -> Gram_coef_UC_CD

rename(.data= Gram_coef_IBD, "IBD PSM" = "Coefficient (IBD PSM)") -> Gram_coef_IBD
rename(.data= Gram_coef_UC_CD, "UC PSM" = "Coefficient (UC PSM)") -> Gram_coef_UC_CD 
rename(.data= Gram_coef_UC_CD, "CD PSM" = "Coefficient (CD PSM)") -> Gram_coef_UC_CD

anno_width = unit(2, "cm")

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")

p1=pheatmap(as.matrix(Gram_coef_IBD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(qval_IBD < 0.1,"*", ""), nrow(qval_IBD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-2.5,0,1.5), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))
p2=pheatmap(as.matrix(Gram_coef_UC_CD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(qval_UC_CD < 0.1,"*", ""), nrow(qval_UC_CD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-2.5,0,1.5), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))

draw(p1+p2, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig))

#Spearman correlation analysis to compare the PSM and main analysis. 
full_join(All_country, JP_MA, by = "ID") -> All_country #This dataframe is originated from the prior main analysis (Ref. Bacteriome.R)

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (IBD JP)", y = All_country$"Coefficient (IBD PSM)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Main analysis vs PSM analysis", subtitle = "IBD", x = "Coefficient value (Main analysis)", y = "Coefficient value (PSM analysis)", tag = "") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_ibd <- sp + stat_cor(method = "spearman", label.x = -1, label.y = 1, cor.coef.name = c("rho"))

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (UC JP)", y = All_country$"Coefficient (UC PSM)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Main analysis vs PSM analysis", subtitle = "UC", x = "Coefficient value (Main analysis)", y = "Coefficient value (PSM analysis)", tag = "") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_uc <- sp + stat_cor(method = "spearman", label.x = -1, label.y = 1, cor.coef.name = c("rho"))

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (CD JP)", y = All_country$"Coefficient (CD PSM)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Main analysis vs PSM analysis", subtitle = "CD", x = "Coefficient value (Main analysis)", y = "Coefficient value (PSM analysis)", tag = "") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_cd <- sp + stat_cor(method = "spearman", label.x = -1, label.y = 1, cor.coef.name = c("rho"))

combine_plots(
  list(sp_ibd, sp_uc, sp_cd),
  plotgrid.args = list(nrow = 1),
  annotation.args = list(
    title = "Correlation between the PSM and main analyses",
    caption = ""
  )
)
