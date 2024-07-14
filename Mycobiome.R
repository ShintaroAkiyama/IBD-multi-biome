library(dplyr)
library(tidyverse)
library(Maaslin2)
library(ComplexHeatmap)
library(ggpubr)
library(ggstatsplot)

#MaAsLin Mycobiome analysis (World data)
X <- 1E-10 # mean abundance > X
Y <- 0.005 # prevalence > Y

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Fungi")
SP_WO3<- read.csv("eukaryote.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = FALSE)
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3

Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")

#He_2017_Crohn (China)
WO3 %>% filter(WO3$study == "He_2017_Crohn" & WO3$timepoint == 0) -> China #timepoint = 0
China %>% distinct(subject_id, .keep_all=TRUE) -> China #Remove duplicates

#Meta
China %>% select(sample_id, disease, Sex, Age, bmi) -> China_meta
column_to_rownames(China_meta, "sample_id") -> China_meta
China_meta$"Sex"[which(China_meta$"Sex" == "male", TRUE)] <- 1
China_meta$"Sex"[which(China_meta$"Sex" == "female", TRUE)] <- 0
China_meta$"disease"[which(China_meta$"disease" == "Crohn's disease", TRUE)] <- 1
China_meta$"disease"[which(China_meta$"disease" == "Control", TRUE)] <- 0
China_meta <-apply(China_meta,c(1:2),as.numeric) %>% as.data.frame 

#SP
China_SP <- China[, -which (colnames(China) %in% c("subject_id", "environment_material", "timepoint", 
                                                   "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                                   "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(China_SP, "sample_id") -> China_SP

China_SP <- na.omit(China_SP) #Omit all NA, metaデータでも同じ操作がいるのかは確認。西嶋先生プログラムエラーに対応するための処置。

keep <- apply(China_SP, 2, mean) > X & apply(China_SP > 0, 2, sum) / nrow(China_SP) > Y # mean abundance > X、Prevalence > Y
China_SP  <- China_SP[, keep]

#Remove non-human data
as.data.frame(t(China_SP)) -> China_SP
rownames_to_column(China_SP, var = "feature") -> China_SP
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(China_SP, Human, by = "feature") -> China_SP
China_SP %>% filter(China_SP$human == 1) -> China_SP # include only fungi or protozoa reported in human 
column_to_rownames(China_SP, var="feature") -> China_SP
China_SP <- China_SP[, -which (colnames(China_SP) %in% c("human"))]
as.data.frame(t(China_SP)) -> China_SP

fit_data = Maaslin2(
  input_data = China_SP, 
  input_metadata = China_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "China_CD", 
  normalization = "NONE",　
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  standardize = FALSE,　
  cores = 4,
  fixed_effects=c("disease,Age,Sex,bmi"),
  reference = c("disease,0")) 

fit_data$results -> China_maaslin
China_maaslin %>% filter(China_maaslin$"metadata" == "disease") -> China_maaslin_CD
China_maaslin_CD <- China_maaslin_CD[, -which (colnames(China_maaslin_CD) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=China_maaslin_CD, "Coefficient (CD China)" = "coef") -> China_maaslin_CD
rename(.data=China_maaslin_CD, "P-value (CD China)" = "pval") -> China_maaslin_CD
rename(.data=China_maaslin_CD, "Q-value (CD China)" = "qval") -> China_maaslin_CD

#ES_MH, MetaHIT (Spain)
#UCvsControl
WO3 %>% filter(WO3$study == "ES_MH, MetaHIT" & WO3$timepoint == 0) -> Spain #timepoint = 0
Spain %>% distinct(subject_id, .keep_all=TRUE) -> Spain #Remove duplicates

Spain %>% filter(Spain$disease == "Ulcerative colitis"| Spain$disease == "Control") -> Spain
Spain %>% select(sample_id, disease, Sex, Age, bmi) -> Spain_meta
column_to_rownames(Spain_meta, "sample_id") -> Spain_meta
Spain_meta$"Sex"[which(Spain_meta$"Sex" == "male", TRUE)] <- 1
Spain_meta$"Sex"[which(Spain_meta$"Sex" == "female", TRUE)] <- 0
Spain_meta$"disease"[which(Spain_meta$"disease" == "Ulcerative colitis", TRUE)] <- 1
Spain_meta$"disease"[which(Spain_meta$"disease" == "Control", TRUE)] <- 0
Spain_meta <- apply(Spain_meta,c(1:2),as.numeric) %>% as.data.frame

#SP
Spain_SP <- Spain[, -which (colnames(Spain) %in% c("subject_id", "environment_material", "timepoint", 
                                                   "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                                   "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(Spain_SP, "sample_id") -> Spain_SP

keep <- apply(Spain_SP, 2, mean) > X & apply(Spain_SP > 0, 2, sum) / nrow(Spain_SP) > Y # mean abundance > X、Prevalence > Y
Spain_SP  <- Spain_SP[, keep]

#Remove non-human data
as.data.frame(t(Spain_SP)) -> Spain_SP
rownames_to_column(Spain_SP, var = "feature") -> Spain_SP
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(Spain_SP, Human, by = "feature") -> Spain_SP
Spain_SP %>% filter(Spain_SP$human == 1) -> Spain_SP # include only fungi or protozoa reported in human 
column_to_rownames(Spain_SP, var="feature") -> Spain_SP
Spain_SP <- Spain_SP[, -which (colnames(Spain_SP) %in% c("human"))]
as.data.frame(t(Spain_SP)) -> Spain_SP

fit_data = Maaslin2(
  input_data = Spain_SP, 
  input_metadata = Spain_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "Spain_UC", 
  normalization = "NONE",
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  standardize = FALSE,　
  cores = 4,
  fixed_effects=c("disease,Age,Sex,bmi"),
  reference = c("disease,0")) 

fit_data$results -> Spain_maaslin
Spain_maaslin %>% filter(Spain_maaslin$"metadata" == "disease") -> Spain_maaslin_UC
Spain_maaslin_UC <- Spain_maaslin_UC[, -which (colnames(Spain_maaslin_UC) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=Spain_maaslin_UC, "Coefficient (UC Spain)" = "coef") -> Spain_maaslin_UC
rename(.data=Spain_maaslin_UC, "P-value (UC Spain)" = "pval") -> Spain_maaslin_UC
rename(.data=Spain_maaslin_UC, "Q-value (UC Spain)" = "qval") -> Spain_maaslin_UC

#CDvsControl
WO3 %>% filter(WO3$study == "ES_MH, MetaHIT" & WO3$timepoint == 0) -> Spain #timepoint =0 
Spain %>% distinct(subject_id, .keep_all=TRUE) -> Spain #Remove duplicates

Spain %>% filter(Spain$disease == "Crohn's disease"| Spain$disease == "Control") -> Spain
Spain %>% select(sample_id, disease, Sex, Age, bmi) -> Spain_meta
column_to_rownames(Spain_meta, "sample_id") -> Spain_meta
Spain_meta$"Sex"[which(Spain_meta$"Sex" == "male", TRUE)] <- 1
Spain_meta$"Sex"[which(Spain_meta$"Sex" == "female", TRUE)] <- 0
Spain_meta$"disease"[which(Spain_meta$"disease" == "Crohn's disease", TRUE)] <- 1
Spain_meta$"disease"[which(Spain_meta$"disease" == "Control", TRUE)] <- 0
Spain_meta <- apply(Spain_meta,c(1:2),as.numeric) %>% as.data.frame

#SP
Spain_SP <- Spain[, -which (colnames(Spain) %in% c("subject_id", "environment_material", "timepoint", 
                                                   "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                                   "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(Spain_SP, "sample_id") -> Spain_SP

keep <- apply(Spain_SP, 2, mean) > X & apply(Spain_SP > 0, 2, sum) / nrow(Spain_SP) > Y # mean abundance > X、Prevalence > Y
Spain_SP  <- Spain_SP[, keep]

#Remove non-human data
as.data.frame(t(Spain_SP)) -> Spain_SP
rownames_to_column(Spain_SP, var = "feature") -> Spain_SP
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(Spain_SP, Human, by = "feature") -> Spain_SP
Spain_SP %>% filter(Spain_SP$human == 1) -> Spain_SP # include only fungi or protozoa reported in human 
column_to_rownames(Spain_SP, var="feature") -> Spain_SP
Spain_SP <- Spain_SP[, -which (colnames(Spain_SP) %in% c("human"))]
as.data.frame(t(Spain_SP)) -> Spain_SP

fit_data = Maaslin2(
  input_data = Spain_SP, 
  input_metadata = Spain_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "Spain_CD", 
  normalization = "NONE",
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  standardize = FALSE,　
  cores = 4,
  fixed_effects=c("disease,Age,Sex,bmi"),
  reference = c("disease,0")) 

fit_data$results -> Spain_maaslin
Spain_maaslin %>% filter(Spain_maaslin$"metadata" == "disease") -> Spain_maaslin_CD
Spain_maaslin_CD <- Spain_maaslin_CD[, -which (colnames(Spain_maaslin_CD) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=Spain_maaslin_CD, "Coefficient (CD Spain)" = "coef") -> Spain_maaslin_CD
rename(.data=Spain_maaslin_CD, "P-value (CD Spain)" = "pval") -> Spain_maaslin_CD
rename(.data=Spain_maaslin_CD, "Q-value (CD Spain)" = "qval") -> Spain_maaslin_CD

#IBDvsControl 
WO3 %>% filter(WO3$study == "ES_MH, MetaHIT" & WO3$timepoint == 0) -> Spain #timepoint = 0
Spain %>% distinct(subject_id, .keep_all=TRUE) -> Spain #Remove duplicates

Spain %>% filter(Spain$disease == "Ulcerative colitis"| Spain$disease == "Control" | Spain$disease == "Crohn's disease") -> Spain
Spain %>% select(sample_id, disease, Sex, Age) -> Spain_meta
column_to_rownames(Spain_meta, "sample_id") -> Spain_meta
Spain_meta$"Sex"[which(Spain_meta$"Sex" == "male", TRUE)] <- 1
Spain_meta$"Sex"[which(Spain_meta$"Sex" == "female", TRUE)] <- 0
Spain_meta$"disease"[which(Spain_meta$"disease" == "Ulcerative colitis" | Spain_meta$"disease" == "Crohn's disease", TRUE)] <- 1
Spain_meta$"disease"[which(Spain_meta$"disease" == "Control", TRUE)] <- 0
Spain_meta <- apply(Spain_meta,c(1:2),as.numeric) %>% as.data.frame

#SP
Spain_SP <- Spain[, -which (colnames(Spain) %in% c("subject_id", "environment_material", "timepoint", 
                                                   "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                                   "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(Spain_SP, "sample_id") -> Spain_SP

keep <- apply(Spain_SP, 2, mean) > X & apply(Spain_SP > 0, 2, sum) / nrow(Spain_SP) > Y # mean abundance > X、Prevalence > Y
Spain_SP  <- Spain_SP[, keep]

#Remove non-human data
as.data.frame(t(Spain_SP)) -> Spain_SP
rownames_to_column(Spain_SP, var = "feature") -> Spain_SP
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(Spain_SP, Human, by = "feature") -> Spain_SP
Spain_SP %>% filter(Spain_SP$human == 1) -> Spain_SP # include only fungi or protozoa reported in human 
column_to_rownames(Spain_SP, var="feature") -> Spain_SP
Spain_SP <- Spain_SP[, -which (colnames(Spain_SP) %in% c("human"))]
as.data.frame(t(Spain_SP)) -> Spain_SP

fit_data = Maaslin2(
  input_data = Spain_SP, 
  input_metadata = Spain_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "Spain_IBD", 
  normalization = "NONE",
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  standardize = FALSE,　
  cores = 4,
  fixed_effects=c("disease,Age,Sex,bmi"),
  reference = c("disease,0")) 

fit_data$results -> Spain_maaslin
Spain_maaslin %>% filter(Spain_maaslin$"metadata" == "disease") -> Spain_maaslin_IBD
Spain_maaslin_IBD <- Spain_maaslin_IBD[, -which (colnames(Spain_maaslin_IBD) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=Spain_maaslin_IBD, "Coefficient (IBD Spain)" = "coef") -> Spain_maaslin_IBD
rename(.data=Spain_maaslin_IBD, "P-value (IBD Spain)" = "pval") -> Spain_maaslin_IBD
rename(.data=Spain_maaslin_IBD, "Q-value (IBD Spain)" = "qval") -> Spain_maaslin_IBD

#Strazar_2021_Tanzania
#UCvsControl
WO3 %>% filter(WO3$study == "Strazar_2021_Tanzania" & WO3$timepoint == 0) -> Tanzania #timepoint = 0 
Tanzania %>% distinct(subject_id, .keep_all=TRUE) -> Tanzania #Remove duplicates

Tanzania %>% filter(Tanzania$disease == "Ulcerative colitis"| Tanzania$disease == "Control") -> Tanzania
Tanzania %>% select(sample_id, disease, Sex, Age, bmi) -> Tanzania_meta
column_to_rownames(Tanzania_meta, "sample_id") -> Tanzania_meta
Tanzania_meta$"Sex"[which(Tanzania_meta$"Sex" == "male", TRUE)] <- 1
Tanzania_meta$"Sex"[which(Tanzania_meta$"Sex" == "female", TRUE)] <- 0
Tanzania_meta$"disease"[which(Tanzania_meta$"disease" == "Ulcerative colitis", TRUE)] <- 1
Tanzania_meta$"disease"[which(Tanzania_meta$"disease" == "Control", TRUE)] <- 0
Tanzania_meta <- apply(Tanzania_meta,c(1:2),as.numeric) %>% as.data.frame

#SP
Tanzania_SP <- Tanzania[, -which (colnames(Tanzania) %in% c("subject_id", "environment_material", "timepoint", 
                                                            "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                                            "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(Tanzania_SP, "sample_id") -> Tanzania_SP

keep <- apply(Tanzania_SP, 2, mean) > X & apply(Tanzania_SP > 0, 2, sum) / nrow(Tanzania_SP) > Y # mean abundance > X、Prevalence > Y
Tanzania_SP  <- Tanzania_SP[, keep]

#Remove non-human data
as.data.frame(t(Tanzania_SP)) -> Tanzania_SP
rownames_to_column(Tanzania_SP, var = "feature") -> Tanzania_SP
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(Tanzania_SP, Human, by = "feature") -> Tanzania_SP
Tanzania_SP %>% filter(Tanzania_SP$human == 1) -> Tanzania_SP # include only fungi or protozoa reported in human 
column_to_rownames(Tanzania_SP, var="feature") -> Tanzania_SP
Tanzania_SP <- Tanzania_SP[, -which (colnames(Tanzania_SP) %in% c("human"))]
as.data.frame(t(Tanzania_SP)) -> Tanzania_SP

fit_data = Maaslin2(
  input_data = Tanzania_SP, 
  input_metadata = Tanzania_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "Tanzania_UC", 
  normalization = "NONE",
  transform = "LOG",
  standardize = FALSE,　
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  cores = 4,
  fixed_effects=c("disease,Age,Sex,bmi"),
  reference = c("disease,0")) 

fit_data$results -> Tanzania_maaslin
Tanzania_maaslin %>% filter(Tanzania_maaslin$"metadata" == "disease") -> Tanzania_maaslin_UC
Tanzania_maaslin_UC <- Tanzania_maaslin_UC[, -which (colnames(Tanzania_maaslin_UC) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=Tanzania_maaslin_UC, "Coefficient (UC Tanzania)" = "coef") -> Tanzania_maaslin_UC
rename(.data=Tanzania_maaslin_UC, "P-value (UC Tanzania)" = "pval") -> Tanzania_maaslin_UC
rename(.data=Tanzania_maaslin_UC, "Q-value (UC Tanzania)" = "qval") -> Tanzania_maaslin_UC

#Lloyd-Price_2019_HMP2IBD (US data)
#UCvsControl
WO3 %>% filter(WO3$study == "Lloyd-Price_2019_HMP2IBD") -> US #No specific timepoint 
US %>% distinct(subject_id, .keep_all=TRUE) -> US #Remove duplicates

US %>% filter(US$disease == "Ulcerative colitis"| US$disease == "Control") -> US
US %>% select(sample_id, disease, Sex, Age) -> US_meta
column_to_rownames(US_meta, "sample_id") -> US_meta
US_meta$"Sex"[which(US_meta$"Sex" == "male", TRUE)] <- 1
US_meta$"Sex"[which(US_meta$"Sex" == "female", TRUE)] <- 0
US_meta$"disease"[which(US_meta$"disease" == "Ulcerative colitis", TRUE)] <- 1
US_meta$"disease"[which(US_meta$"disease" == "Control", TRUE)] <- 0
US_meta <- apply(US_meta,c(1:2),as.numeric) %>% as.data.frame

#SP
US_SP <- US[, -which (colnames(US) %in% c("subject_id", "environment_material", "timepoint", 
                                          "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                          "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US_SP, "sample_id") -> US_SP

keep <- apply(US_SP, 2, mean) > X & apply(US_SP > 0, 2, sum) / nrow(US_SP) > Y # mean abundance > X、Prevalence > Y
US_SP  <- US_SP[, keep]

#Remove non-human data
as.data.frame(t(US_SP)) -> US_SP
rownames_to_column(US_SP, var = "feature") -> US_SP
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(US_SP, Human, by = "feature") -> US_SP
US_SP %>% filter(US_SP$human == 1) -> US_SP # include only fungi or protozoa reported in human 
column_to_rownames(US_SP, var="feature") -> US_SP
US_SP <- US_SP[, -which (colnames(US_SP) %in% c("human"))]
as.data.frame(t(US_SP)) -> US_SP

fit_data = Maaslin2(
  input_data = US_SP, 
  input_metadata = US_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "US_UC", 
  normalization = "NONE",
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  standardize = FALSE,　
  cores = 4,
  fixed_effects=c("disease,Age,Sex"),
  reference = c("disease,0")) 

fit_data$results -> US_maaslin
US_maaslin %>% filter(US_maaslin$"metadata" == "disease") -> US_maaslin_UC
US_maaslin_UC <- US_maaslin_UC[, -which (colnames(US_maaslin_UC) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=US_maaslin_UC, "Coefficient (UC US)" = "coef") -> US_maaslin_UC
rename(.data=US_maaslin_UC, "P-value (UC US)" = "pval") -> US_maaslin_UC
rename(.data=US_maaslin_UC, "Q-value (UC US)" = "qval") -> US_maaslin_UC

#CDvsControl
WO3 %>% filter(WO3$study == "Lloyd-Price_2019_HMP2IBD") -> US #No specific timepoint 
US %>% distinct(subject_id, .keep_all=TRUE) -> US #Remove duplicates

US %>% filter(US$disease == "Crohn's disease"| US$disease == "Control") -> US
US %>% select(sample_id, disease, Sex, Age) -> US_meta
column_to_rownames(US_meta, "sample_id") -> US_meta
US_meta$"Sex"[which(US_meta$"Sex" == "male", TRUE)] <- 1
US_meta$"Sex"[which(US_meta$"Sex" == "female", TRUE)] <- 0
US_meta$"disease"[which(US_meta$"disease" == "Crohn's disease", TRUE)] <- 1
US_meta$"disease"[which(US_meta$"disease" == "Control", TRUE)] <- 0
US_meta <- apply(US_meta,c(1:2),as.numeric) %>% as.data.frame

#SP
US_SP <- US[, -which (colnames(US) %in% c("subject_id", "environment_material", "timepoint", 
                                          "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                          "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US_SP, "sample_id") -> US_SP

keep <- apply(US_SP, 2, mean) > X & apply(US_SP > 0, 2, sum) / nrow(US_SP) > Y # mean abundance > X、Prevalence > Y
US_SP  <- US_SP[, keep]

#Remove non-human data
as.data.frame(t(US_SP)) -> US_SP
rownames_to_column(US_SP, var = "feature") -> US_SP
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(US_SP, Human, by = "feature") -> US_SP
US_SP %>% filter(US_SP$human == 1) -> US_SP # include only fungi or protozoa reported in human 
column_to_rownames(US_SP, var="feature") -> US_SP
US_SP <- US_SP[, -which (colnames(US_SP) %in% c("human"))]
as.data.frame(t(US_SP)) -> US_SP

fit_data = Maaslin2(
  input_data = US_SP, 
  input_metadata = US_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "US_CD", 
  normalization = "NONE",
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  standardize = FALSE,　
  cores = 4,
  fixed_effects=c("disease,Age,Sex"),
  reference = c("disease,0")) 

fit_data$results -> US_maaslin
US_maaslin %>% filter(US_maaslin$"metadata" == "disease") -> US_maaslin_CD
US_maaslin_CD <- US_maaslin_CD[, -which (colnames(US_maaslin_CD) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=US_maaslin_CD, "Coefficient (CD US)" = "coef") -> US_maaslin_CD
rename(.data=US_maaslin_CD, "P-value (CD US)" = "pval") -> US_maaslin_CD
rename(.data=US_maaslin_CD, "Q-value (CD US)" = "qval") -> US_maaslin_CD

#IBDvsControl 
WO3 %>% filter(WO3$study == "Lloyd-Price_2019_HMP2IBD") -> US #No specific timepoint
US %>% distinct(subject_id, .keep_all=TRUE) -> US #Remove duplicates

US %>% filter(US$disease == "Ulcerative colitis"| US$disease == "Control" | US$disease == "Crohn's disease") -> US
US %>% select(sample_id, disease, Sex, Age) -> US_meta
column_to_rownames(US_meta, "sample_id") -> US_meta
US_meta$"Sex"[which(US_meta$"Sex" == "male", TRUE)] <- 1
US_meta$"Sex"[which(US_meta$"Sex" == "female", TRUE)] <- 0
US_meta$"disease"[which(US_meta$"disease" == "Ulcerative colitis" | US_meta$"disease" == "Crohn's disease", TRUE)] <- 1
US_meta$"disease"[which(US_meta$"disease" == "Control", TRUE)] <- 0
US_meta <- apply(US_meta,c(1:2),as.numeric) %>% as.data.frame

#SP
US_SP <- US[, -which (colnames(US) %in% c("subject_id", "environment_material", "timepoint", 
                                          "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                          "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US_SP, "sample_id") -> US_SP

keep <- apply(US_SP, 2, mean) > X & apply(US_SP > 0, 2, sum) / nrow(US_SP) > Y # mean abundance > X、Prevalence > Y
US_SP  <- US_SP[, keep]

#Remove non-human data
as.data.frame(t(US_SP)) -> US_SP
rownames_to_column(US_SP, var = "feature") -> US_SP
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(US_SP, Human, by = "feature") -> US_SP
US_SP %>% filter(US_SP$human == 1) -> US_SP # include only fungi or protozoa reported in human
column_to_rownames(US_SP, var="feature") -> US_SP
US_SP <- US_SP[, -which (colnames(US_SP) %in% c("human"))]
as.data.frame(t(US_SP)) -> US_SP

fit_data = Maaslin2(
  input_data = US_SP, 
  input_metadata = US_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "US_IBD", 
  normalization = "NONE",
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  standardize = FALSE,
  cores = 4,
  fixed_effects=c("disease,Age,Sex"),
  reference = c("disease,0")) 

fit_data$results -> US_maaslin
US_maaslin %>% filter(US_maaslin$"metadata" == "disease") -> US_maaslin_IBD
US_maaslin_IBD <- US_maaslin_IBD[, -which (colnames(US_maaslin_IBD) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=US_maaslin_IBD, "Coefficient (IBD US)" = "coef") -> US_maaslin_IBD
rename(.data=US_maaslin_IBD, "P-value (IBD US)" = "pval") -> US_maaslin_IBD
rename(.data=US_maaslin_IBD, "Q-value (IBD US)" = "qval") -> US_maaslin_IBD

#Franzosa_2018_IBD (US cohort)
#UCvsControl
WO3 %>% filter(WO3$study == "Franzosa_2018_IBD" & WO3$geographic_location == "USA") -> US2 #timepoints are all zero in this study 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #Remove duplicates

US2 %>% filter(US2$disease == "Ulcerative colitis"| US2$disease == "Control") -> US2
US2 %>% select(sample_id, disease, Sex, Age) -> US2_meta
column_to_rownames(US2_meta, "sample_id") -> US2_meta
US2_meta$"disease"[which(US2_meta$"disease" == "Ulcerative colitis", TRUE)] <- 1
US2_meta$"disease"[which(US2_meta$"disease" == "Control", TRUE)] <- 0
US2_meta <- apply(US2_meta,c(1:2),as.numeric) %>% as.data.frame

#SP
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP

keep <- apply(US2_SP, 2, mean) > X & apply(US2_SP > 0, 2, sum) / nrow(US2_SP) > Y # mean abundance > X、Prevalence > Y
US2_SP  <- US2_SP[, keep]

#Remove non-human data
as.data.frame(t(US2_SP)) -> US2_SP
rownames_to_column(US2_SP, var = "feature") -> US2_SP
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(US2_SP, Human, by = "feature") -> US2_SP
US2_SP %>% filter(US2_SP$human == 1) -> US2_SP # include only fungi or protozoa reported in human 
column_to_rownames(US2_SP, var="feature") -> US2_SP
US2_SP <- US2_SP[, -which (colnames(US2_SP) %in% c("human"))]
as.data.frame(t(US2_SP)) -> US2_SP

fit_data = Maaslin2(
  input_data = US2_SP, 
  input_metadata = US2_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "US2_UC", 
  normalization = "NONE",
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  standardize = FALSE,
  cores = 4,
  fixed_effects=c("disease,Age"),
  reference = c("disease,0")) 

fit_data$results -> US2_maaslin
US2_maaslin %>% filter(US2_maaslin$"metadata" == "disease") -> US2_maaslin_UC
US2_maaslin_UC <- US2_maaslin_UC[, -which (colnames(US2_maaslin_UC) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=US2_maaslin_UC, "Coefficient (UC US2)" = "coef") -> US2_maaslin_UC
rename(.data=US2_maaslin_UC, "P-value (UC US2)" = "pval") -> US2_maaslin_UC
rename(.data=US2_maaslin_UC, "Q-value (UC US2)" = "qval") -> US2_maaslin_UC

#CDvsControl 
WO3 %>% filter(WO3$study == "Franzosa_2018_IBD" & WO3$geographic_location == "USA") -> US2 #timepoints are all zero in this study 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #Remove duplicates

US2 %>% filter(US2$disease == "Crohn's disease"| US2$disease == "Control") -> US2
US2 %>% select(sample_id, disease, Sex, Age) -> US2_meta
column_to_rownames(US2_meta, "sample_id") -> US2_meta
US2_meta$"disease"[which(US2_meta$"disease" == "Crohn's disease", TRUE)] <- 1
US2_meta$"disease"[which(US2_meta$"disease" == "Control", TRUE)] <- 0
US2_meta <- apply(US2_meta,c(1:2),as.numeric) %>% as.data.frame

#SP
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP

keep <- apply(US2_SP, 2, mean) > X & apply(US2_SP > 0, 2, sum) / nrow(US2_SP) > Y # mean abundance > X、Prevalence > Y
US2_SP  <- US2_SP[, keep]

#Remove non-human data
as.data.frame(t(US2_SP)) -> US2_SP
rownames_to_column(US2_SP, var = "feature") -> US2_SP
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(US2_SP, Human, by = "feature") -> US2_SP
US2_SP %>% filter(US2_SP$human == 1) -> US2_SP # include only fungi or protozoa reported in human 
column_to_rownames(US2_SP, var="feature") -> US2_SP
US2_SP <- US2_SP[, -which (colnames(US2_SP) %in% c("human"))]
as.data.frame(t(US2_SP)) -> US2_SP

fit_data = Maaslin2(
  input_data = US2_SP, 
  input_metadata = US2_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "US2_CD", 
  normalization = "NONE",
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  standardize = FALSE,
  cores = 4,
  fixed_effects=c("disease,Age"),
  reference = c("disease,0")) 

fit_data$results -> US2_maaslin
US2_maaslin %>% filter(US2_maaslin$"metadata" == "disease") -> US2_maaslin_CD
US2_maaslin_CD <- US2_maaslin_CD[, -which (colnames(US2_maaslin_CD) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=US2_maaslin_CD, "Coefficient (CD US2)" = "coef") -> US2_maaslin_CD
rename(.data=US2_maaslin_CD, "P-value (CD US2)" = "pval") -> US2_maaslin_CD
rename(.data=US2_maaslin_CD, "Q-value (CD US2)" = "qval") -> US2_maaslin_CD

#IBDvsControl 
WO3 %>% filter(WO3$study == "Franzosa_2018_IBD" & WO3$geographic_location == "USA") -> US2 #timepoint=0
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #Remove duplicates

US2 %>% filter(US2$disease == "Ulcerative colitis"| US2$disease == "Control" | US2$disease == "Crohn's disease") -> US2
US2 %>% select(sample_id, disease, Sex, Age) -> US2_meta
column_to_rownames(US2_meta, "sample_id") -> US2_meta
US2_meta$"disease"[which(US2_meta$"disease" == "Ulcerative colitis" | US2_meta$"disease" == "Crohn's disease", TRUE)] <- 1
US2_meta$"disease"[which(US2_meta$"disease" == "Control", TRUE)] <- 0
US2_meta <- apply(US2_meta,c(1:2),as.numeric) %>% as.data.frame

#SP
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP

keep <- apply(US2_SP, 2, mean) > X & apply(US2_SP > 0, 2, sum) / nrow(US2_SP) > Y # mean abundance > X、Prevalence > Y
US2_SP <- US2_SP[, keep]

#Remove non-human data
as.data.frame(t(US2_SP)) -> US2_SP
rownames_to_column(US2_SP, var = "feature") -> US2_SP
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(US2_SP, Human, by = "feature") -> US2_SP
US2_SP %>% filter(US2_SP$human == 1) -> US2_SP # include only fungi or protozoa reported in human 
column_to_rownames(US2_SP, var="feature") -> US2_SP
US2_SP <- US2_SP[, -which (colnames(US2_SP) %in% c("human"))]
as.data.frame(t(US2_SP)) -> US2_SP

fit_data = Maaslin2(
  input_data = US2_SP, 
  input_metadata = US2_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "US2_IBD", 
  normalization = "NONE",
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  standardize = FALSE,
  cores = 4,
  fixed_effects=c("disease,Age"),
  reference = c("disease,0")) 

fit_data$results -> US2_maaslin
US2_maaslin %>% filter(US2_maaslin$"metadata" == "disease") -> US2_maaslin_IBD
US2_maaslin_IBD <- US2_maaslin_IBD[, -which (colnames(US2_maaslin_IBD) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=US2_maaslin_IBD, "Coefficient (IBD US2)" = "coef") -> US2_maaslin_IBD
rename(.data=US2_maaslin_IBD, "P-value (IBD US2)" = "pval") -> US2_maaslin_IBD
rename(.data=US2_maaslin_IBD, "Q-value (IBD US2)" = "qval") -> US2_maaslin_IBD

#Franzosa_2018_IBD (Netherlands cohort)
#UCvsControls
WO3 %>% filter(WO3$study == "Franzosa_2018_IBD" & WO3$geographic_location == "Netherlands") -> Netherlands #timepoints are all zero in this study 
Netherlands %>% distinct(subject_id, .keep_all=TRUE) -> Netherlands #Remove duplicates

Netherlands %>% filter(Netherlands$disease == "Ulcerative colitis"| Netherlands$disease == "Control") -> Netherlands
Netherlands %>% select(sample_id, disease, Sex, Age) -> Netherlands_meta
column_to_rownames(Netherlands_meta, "sample_id") -> Netherlands_meta
Netherlands_meta$"disease"[which(Netherlands_meta$"disease" == "Ulcerative colitis", TRUE)] <- 1
Netherlands_meta$"disease"[which(Netherlands_meta$"disease" == "Control", TRUE)] <- 0
Netherlands_meta <- apply(Netherlands_meta,c(1:2),as.numeric) %>% as.data.frame

#SP
Netherlands_SP <- Netherlands[, -which (colnames(Netherlands) %in% c("subject_id", "environment_material", "timepoint", 
                                                                     "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                                                     "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(Netherlands_SP, "sample_id") -> Netherlands_SP

keep <- apply(Netherlands_SP, 2, mean) > X & apply(Netherlands_SP > 0, 2, sum) / nrow(Netherlands_SP) > Y # mean abundance > X、Prevalence > Y
Netherlands_SP <- Netherlands_SP[, keep]

#Remove non-human data
as.data.frame(t(Netherlands_SP)) -> Netherlands_SP
rownames_to_column(Netherlands_SP, var = "feature") -> Netherlands_SP
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(Netherlands_SP, Human, by = "feature") -> Netherlands_SP
Netherlands_SP %>% filter(Netherlands_SP$human == 1) -> Netherlands_SP # include only fungi or protozoa reported in human 
column_to_rownames(Netherlands_SP, var="feature") -> Netherlands_SP
Netherlands_SP <- Netherlands_SP[, -which (colnames(Netherlands_SP) %in% c("human"))]
as.data.frame(t(Netherlands_SP)) -> Netherlands_SP

fit_data = Maaslin2(
  input_data = Netherlands_SP, 
  input_metadata = Netherlands_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "Netherlands_UC", 
  normalization = "NONE",
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  standardize = FALSE,
  cores = 4,
  fixed_effects=c("disease,Age"),
  reference = c("disease,0")) 

fit_data$results -> Netherlands_maaslin
Netherlands_maaslin %>% filter(Netherlands_maaslin$"metadata" == "disease") -> Netherlands_maaslin_UC
Netherlands_maaslin_UC <- Netherlands_maaslin_UC[, -which (colnames(Netherlands_maaslin_UC) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=Netherlands_maaslin_UC, "Coefficient (UC Netherlands)" = "coef") -> Netherlands_maaslin_UC
rename(.data=Netherlands_maaslin_UC, "P-value (UC Netherlands)" = "pval") -> Netherlands_maaslin_UC
rename(.data=Netherlands_maaslin_UC, "Q-value (UC Netherlands)" = "qval") -> Netherlands_maaslin_UC

#CDvsControl 
WO3 %>% filter(WO3$study == "Franzosa_2018_IBD" & WO3$geographic_location == "Netherlands") -> Netherlands  #timepoints are all zero in this study 
Netherlands %>% distinct(subject_id, .keep_all=TRUE) -> Netherlands #Remove duplicates

Netherlands %>% filter(Netherlands$disease == "Crohn's disease"| Netherlands$disease == "Control") -> Netherlands
Netherlands %>% select(sample_id, disease, Sex, Age) -> Netherlands_meta
column_to_rownames(Netherlands_meta, "sample_id") -> Netherlands_meta
Netherlands_meta$"disease"[which(Netherlands_meta$"disease" == "Crohn's disease", TRUE)] <- 1
Netherlands_meta$"disease"[which(Netherlands_meta$"disease" == "Control", TRUE)] <- 0
Netherlands_meta <- apply(Netherlands_meta,c(1:2),as.numeric) %>% as.data.frame

#SP
Netherlands_SP <- Netherlands[, -which (colnames(Netherlands) %in% c("subject_id", "environment_material", "timepoint", 
                                                                     "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                                                     "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(Netherlands_SP, "sample_id") -> Netherlands_SP

keep <- apply(Netherlands_SP, 2, mean) > X & apply(Netherlands_SP > 0, 2, sum) / nrow(Netherlands_SP) > Y # mean abundance > X、Prevalence > Y
Netherlands_SP <- Netherlands_SP[, keep]

#Remove non-human data
as.data.frame(t(Netherlands_SP)) -> Netherlands_SP
rownames_to_column(Netherlands_SP, var = "feature") -> Netherlands_SP
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(Netherlands_SP, Human, by = "feature") -> Netherlands_SP
Netherlands_SP %>% filter(Netherlands_SP$human == 1) -> Netherlands_SP # include only fungi or protozoa reported in human 
column_to_rownames(Netherlands_SP, var="feature") -> Netherlands_SP
Netherlands_SP <- Netherlands_SP[, -which (colnames(Netherlands_SP) %in% c("human"))]
as.data.frame(t(Netherlands_SP)) -> Netherlands_SP

fit_data = Maaslin2(
  input_data = Netherlands_SP, 
  input_metadata = Netherlands_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "Netherlands_CD", 
  normalization = "NONE",
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE, 
  standardize = FALSE,
  cores = 4,
  fixed_effects=c("disease,Age"),
  reference = c("disease,0")) 

fit_data$results -> Netherlands_maaslin
Netherlands_maaslin %>% filter(Netherlands_maaslin$"metadata" == "disease") -> Netherlands_maaslin_CD
Netherlands_maaslin_CD <- Netherlands_maaslin_CD[, -which (colnames(Netherlands_maaslin_CD) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=Netherlands_maaslin_CD, "Coefficient (CD Netherlands)" = "coef") -> Netherlands_maaslin_CD
rename(.data=Netherlands_maaslin_CD, "P-value (CD Netherlands)" = "pval") -> Netherlands_maaslin_CD
rename(.data=Netherlands_maaslin_CD, "Q-value (CD Netherlands)" = "qval") -> Netherlands_maaslin_CD

#IBDvsControl 
WO3 %>% filter(WO3$study == "Franzosa_2018_IBD" & WO3$geographic_location == "Netherlands") -> Netherlands #timepoints are all zero in this study 
Netherlands %>% distinct(subject_id, .keep_all=TRUE) -> Netherlands #Remove duplicates

Netherlands %>% filter(Netherlands$disease == "Ulcerative colitis"| Netherlands$disease == "Control" | Netherlands$disease == "Crohn's disease") -> Netherlands
Netherlands %>% select(sample_id, disease, Sex, Age) -> Netherlands_meta
column_to_rownames(Netherlands_meta, "sample_id") -> Netherlands_meta
Netherlands_meta$"disease"[which(Netherlands_meta$"disease" == "Ulcerative colitis" | Netherlands_meta$"disease" == "Crohn's disease", TRUE)] <- 1
Netherlands_meta$"disease"[which(Netherlands_meta$"disease" == "Control", TRUE)] <- 0
Netherlands_meta <- apply(Netherlands_meta,c(1:2),as.numeric) %>% as.data.frame

#SP
Netherlands_SP <- Netherlands[, -which (colnames(Netherlands) %in% c("subject_id", "environment_material", "timepoint", 
                                                                     "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                                                     "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(Netherlands_SP, "sample_id") -> Netherlands_SP

keep <- apply(Netherlands_SP, 2, mean) > X & apply(Netherlands_SP > 0, 2, sum) / nrow(Netherlands_SP) > Y # mean abundance > X、Prevalence > Y
Netherlands_SP <- Netherlands_SP[, keep]

#Remove non-human data
as.data.frame(t(Netherlands_SP)) -> Netherlands_SP
rownames_to_column(Netherlands_SP, var = "feature") -> Netherlands_SP
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(Netherlands_SP, Human, by = "feature") -> Netherlands_SP
Netherlands_SP %>% filter(Netherlands_SP$human == 1) -> Netherlands_SP # include only fungi or protozoa reported in human 
column_to_rownames(Netherlands_SP, var="feature") -> Netherlands_SP
Netherlands_SP <- Netherlands_SP[, -which (colnames(Netherlands_SP) %in% c("human"))]
as.data.frame(t(Netherlands_SP)) -> Netherlands_SP

fit_data = Maaslin2(
  input_data = Netherlands_SP, 
  input_metadata = Netherlands_meta, 
  min_abundance = 0, min_prevalence = 0, 
  output = "Netherlands_IBD", 
  normalization = "NONE",
  transform = "LOG",
  plot_scatter = FALSE, 
  plot_heatmap = FALSE,
  standardize = FALSE, 
  cores = 4,
  fixed_effects=c("disease,Age"),
  reference = c("disease,0")) 

fit_data$results -> Netherlands_maaslin
Netherlands_maaslin %>% filter(Netherlands_maaslin$"metadata" == "disease") -> Netherlands_maaslin_IBD
Netherlands_maaslin_IBD <- Netherlands_maaslin_IBD[, -which (colnames(Netherlands_maaslin_IBD) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=Netherlands_maaslin_IBD, "Coefficient (IBD Netherlands)" = "coef") -> Netherlands_maaslin_IBD
rename(.data=Netherlands_maaslin_IBD, "P-value (IBD Netherlands)" = "pval") -> Netherlands_maaslin_IBD
rename(.data=Netherlands_maaslin_IBD, "Q-value (IBD Netherlands)" = "qval") -> Netherlands_maaslin_IBD

#4D Japanese MaAsLin Fungi解析
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Fungi")
#UCvsControls 
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% filter(Meta$UC_Dx == 1| Meta$UC_Dx ==0) -> Meta
Meta %>% select(ID, UC_Dx, Age, sex, BMI) -> JP_meta

SP_input <- read.csv("feces_230925.eukaryote.csv", header = TRUE, na.strings = c(NA, ''), row.names=1, check.names = FALSE) 
as.data.frame(t(SP_input)) -> SP_input
rownames_to_column(SP_input, var = "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(JP_meta, SP_input, by ="ID") -> SP_input
column_to_rownames(SP_input, var="ID") -> SP_input
rownames(JP_meta) <- NULL 
column_to_rownames(JP_meta, var="ID") -> JP_meta

SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("UC_Dx", "Age", "sex", "BMI"))]
keep <- apply(SP_input2, 2, mean) > X & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > Y # mean abundance > X、Prevalence > Y
SP_input2  <- SP_input2[, keep]

#Remove non-human data
as.data.frame(t(SP_input2)) -> SP_input2
rownames_to_column(SP_input2, var = "feature") -> SP_input2
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(SP_input2, Human, by = "feature") -> SP_input2
SP_input2 %>% filter(SP_input2$human == 1) -> SP_input2 # include only fungi or protozoa reported in human 
column_to_rownames(SP_input2, var="feature") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("human"))]
as.data.frame(t(SP_input2)) -> SP_input2

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

fit_data$results -> JP_maaslin_SP
JP_maaslin_SP %>% filter(JP_maaslin_SP$"metadata" == "UC_Dx") -> JP_maaslin_UC_SP
JP_maaslin_UC_SP <- JP_maaslin_UC_SP[, -which (colnames(JP_maaslin_UC_SP) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))]
rename(.data=JP_maaslin_UC_SP, "Coefficient (UC JP)" = "coef") -> JP_maaslin_UC_SP
rename(.data=JP_maaslin_UC_SP, "P-value (UC JP)" = "pval") -> JP_maaslin_UC_SP
rename(.data=JP_maaslin_UC_SP, "Q-value (UC JP)" = "qval") -> JP_maaslin_UC_SP

#CDvsControl
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% filter(Meta$CD_Dx == 1| Meta$CD_Dx ==0) -> Meta
Meta %>% select(ID, CD_Dx, Age, sex, BMI) -> JP_meta

SP_input <- read.csv("feces_230925.eukaryote.csv", header = TRUE, na.strings = c(NA, ''), row.names=1, check.names = FALSE) #abundance data
as.data.frame(t(SP_input)) -> SP_input
rownames_to_column(SP_input, var = "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(JP_meta, SP_input, by ="ID") -> SP_input
column_to_rownames(SP_input, var="ID") -> SP_input
rownames(JP_meta) <- NULL 
column_to_rownames(JP_meta, var="ID") -> JP_meta

SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("CD_Dx", "Age", "sex", "BMI"))]
keep <- apply(SP_input2, 2, mean) > X & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > Y # mean abundance > X、Prevalence > Y
SP_input2  <- SP_input2[, keep]

#Remove non-human data
as.data.frame(t(SP_input2)) -> SP_input2
rownames_to_column(SP_input2, var = "feature") -> SP_input2
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(SP_input2, Human, by = "feature") -> SP_input2
SP_input2 %>% filter(SP_input2$human == 1) -> SP_input2 # include only fungi or protozoa reported in human 
column_to_rownames(SP_input2, var="feature") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("human"))]
as.data.frame(t(SP_input2)) -> SP_input2

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

fit_data$results -> JP_maaslin_SP
JP_maaslin_SP %>% filter(JP_maaslin_SP$"metadata" == "CD_Dx") -> JP_maaslin_CD_SP
JP_maaslin_CD_SP <- JP_maaslin_CD_SP[, -which (colnames(JP_maaslin_CD_SP) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))]
rename(.data=JP_maaslin_CD_SP, "Coefficient (CD JP)" = "coef") -> JP_maaslin_CD_SP
rename(.data=JP_maaslin_CD_SP, "P-value (CD JP)" = "pval") -> JP_maaslin_CD_SP
rename(.data=JP_maaslin_CD_SP, "Q-value (CD JP)" = "qval") -> JP_maaslin_CD_SP

#IBDvsControl
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% filter(Meta$IBD == 1| Meta$IBD ==0) -> Meta
Meta %>% select(ID, IBD, Age, sex, BMI) -> JP_meta

SP_input <- read.csv("feces_230925.eukaryote.csv", header = TRUE, na.strings = c(NA, ''), row.names=1, check.names = FALSE) #abundance data
as.data.frame(t(SP_input)) -> SP_input
rownames_to_column(SP_input, var = "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(JP_meta, SP_input, by ="ID") -> SP_input
column_to_rownames(SP_input, var="ID") -> SP_input
rownames(JP_meta) <- NULL 
column_to_rownames(JP_meta, var="ID") -> JP_meta

SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("IBD", "Age", "sex", "BMI"))]
keep <- apply(SP_input2, 2, mean) > X & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > Y # mean abundance > X、Prevalence > Y
SP_input2  <- SP_input2[, keep]

#Remove non-human data
as.data.frame(t(SP_input2)) -> SP_input2
rownames_to_column(SP_input2, var = "feature") -> SP_input2
Human<- read.csv("Human.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #This list includes fungi or protozoa reported in human stool.
rownames_to_column(Human, var = "feature") -> Human
left_join(SP_input2, Human, by = "feature") -> SP_input2
SP_input2 %>% filter(SP_input2$human == 1) -> SP_input2 # include only fungi or protozoa reported in human 
column_to_rownames(SP_input2, var="feature") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("human"))]
as.data.frame(t(SP_input2)) -> SP_input2

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

fit_data$results -> JP_maaslin_SP
JP_maaslin_SP %>% filter(JP_maaslin_SP$"metadata" == "IBD") -> JP_maaslin_IBD_SP
JP_maaslin_IBD_SP <- JP_maaslin_IBD_SP[, -which (colnames(JP_maaslin_IBD_SP) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))]
rename(.data=JP_maaslin_IBD_SP, "Coefficient (IBD JP)" = "coef") -> JP_maaslin_IBD_SP
rename(.data=JP_maaslin_IBD_SP, "P-value (IBD JP)" = "pval") -> JP_maaslin_IBD_SP
rename(.data=JP_maaslin_IBD_SP, "Q-value (IBD JP)" = "qval") -> JP_maaslin_IBD_SP

#Data integration
US_MA <- full_join(full_join(US_maaslin_IBD, US_maaslin_UC, by = "feature"), US_maaslin_CD, by = "feature")
US2_MA <- full_join(full_join(US2_maaslin_IBD, US2_maaslin_UC, by = "feature"), US2_maaslin_CD, by = "feature")
US_MA2 <-full_join(US_MA, US2_MA, by = "feature")

Spain_MA <- full_join(full_join(Spain_maaslin_IBD, Spain_maaslin_UC, by = "feature"), Spain_maaslin_CD, by = "feature")
Netherlands_MA <- full_join(full_join(Netherlands_maaslin_IBD, Netherlands_maaslin_UC, by = "feature"), Netherlands_maaslin_CD, by = "feature")
EU_MA <- full_join(Spain_MA, Netherlands_MA, by = "feature")

All_OTHER <- full_join(full_join(full_join(US_MA2, EU_MA, by = "feature"), Tanzania_maaslin_UC, by = "feature"), China_maaslin_CD, by = "feature")
JP_MA <- full_join(full_join(JP_maaslin_IBD_SP, JP_maaslin_UC_SP, by = "feature"), JP_maaslin_CD_SP, by = "feature")

All_country <- full_join(JP_MA, All_OTHER, by ="feature")

#For network analysis
JP_MA %>% filter(JP_MA$"P-value (UC JP)"< 0.05 | JP_MA$"P-value (CD JP)" < 0.05) -> Fungi_map
Fungi_map %>% arrange(-Fungi_map$"Coefficient (IBD JP)") ->Fungi_map

#Heatmap creation
#Select species with P-value<0.2 in UC or CD 
All_country %>% filter(All_country$"P-value (CD JP)" <0.2 | All_country$"P-value (UC JP)"<0.2) -> SP_input
SP_input %>% arrange(-SP_input$"Coefficient (IBD JP)") -> SP_input2 

SP_input2$feature <- str_replace_all(SP_input2$feature, pattern="[.]", replacement=" ") #Remove dot in the fungal name
column_to_rownames(SP_input2, "feature") ->SP_input2

SP_input2 %>% select("Coefficient (IBD JP)", "P-value (IBD JP)", "Q-value (IBD JP)", "Coefficient (UC JP)", "P-value (UC JP)", "Q-value (UC JP)", "Coefficient (CD JP)", "P-value (CD JP)", "Q-value (CD JP)", 
                                "Coefficient (IBD US)", "P-value (IBD US)", "Q-value (IBD US)", "Coefficient (UC US)", "P-value (UC US)", "Q-value (UC US)", "Coefficient (CD US)", "P-value (CD US)", "Q-value (CD US)",
                                "Coefficient (IBD US2)", "P-value (IBD US2)", "Q-value (IBD US2)", "Coefficient (UC US2)", "P-value (UC US2)", "Q-value (UC US2)", "Coefficient (CD US2)", "P-value (CD US2)", "Q-value (CD US2)", 
                                "Coefficient (IBD Spain)", "P-value (IBD Spain)", "Q-value (IBD Spain)", "Coefficient (UC Spain)", "P-value (UC Spain)", "Q-value (UC Spain)", "Coefficient (CD Spain)", "P-value (CD Spain)", "Q-value (CD Spain)", 
                                "Coefficient (IBD Netherlands)", "P-value (IBD Netherlands)", "Q-value (IBD Netherlands)", "Coefficient (UC Netherlands)", "P-value (UC Netherlands)", "Q-value (UC Netherlands)", "Coefficient (CD Netherlands)", "P-value (CD Netherlands)", "Q-value (CD Netherlands)", 
                                "Coefficient (UC Tanzania)", "P-value (UC Tanzania)", "Q-value (UC Tanzania)", "Coefficient (CD China)", "P-value (CD China)", "Q-value (CD China)") ->  Fig1 

Fig1 %>% select("P-value (IBD JP)", "P-value (IBD US)", "P-value (IBD US2)", "P-value (IBD Spain)", "P-value (IBD Netherlands)") -> pval_IBD
Fig1 %>% select("P-value (UC JP)", "P-value (UC US)",  "P-value (UC US2)", "P-value (UC Spain)", "P-value (UC Netherlands)", "P-value (UC Tanzania)") -> pval_UC
Fig1 %>% select("P-value (CD JP)", "P-value (CD US)", "P-value (CD US2)", "P-value (CD Spain)", "P-value (CD Netherlands)", "P-value (CD China)") -> pval_CD

pval_IBD[is.na(pval_IBD)] <- 1
pval_UC[is.na(pval_UC)] <- 1
pval_CD[is.na(pval_CD)] <- 1

Fig1 %>% select("Coefficient (IBD JP)", "Coefficient (IBD US)", "Coefficient (IBD US2)", "Coefficient (IBD Spain)", "Coefficient (IBD Netherlands)") -> Gram_coef_IBD
Fig1 %>% select("Coefficient (UC JP)", "Coefficient (UC US)",  "Coefficient (UC US2)", "Coefficient (UC Spain)", "Coefficient (UC Netherlands)", "Coefficient (UC Tanzania)") -> Gram_coef_UC
Fig1 %>% select("Coefficient (CD JP)", "Coefficient (CD US)", "Coefficient (CD US2)", "Coefficient (CD Spain)", "Coefficient (CD Netherlands)", "Coefficient (CD China)") -> Gram_coef_CD

anno_width = unit(2, "cm")

rename(.data= Gram_coef_IBD, "IBD (Japan)" = "Coefficient (IBD JP)") -> Gram_coef_IBD
rename(.data= Gram_coef_IBD, "IBD (United States, Lloyd-Price_2019)" = "Coefficient (IBD US)") -> Gram_coef_IBD
rename(.data= Gram_coef_IBD, "IBD (United States, Franzosa_2018)" = "Coefficient (IBD US2)") -> Gram_coef_IBD
rename(.data= Gram_coef_IBD, "IBD (Netherlands)" = "Coefficient (IBD Netherlands)") -> Gram_coef_IBD
rename(.data= Gram_coef_IBD, "IBD (Spain)" = "Coefficient (IBD Spain)") -> Gram_coef_IBD

rename(.data= Gram_coef_UC, "UC (Japan)" = "Coefficient (UC JP)") -> Gram_coef_UC
rename(.data= Gram_coef_UC, "UC (United States, Lloyd-Price_2019)" = "Coefficient (UC US)") -> Gram_coef_UC
rename(.data= Gram_coef_UC, "UC (United States, Franzosa_2018)" = "Coefficient (UC US2)") -> Gram_coef_UC
rename(.data= Gram_coef_UC, "UC (Netherlands)" = "Coefficient (UC Netherlands)") -> Gram_coef_UC
rename(.data= Gram_coef_UC, "UC (Spain)" = "Coefficient (UC Spain)") -> Gram_coef_UC
rename(.data= Gram_coef_UC, "UC (Tanzania)" = "Coefficient (UC Tanzania)") -> Gram_coef_UC

rename(.data= Gram_coef_CD, "CD (Japan)" = "Coefficient (CD JP)") -> Gram_coef_CD
rename(.data= Gram_coef_CD, "CD (United States, Lloyd-Price_2019)" = "Coefficient (CD US)") -> Gram_coef_CD
rename(.data= Gram_coef_CD, "CD (United States, Franzosa_2018)" = "Coefficient (CD US2)") -> Gram_coef_CD
rename(.data= Gram_coef_CD, "CD (Netherlands)" = "Coefficient (CD Netherlands)") -> Gram_coef_CD
rename(.data= Gram_coef_CD, "CD (Spain)" = "Coefficient (CD Spain)") -> Gram_coef_CD
rename(.data= Gram_coef_CD, "CD (China)" = "Coefficient (CD China)") -> Gram_coef_CD

lgd_sig = Legend(pch = "*", type = "points", labels = "P < 0.05")

p1=pheatmap(as.matrix(Gram_coef_IBD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(pval_IBD < 0.05,"*", ""), nrow(pval_IBD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.4,0,0.4), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))
p2=pheatmap(as.matrix(Gram_coef_UC), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(pval_UC < 0.05,"*", ""), nrow(pval_UC)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.4,0,0.4), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))
p3=pheatmap(as.matrix(Gram_coef_CD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(pval_CD < 0.05,"*", ""), nrow(pval_CD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.4,0,0.4), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))

draw(p1+p2+p3, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig))

#Removed Tanzania & Lloyd Price (True candidate for fig)
#Select species with P<0.2 in UC or CD 
All_country %>% filter(All_country$"P-value (CD JP)" < 0.2 | All_country$"P-value (UC JP)" < 0.2) -> SP_input
SP_input %>% arrange(-SP_input$"Coefficient (IBD JP)") -> SP_input2 
#column_to_rownames(SP_input2, "Feature") -> SP_input2_bacteriome

SP_input2$feature <- str_replace_all(SP_input2$feature, pattern="[.]", replacement=" ") #Remove dot in the fungal name
column_to_rownames(SP_input2, "feature") ->SP_input2

SP_input2 %>% select("Coefficient (IBD JP)", "P-value (IBD JP)", "Q-value (IBD JP)", "Coefficient (UC JP)", "P-value (UC JP)", "Q-value (UC JP)", "Coefficient (CD JP)", "P-value (CD JP)", "Q-value (CD JP)", 
                                "Coefficient (IBD US2)", "P-value (IBD US2)", "Q-value (IBD US2)", "Coefficient (UC US2)", "P-value (UC US2)",  "Q-value (UC US2)", "Coefficient (CD US2)", "P-value (CD US2)", "Q-value (CD US2)", 
                                "Coefficient (IBD Spain)", "P-value (IBD Spain)", "Q-value (IBD Spain)", "Coefficient (UC Spain)", "P-value (UC Spain)", "Q-value (UC Spain)", "Coefficient (CD Spain)", "P-value (CD Spain)", "Q-value (CD Spain)", 
                                "Coefficient (IBD Netherlands)", "P-value (IBD Netherlands)", "Q-value (IBD Netherlands)", "Coefficient (UC Netherlands)", "P-value (UC Netherlands)", "Q-value (UC Netherlands)", "Coefficient (CD Netherlands)", "P-value (CD Netherlands)", "Q-value (CD Netherlands)",  
                                "Coefficient (CD China)", "P-value (CD China)", "Q-value (CD China)") ->  Fig1 

Fig1 %>% select("P-value (IBD JP)", "P-value (IBD US2)", "P-value (IBD Spain)", "P-value (IBD Netherlands)") -> pval_IBD
Fig1 %>% select("P-value (UC JP)",  "P-value (UC US2)", "P-value (UC Spain)", "P-value (UC Netherlands)") -> pval_UC
Fig1 %>% select("P-value (CD JP)", "P-value (CD US2)", "P-value (CD Spain)", "P-value (CD Netherlands)", "P-value (CD China)") -> pval_CD

pval_IBD[is.na(pval_IBD)] <- 1
pval_UC[is.na(pval_UC)] <- 1
pval_CD[is.na(pval_CD)] <- 1

Fig1 %>% select("Coefficient (IBD JP)", "Coefficient (IBD US2)", "Coefficient (IBD Spain)", "Coefficient (IBD Netherlands)") -> Gram_coef_IBD
Fig1 %>% select("Coefficient (UC JP)",  "Coefficient (UC US2)", "Coefficient (UC Spain)", "Coefficient (UC Netherlands)") -> Gram_coef_UC
Fig1 %>% select("Coefficient (CD JP)", "Coefficient (CD US2)", "Coefficient (CD Spain)", "Coefficient (CD Netherlands)", "Coefficient (CD China)") -> Gram_coef_CD

anno_width = unit(2, "cm")

rename(.data= Gram_coef_IBD, "IBD (Japan)" = "Coefficient (IBD JP)") -> Gram_coef_IBD
rename(.data= Gram_coef_IBD, "IBD (United States)" = "Coefficient (IBD US2)") -> Gram_coef_IBD #Franzosa_2018
rename(.data= Gram_coef_IBD, "IBD (Netherlands)" = "Coefficient (IBD Netherlands)") -> Gram_coef_IBD
rename(.data= Gram_coef_IBD, "IBD (Spain)" = "Coefficient (IBD Spain)") -> Gram_coef_IBD

rename(.data= Gram_coef_UC, "UC (Japan)" = "Coefficient (UC JP)") -> Gram_coef_UC
rename(.data= Gram_coef_UC, "UC (United States)" = "Coefficient (UC US2)") -> Gram_coef_UC #Franzosa_2018
rename(.data= Gram_coef_UC, "UC (Netherlands)" = "Coefficient (UC Netherlands)") -> Gram_coef_UC
rename(.data= Gram_coef_UC, "UC (Spain)" = "Coefficient (UC Spain)") -> Gram_coef_UC

rename(.data= Gram_coef_CD, "CD (Japan)" = "Coefficient (CD JP)") -> Gram_coef_CD
rename(.data= Gram_coef_CD, "CD (United States)" = "Coefficient (CD US2)") -> Gram_coef_CD #Franzosa_2018
rename(.data= Gram_coef_CD, "CD (Netherlands)" = "Coefficient (CD Netherlands)") -> Gram_coef_CD
rename(.data= Gram_coef_CD, "CD (Spain)" = "Coefficient (CD Spain)") -> Gram_coef_CD
rename(.data= Gram_coef_CD, "CD (China)" = "Coefficient (CD China)") -> Gram_coef_CD

lgd_sig = Legend(pch = "*", type = "points", labels = "P < 0.05")

p1=pheatmap(as.matrix(Gram_coef_IBD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(pval_IBD < 0.05,"*", ""), nrow(pval_IBD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.4,0,0.4), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))
p2=pheatmap(as.matrix(Gram_coef_UC), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(pval_UC < 0.05,"*", ""), nrow(pval_UC)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.4,0,0.4), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))
p3=pheatmap(as.matrix(Gram_coef_CD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(pval_CD < 0.05,"*", ""), nrow(pval_CD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.4,0,0.4), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))

draw(p1+p2+p3, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig))

#Spearman correlation analysis
#USvsJP 
sp <- ggplot(All_country, aes(x = All_country$"Coefficient (IBD JP)", y = All_country$"Coefficient (IBD US)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japan vs United States (Lloyd-Price_2019)", subtitle = "IBD", x = "Coefficient value (Japan)", y = "Coefficient value (US)", tag = "A") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_us_jp_ibd <- sp +  stat_cor(method = "spearman", label.x = -0.25, label.y = 1, cor.coef.name = c("rho"))

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (UC JP)", y = All_country$"Coefficient (UC US)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japan vs United States (Lloyd-Price_2019)", subtitle = "UC", x = "Coefficient value (Japan)", y = "Coefficient value (US)", tag = "B") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_us_jp_uc <- sp + stat_cor(method = "spearman", label.x = -0.25, label.y = 1, cor.coef.name = c("rho"))

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (CD JP)", y = All_country$"Coefficient (CD US)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japan vs United States (Lloyd-Price_2019)", subtitle = "CD", x = "Coefficient value (Japan)", y = "Coefficient value (US)", tag = "C") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_us_jp_cd <- sp + stat_cor(method = "spearman", label.x = -0.25, label.y = 1, cor.coef.name = c("rho"))

#US2vsJP 
sp <- ggplot(All_country, aes(x = All_country$"Coefficient (IBD JP)", y = All_country$"Coefficient (IBD US2)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japan vs United State (Franzosa_2018)", subtitle = "IBD", x = "Coefficient value (Japan)", y = "Coefficient value (US)", tag = "D") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_US2_jp_ibd <- sp +  stat_cor(method = "spearman", label.x = -0.25, label.y = 1, cor.coef.name = c("rho"))

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (UC JP)", y = All_country$"Coefficient (UC US2)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japan vs United State (Franzosa_2018)", subtitle = "UC", x = "Coefficient value (Japan)", y = "Coefficient value (US)", tag = "E") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_US2_jp_uc <- sp + stat_cor(method = "spearman", label.x = -0.25, label.y = 1, cor.coef.name = c("rho"))

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (CD JP)", y = All_country$"Coefficient (CD US2)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japan vs United State (Franzosa_2018)", subtitle = "CD", x = "Coefficient value (Japan)", y = "Coefficient value (US)", tag = "F") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_US2_jp_cd <- sp + stat_cor(method = "spearman", label.x = -0.25, label.y = 2, cor.coef.name = c("rho"))

combine_plots(
  list(sp_us_jp_ibd, sp_us_jp_uc, sp_us_jp_cd, sp_US2_jp_ibd, sp_US2_jp_uc, sp_US2_jp_cd),
  plotgrid.args = list(nrow = 2),
  annotation.args = list(
    title = "Fungi between Japan and United States",
    caption = ""
  )
)

#SpainvsJP
sp <- ggplot(All_country, aes(x = All_country$"Coefficient (IBD JP)", y = All_country$"Coefficient (IBD Spain)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japan vs Spain", subtitle = "IBD", x = "Coefficient value (Japan)", y = "Coefficient value (Spain)", tag = "a") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_Spain_jp_ibd <- sp + stat_cor(method = "spearman", label.x = -0.25, label.y = 1, cor.coef.name = c("rho"))

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (UC JP)", y = All_country$"Coefficient (UC Spain)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japan vs Spain", subtitle = "UC", x = "Coefficient value (Japan)", y = "Coefficient value (Spain)", tag = "b") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_Spain_jp_uc <- sp + stat_cor(method = "spearman", label.x = -0.25, label.y = 1, cor.coef.name = c("rho"))

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (CD JP)", y = All_country$"Coefficient (CD Spain)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japan vs Spain", subtitle = "CD", x = "Coefficient value (Japan)", y = "Coefficient value (Spain)", tag = "c") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_Spain_jp_cd <- sp + stat_cor(method = "spearman", label.x = -0.25, label.y = 1, cor.coef.name = c("rho"))

#TanzaniavsJP
sp <- ggplot(All_country, aes(x = All_country$"Coefficient (UC JP)", y = All_country$"Coefficient (UC Tanzania)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japan vs Tanzania", subtitle = "UC", x = "Coefficient value (Japan)", y = "Coefficient value (Tanzania)", tag = "A") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_tan_jp_uc <- sp + stat_cor(method = "spearman", label.x = -0.25, label.y = 1, cor.coef.name = c("rho"))

#NetherlandsvsJP 
sp <- ggplot(All_country, aes(x = All_country$"Coefficient (IBD JP)", y = All_country$"Coefficient (IBD Netherlands)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japan vs Netherlands", subtitle = "IBD", x = "Coefficient value (Japan)", y = "Coefficient value (Netherlands)", tag = "d") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_Netherlands_jp_ibd <- sp +  stat_cor(method = "spearman", label.x = -0.25, label.y = 1, cor.coef.name = c("rho"))

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (UC JP)", y = All_country$"Coefficient (UC Netherlands)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japan vs Netherlands", subtitle = "UC", x = "Coefficient value (Japan)", y = "Coefficient value (Netherlands)", tag = "e") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_Netherlands_jp_uc <- sp + stat_cor(method = "spearman", label.x = -0.25, label.y = 1, cor.coef.name = c("rho"))

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (CD JP)", y = All_country$"Coefficient (CD Netherlands)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japan vs Netherlands", subtitle = "CD", x = "Coefficient value (Japan)", y = "Coefficient value (Netherlands)", tag = "f") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_Netherlands_jp_cd <- sp + stat_cor(method = "spearman", label.x = -0.25, label.y = 1, cor.coef.name = c("rho"))

#ChinavsJP
sp <- ggplot(All_country, aes(x = All_country$"Coefficient (CD JP)", y = All_country$"Coefficient (CD China)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japan vs China", subtitle = "CD", x = "Coefficient value (Japan)", y = "Coefficient value (China)", tag = "B") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_china_jp_cd <- sp + stat_cor(method = "spearman", label.x = -0.25, label.y = 1, cor.coef.name = c("rho"))

combine_plots(
  list(sp_Spain_jp_ibd, sp_Spain_jp_uc, sp_Spain_jp_cd, sp_Netherlands_jp_ibd, sp_Netherlands_jp_uc, sp_Netherlands_jp_cd),
  plotgrid.args = list(nrow = 2),
  annotation.args = list(
    title = "Fungi between Japan and EU",
    caption = ""
  )
)

combine_plots(
  list(sp_tan_jp_uc, sp_china_jp_cd),
  plotgrid.args = list(nrow = 1),
  annotation.args = list(
    title = "Species between Japan and Other Countries (mOTU3)",
    caption = ""
  )
)

#True candidate for figure (Only Franzosa for US and Llyoid-Price and Tanzanian data are excluded)
sp <- ggplot(All_country, aes(x = All_country$"Coefficient (CD JP)", y = All_country$"Coefficient (CD China)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japan vs China", subtitle = "CD", x = "Coefficient value (Japan)", y = "Coefficient value (China)", tag = "g") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_china_jp_cd <- sp + stat_cor(method = "spearman", label.x = -0.25, label.y = 2, cor.coef.name = c("rho"))

combine_plots(
  list(sp_Spain_jp_ibd, sp_Spain_jp_uc, sp_Spain_jp_cd, sp_Netherlands_jp_ibd, sp_Netherlands_jp_uc, sp_Netherlands_jp_cd, sp_china_jp_cd),
  plotgrid.args = list(nrow = 3),
  annotation.args = list(
    title = "Species between Japan and EU & China (mOTU3)",
    caption = ""
  )
)

