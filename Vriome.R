library(dplyr)
library(tidyverse)
library(Maaslin2)
library(ComplexHeatmap)
library(ggpubr)
library(ggstatsplot)
library(reshape2)

#MaAsLin Virome analysis (World data, Figure 4e)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
SP_WO3<- read.csv("phage.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = FALSE) 
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3

Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")

#He_2017_Crohn (China)
WO3 %>% filter(WO3$study == "He_2017_Crohn" & WO3$timepoint == 0) -> China #timepoint = 0 abundance data (RPKM)
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

keep <- apply(China_SP, 2, mean) > 1E-4 & apply(China_SP > 0, 2, sum) / nrow(China_SP) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
China_SP <- China_SP[, keep]

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
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
SP_WO3<- read.csv("phage.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1) #check.names deleted as Spain ID has some unique ID modified by R
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3

Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")

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
Spain_SP <- na.omit(Spain_SP) #Spain data has NA for phages and fungi

keep <- apply(Spain_SP, 2, mean) > 1E-4 & apply(Spain_SP > 0, 2, sum) / nrow(Spain_SP) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
Spain_SP <- Spain_SP[, keep]

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
Spain_SP <- na.omit(Spain_SP)

keep <- apply(Spain_SP, 2, mean) > 1E-4 & apply(Spain_SP > 0, 2, sum) / nrow(Spain_SP) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
Spain_SP <- Spain_SP[, keep]

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
Spain %>% select(sample_id, disease, Sex, Age, bmi) -> Spain_meta
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
Spain_SP <- na.omit(Spain_SP)

keep <- apply(Spain_SP, 2, mean) > 1E-4 & apply(Spain_SP > 0, 2, sum) / nrow(Spain_SP) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
Spain_SP <- Spain_SP[, keep]

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

#Franzosa_2018_IBD (US cohort)
#UCvsControl
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
SP_WO3<- read.csv("phage.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = FALSE) #check.names FALSE again after Spain data
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3

Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")

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

keep <- apply(US2_SP, 2, mean) > 1E-4 & apply(US2_SP > 0, 2, sum) / nrow(US2_SP) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
US2_SP <- US2_SP[, keep]

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

keep <- apply(US2_SP, 2, mean) > 1E-4 & apply(US2_SP > 0, 2, sum) / nrow(US2_SP) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
US2_SP <- US2_SP[, keep]

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

keep <- apply(US2_SP, 2, mean) > 1E-4 & apply(US2_SP > 0, 2, sum) / nrow(US2_SP) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
US2_SP <- US2_SP[, keep]

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

keep <- apply(Netherlands_SP, 2, mean) > 1E-4 & apply(Netherlands_SP > 0, 2, sum) / nrow(Netherlands_SP) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
Netherlands_SP <- Netherlands_SP[, keep]

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

keep <- apply(Netherlands_SP, 2, mean) > 1E-4 & apply(Netherlands_SP > 0, 2, sum) / nrow(Netherlands_SP) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
Netherlands_SP <- Netherlands_SP[, keep]

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

keep <- apply(Netherlands_SP, 2, mean) > 1E-4 & apply(Netherlands_SP > 0, 2, sum) / nrow(Netherlands_SP) > 0.1 # mean abundance > 1E-4、Prevalence > 0.1
Netherlands_SP <- Netherlands_SP[, keep]

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

#4D JP cohort
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% filter(Meta$UC_Dx == 1| Meta$UC_Dx ==0) -> Meta

SP_input <- read.csv("vOTU.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(Meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("IBD_HC2", "Age", "sex", "BMI", "IBD_HC2.1", "UC", "Crohn", "UC_Dx", "CD_Dx", "IBD", "IBD_HC"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2

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

SP_input <- read.csv("vOTU.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(Meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("IBD_HC2", "Age", "sex", "BMI", "IBD_HC2.1", "UC", "Crohn", "UC_Dx", "CD_Dx", "IBD", "IBD_HC"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2

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

SP_input <- read.csv("vOTU.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(Meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("IBD_HC2", "Age", "sex", "BMI", "IBD_HC2.1", "UC", "Crohn", "UC_Dx", "CD_Dx", "IBD", "IBD_HC"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2

rownames(Meta) <- NULL　#行名のリセット
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

#Data integration
US2_MA <- full_join(full_join(US2_maaslin_IBD, US2_maaslin_UC, by = "feature"), US2_maaslin_CD, by = "feature")

Spain_MA <- full_join(full_join(Spain_maaslin_IBD, Spain_maaslin_UC, by = "feature"), Spain_maaslin_CD, by = "feature")
Netherlands_MA <- full_join(full_join(Netherlands_maaslin_IBD, Netherlands_maaslin_UC, by = "feature"), Netherlands_maaslin_CD, by = "feature")
EU_MA <- full_join(Spain_MA, Netherlands_MA, by = "feature")

All_OTHER <- full_join(full_join(US2_MA, EU_MA, by = "feature"), China_maaslin_CD, by = "feature")
JP_MA <- full_join(full_join(JP_maaslin_IBD, JP_maaslin_UC, by = "feature"), JP_maaslin_CD, by = "feature")

All_country <- full_join(JP_MA, All_OTHER, by ="feature")

#For network analysis
vOTU_ID <- read.csv("vOTU_ID.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
left_join(JP_MA, vOTU_ID, by = "feature") -> JP_MA
Ref <- read.csv("Reference.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
left_join(JP_MA, Ref, by = "vOTU") -> JP_MA
JP_MA %>% filter(abs(JP_MA$"Coefficient (CD JP)") > 1 & JP_MA$"Q-value (CD JP)" < 0.1 | abs(JP_MA$"Coefficient (UC JP)") > 1 & JP_MA$"Q-value (UC JP)" < 0.1) -> vOTU_Heat_qval_fig  
vOTU_Heat_qval_fig %>% arrange(-vOTU_Heat_qval_fig$"Coefficient (IBD JP)") -> vOTU_Heat_qval_fig 
column_to_rownames(vOTU_Heat_qval_fig, "feature2") -> vOTU_Heat_qval_fig

#Heatmap creation (Figure 4e)
vOTU_ID <- read.csv("vOTU_ID.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
left_join(All_country, vOTU_ID, by = "feature") -> All_country2 #add vOTU number

Ref <- read.csv("Reference.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
left_join(All_country2, Ref, by = "vOTU") -> All_country2 #add phage reference data

#Select vOTU with abs coef >1 and FDR <0.25 in UC or CD
All_country2 %>% filter(abs(All_country2$"Coefficient (CD JP)") > 1 & All_country2$"Q-value (CD JP)" < 0.25 | abs(All_country2$"Coefficient (UC JP)") > 1 & All_country2$"Q-value (UC JP)" < 0.25) -> SP_input  
SP_input %>% arrange(-SP_input$"Coefficient (IBD JP)") -> SP_input2 

column_to_rownames(SP_input2, "feature2") -> SP_input2 
SP_input2 %>% select(c(55:64)) -> species_table_joined #Extract predicted host data (need to check the location of host data)

SP_input2 %>% select(Host..phylum., Lytic_Lysogenic, Host) -> annotation
rename(.data= annotation, "Host (phylum)" = "Host..phylum.") -> annotation
rename(.data= annotation, "Lytic/Lysogenic" = "Lytic_Lysogenic") -> annotation
rename(.data= annotation, "Specialist/Generalist" = "Host") -> annotation
annotation$"Host (phylum)"[is.na(annotation$"Host (phylum)")] <- "Undefined"

rownames_to_column(species_table_joined, "vOTU") -> species_table_joined
species_table_joined %>% select(vOTU) -> Host_list
species_table_joined <- species_table_joined[, -which (colnames(species_table_joined) %in% c("vOTU"))]  #Remove vOTU

species_count <- species_table_joined  %>% as.matrix() %>% melt() %>% select(3) %>% table() %>% as.data.frame()
species_list <- as.character(species_count[order(species_count$Freq, decreasing = TRUE), 1])

annotation_species <- Host_list %>% as.data.frame()
for (i in 1:length(species_list)){
  for (j in 1:nrow(Host_list)){
    target <- species_table_joined[j, 1:10]
    search <- species_list[i]
    if (search %in% target){
      annotation_species[j, i+1] <- "infect"
    } else{
      annotation_species[j, i+1] <- "not infect"
    }
  }
}

colnames(annotation_species) <- c("feature", species_list)
col_species = c("infect" = "#F59B6A", "not infect" = "#E3EABD")

Legend_list <- list(`Ruminococcus` = list(title = "Predicted host (genus)", at = c("infect", "not infect"), labels = c("YES", "NO")))
annotation_species %>% select(-1) -> annotation_species

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")

SP_input2 %>% select("Coefficient (IBD JP)", "P-value (IBD JP)", "Q-value (IBD JP)", "Coefficient (UC JP)", "P-value (UC JP)", "Q-value (UC JP)", "Coefficient (CD JP)", "P-value (CD JP)", "Q-value (CD JP)", 
                                "Coefficient (IBD US2)", "P-value (IBD US2)", "Q-value (IBD US2)", "Coefficient (UC US2)", "P-value (UC US2)", "Q-value (UC US2)", "Coefficient (CD US2)", "P-value (CD US2)", "Q-value (CD US2)", 
                                "Coefficient (IBD Spain)", "P-value (IBD Spain)", "Q-value (IBD Spain)", "Coefficient (UC Spain)", "P-value (UC Spain)", "Q-value (UC Spain)", "Coefficient (CD Spain)", "P-value (CD Spain)", "Q-value (CD Spain)", 
                                "Coefficient (IBD Netherlands)", "P-value (IBD Netherlands)", "Q-value (IBD Netherlands)", "Coefficient (UC Netherlands)", "P-value (UC Netherlands)", "Q-value (UC Netherlands)", "Coefficient (CD Netherlands)", "P-value (CD Netherlands)", "Q-value (CD Netherlands)", 
                                "Coefficient (CD China)", "P-value (CD China)", "Q-value (CD China)") ->  Fig1 

Fig1 %>% select("Q-value (IBD JP)", "Q-value (IBD US2)", "Q-value (IBD Spain)", "Q-value (IBD Netherlands)") -> qval_IBD
Fig1 %>% select("Q-value (UC JP)",  "Q-value (UC US2)", "Q-value (UC Spain)", "Q-value (UC Netherlands)") -> qval_UC
Fig1 %>% select("Q-value (CD JP)", "Q-value (CD US2)", "Q-value (CD Spain)", "Q-value (CD Netherlands)", "Q-value (CD China)") -> qval_CD

qval_IBD[is.na(qval_IBD)] <- 1
qval_UC[is.na(qval_UC)] <- 1
qval_CD[is.na(qval_CD)] <- 1

Fig1 %>% select("Coefficient (IBD JP)", "Coefficient (IBD US2)", "Coefficient (IBD Spain)", "Coefficient (IBD Netherlands)") -> Gram_coef_IBD
Fig1 %>% select("Coefficient (UC JP)",  "Coefficient (UC US2)", "Coefficient (UC Spain)", "Coefficient (UC Netherlands)") -> Gram_coef_UC
Fig1 %>% select("Coefficient (CD JP)", "Coefficient (CD US2)", "Coefficient (CD Spain)", "Coefficient (CD Netherlands)", "Coefficient (CD China)") -> Gram_coef_CD

anno_width = unit(2, "cm")

rename(.data= Gram_coef_IBD, "Japanese 4D cohort" = "Coefficient (IBD JP)") -> Gram_coef_IBD
rename(.data= Gram_coef_IBD, "US cohort" = "Coefficient (IBD US2)") -> Gram_coef_IBD 
rename(.data= Gram_coef_IBD, "NL cohort" = "Coefficient (IBD Netherlands)") -> Gram_coef_IBD
rename(.data= Gram_coef_IBD, "ES cohort" = "Coefficient (IBD Spain)") -> Gram_coef_IBD

rename(.data= Gram_coef_UC, "Japanese 4D cohort" = "Coefficient (UC JP)") -> Gram_coef_UC
rename(.data= Gram_coef_UC, "US cohort" = "Coefficient (UC US2)") -> Gram_coef_UC 
rename(.data= Gram_coef_UC, "NL cohort" = "Coefficient (UC Netherlands)") -> Gram_coef_UC
rename(.data= Gram_coef_UC, "ES cohort" = "Coefficient (UC Spain)") -> Gram_coef_UC

rename(.data= Gram_coef_CD, "Japanese 4D cohort" = "Coefficient (CD JP)") -> Gram_coef_CD
rename(.data= Gram_coef_CD, "US cohort" = "Coefficient (CD US2)") -> Gram_coef_CD 
rename(.data= Gram_coef_CD, "NL cohort" = "Coefficient (CD Netherlands)") -> Gram_coef_CD
rename(.data= Gram_coef_CD, "ES cohort" = "Coefficient (CD Spain)") -> Gram_coef_CD
rename(.data= Gram_coef_CD, "CN cohort" = "Coefficient (CD China)") -> Gram_coef_CD

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")

ra3 = rowAnnotation(df=annotation_species, col = list('Ruminococcus' = col_species, `Bacteroides` = col_species,`Blautia` = col_species,`Streptococcus` = col_species,`Faecalibacterium` = col_species, `Clostridium` = col_species, `Eubacterium` = col_species, `Coprococcus` = col_species, `Roseburia` = col_species, `Parabacteroides` = col_species, `Fusicatenibacter` = col_species, `Porphyromonas` = col_species, `Unclassified Firmicutes` = col_species, `Lactobacillus` = col_species, `Anaerostipes` = col_species, `Dorea` = col_species, `Unclassified Lachnospiraceae` = col_species, `Alistipes` = col_species, `Bifidobacterium` = col_species, `Anaerotruncus` = col_species, `Bacillus` = col_species,`Prevotella` = col_species,`Akkermansia` = col_species,`Collinsella` = col_species,`Eggerthella` = col_species,`Escherichia` = col_species,`Flavonifractor` = col_species,`Gemmiger` = col_species,`Lachnoclostridium` = col_species, `Odoribacter` = col_species,`Petroclostridium` = col_species,`Tyzzerella` = col_species,`Unclassified Clostridia` = col_species,`Unclassified Clostridiales` = col_species,'Shigella'= col_species, 'Achromobacter'=col_species, 'Senegalimassilia'=col_species,'Veillonella'=col_species, 'Others'= col_species), gp = gpar(col = "black"), simple_anno_size = unit(2, "mm"),  annotation_name_gp= gpar(fontsize = 7),  show_legend = c(FALSE, TRUE,FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,FALSE, FALSE, FALSE, FALSE, FALSE,FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), annotation_legend_param = Legend_list)　

ra_spe = rowAnnotation("Lytic/Lysogenic" = annotation$"Lytic/Lysogenic","Specialist/Generalist" = annotation$"Specialist/Generalist", 
                       col = list("Lytic/Lysogenic"=c("Lysogenic"="#d8f255", "Lytic" = "purple"), "Specialist/Generalist"=c("Specialist" = "#ffcabf", "Generalist" = "#ffff80", "Host unknown" = "#c8c8cb")), 
                       simple_anno_size = unit(2, "mm"), gp = gpar(col = "black"), annotation_name_gp= gpar(fontsize = 7))

ra_host = rowAnnotation("Host (phylum)"=annotation$"Host (phylum)", col=list("Host (phylum)"=c("Actinobacteria"="#03af7a", "Bacteroidetes" = "purple", "Bacteroidetes, Verrucomicrobia" = "#f6aa00", "Firmicutes" = "#ff8082", "Firmicutes, Bacteroidetes" = "turquoise3", "Proteobacteria" = "#804000", "Undefined" = "#c8c8cb")), simple_anno_size = unit(2, "mm"), gp = gpar(col = "black"), annotation_name_gp= gpar(fontsize = 7))

p1=pheatmap(as.matrix(Gram_coef_IBD), fontsize = 7, cellwidth = 6, cellheight = 6, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 7, fontsize_row = 7, display_numbers = matrix(ifelse(qval_IBD < 0.1,"*", ""), nrow(qval_IBD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-2.5,0,1.5), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))
p2=pheatmap(as.matrix(Gram_coef_UC), fontsize = 7, cellwidth = 6, cellheight = 6, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 7, fontsize_row = 7, display_numbers = matrix(ifelse(qval_UC < 0.1,"*", ""), nrow(qval_UC)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-2.5,0,1.5), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))
p3=pheatmap(as.matrix(Gram_coef_CD), fontsize = 7, cellwidth = 6, cellheight = 6, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 7, fontsize_row = 7, display_numbers = matrix(ifelse(qval_CD < 0.1,"*", ""), nrow(qval_CD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-2.5,0,1.5), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))

draw(ra_host+ra3+ra_spe+p1+p2+p3, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig))

#Spearman correlation analysis (Supplementary Figure 10)
#US2(Franzosa)vsJP 
sp <- ggplot(All_country, aes(x = All_country$"Coefficient (IBD JP)", y = All_country$"Coefficient (IBD US2)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japanese 4D cohort vs US cohort", subtitle = "IBD", x = "Coefficient value (Japanese 4D cohort)", y = "Coefficient value (US cohort)", tag = "a") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_US2_jp_ibd <- sp +  stat_cor(method = "spearman", label.x = -2, label.y = 2, cor.coef.name = c("rho"))

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (UC JP)", y = All_country$"Coefficient (UC US2)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japanese 4D cohort vs US cohort", subtitle = "UC", x = "Coefficient value (Japanese 4D cohort)", y = "Coefficient value (US cohort)", tag = "b") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_US2_jp_uc <- sp + stat_cor(method = "spearman", label.x = -2, label.y = 2, cor.coef.name = c("rho"))

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (CD JP)", y = All_country$"Coefficient (CD US2)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japanese 4D cohort vs US cohort", subtitle = "CD", x = "Coefficient value (Japanese 4D cohort)", y = "Coefficient value (US cohort)", tag = "c") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_US2_jp_cd <- sp + stat_cor(method = "spearman", label.x = -3, label.y = 2.5, cor.coef.name = c("rho"))

#SpainvsJP
sp <- ggplot(All_country, aes(x = All_country$"Coefficient (IBD JP)", y = All_country$"Coefficient (IBD Spain)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japanese 4D cohort vs ES cohort", subtitle = "IBD", x = "Coefficient value (Japanese 4D cohort)", y = "Coefficient value (ES cohort)", tag = "d") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_Spain_jp_ibd <- sp + stat_cor(method = "spearman", label.x = -2, label.y = 2, cor.coef.name = c("rho"))

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (UC JP)", y = All_country$"Coefficient (UC Spain)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japanese 4D cohort vs ES cohort", subtitle = "UC", x = "Coefficient value (Japanese 4D cohort)", y = "Coefficient value (ES cohort)", tag = "e") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_Spain_jp_uc <- sp + stat_cor(method = "spearman", label.x = -2, label.y = 2, cor.coef.name = c("rho"))

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (CD JP)", y = All_country$"Coefficient (CD Spain)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japanese 4D cohort vs ES cohort", subtitle = "CD", x = "Coefficient value (Japanese 4D cohort)", y = "Coefficient value (ES cohort)", tag = "f") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_Spain_jp_cd <- sp + stat_cor(method = "spearman", label.x = -2, label.y = 4, cor.coef.name = c("rho"))

#NetherlandsvsJP 
sp <- ggplot(All_country, aes(x = All_country$"Coefficient (IBD JP)", y = All_country$"Coefficient (IBD Netherlands)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japanese 4D cohort vs NL cohort", subtitle = "IBD", x = "Coefficient value (Japanese 4D cohort)", y = "Coefficient value (NL cohort)", tag = "g") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_Netherlands_jp_ibd <- sp +  stat_cor(method = "spearman", label.x = -2, label.y = 2, cor.coef.name = c("rho"))

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (UC JP)", y = All_country$"Coefficient (UC Netherlands)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japanese 4D cohort vs NL cohort", subtitle = "UC", x = "Coefficient value (Japanese 4D cohort)", y = "Coefficient value (NL cohort)", tag = "h") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_Netherlands_jp_uc <- sp + stat_cor(method = "spearman", label.x = -2, label.y = 2, cor.coef.name = c("rho"))

sp <- ggplot(All_country, aes(x = All_country$"Coefficient (CD JP)", y = All_country$"Coefficient (CD Netherlands)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japanese 4D cohort vs NL cohort", subtitle = "CD", x = "Coefficient value (Japanese 4D cohort)", y = "Coefficient value (NL cohort)", tag = "i") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_Netherlands_jp_cd <- sp + stat_cor(method = "spearman", label.x = -2.5, label.y = 2.5, cor.coef.name = c("rho"))

#ChinavsJP
sp <- ggplot(All_country, aes(x = All_country$"Coefficient (CD JP)", y = All_country$"Coefficient (CD China)")) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(title = "Japanese 4D cohort vs CN cohort", subtitle = "CD", x = "Coefficient value (Japanese 4D cohort)", y = "Coefficient value (CN cohort)", tag = "j") +
  theme(plot.title = element_text(face = "bold", color = "black"), plot.tag  = element_text(face = "bold", color = "black", size =24))
sp_china_jp_cd <- sp + stat_cor(method = "spearman", label.x = -3, label.y = 3, cor.coef.name = c("rho"))

combine_plots(
  list(sp_US2_jp_ibd, sp_US2_jp_uc, sp_US2_jp_cd, sp_Spain_jp_ibd, sp_Spain_jp_uc, sp_Spain_jp_cd, sp_Netherlands_jp_ibd, sp_Netherlands_jp_uc, sp_Netherlands_jp_cd, sp_china_jp_cd),
  plotgrid.args = list(nrow = 4),
  annotation.args = list(
    title = "Phages between Japanese 4D cohort and external cohorts",
    caption = ""
  )
)
