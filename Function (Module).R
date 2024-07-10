#Function (Module) analysis of bacteriome
library(dplyr)
library(tidyverse)
library(Maaslin2)
library(ggplot2)
library(ggrepel)

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Module")
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% filter(Meta$UC_Dx == 1| Meta$UC_Dx ==0) -> Meta

SP_input <- read.csv("MO.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(Meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("IBD_HC2", "Age", "sex", "BMI", "IBD_HC2.1", "UC", "Crohn", "UC_Dx", "CD_Dx", "IBD"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2

rownames(Meta) <- NULL　
column_to_rownames(Meta, "ID") -> Meta
Meta %>% dplyr::select(UC_Dx, Age, sex, BMI) -> JP_meta

keep <- apply(SP_input2, 2, mean) > 1E-6 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-6, prevalence > 0.1
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

SP_input <- read.csv("MO.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(Meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("IBD_HC2", "Age", "sex", "BMI", "IBD_HC2.1", "UC", "Crohn", "UC_Dx", "CD_Dx", "IBD"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2

rownames(Meta) <- NULL　
column_to_rownames(Meta, "ID") -> Meta
Meta %>% dplyr::select(CD_Dx, Age, sex, BMI) -> JP_meta

keep <- apply(SP_input2, 2, mean) > 1E-6 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-6, prevalence > 0.1
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

SP_input <- read.csv("MO.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(Meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("IBD_HC2", "Age", "sex", "BMI", "IBD_HC2.1", "UC", "Crohn", "UC_Dx", "CD_Dx", "IBD"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2

rownames(Meta) <- NULL　
column_to_rownames(Meta, "ID") -> Meta
Meta %>% dplyr::select(IBD, Age, sex, BMI) -> JP_meta

keep <- apply(SP_input2, 2, mean) > 1E-6 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-6, prevalence > 0.1
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

JP_MA <- full_join(full_join(JP_maaslin_IBD, JP_maaslin_UC, by = "feature"), JP_maaslin_CD, by = "feature")

Feature<- read.csv("Feature.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
rownames_to_column(Feature, "feature") -> Feature
inner_join(JP_MA, Feature, by = "feature") -> volcano

#Volcano plot creation for UC
volcano %>% mutate(log10q_IBD = -log10(volcano$"Q-value (IBD JP)"), log10q_UC = -log10(volcano$"Q-value (UC JP)"), log10q_CD = -log10(volcano$"Q-value (CD JP)"))-> volcano_q
Up_Hits = head(arrange(volcano_q, -volcano_q$"Coefficient (UC JP)"),10)
Down_Hits = head(arrange(volcano_q, volcano_q$"Coefficient (UC JP)"),10)
Hits=rbind(Up_Hits, Down_Hits)  
volcano_q$label = if_else(volcano_q$feature %in% Hits$feature, volcano_q$feature, "") #Coef included in Top 10 and Bottom 10 are labeled.

X <- 0.9
Y <- 6 

volcano_q %>% mutate(class = case_when(volcano_q$"Coefficient (UC JP)" >= X & log10q_UC >= 1 ~ "Increase", volcano_q$"Coefficient (UC JP)" <= -X & log10q_UC >= 1 ~ "Decrease", TRUE ~ "None")) -> volcano_q 
volcano_q %>% mutate(class2 = case_when(volcano_q$"Coefficient (UC JP)" >= X & log10q_UC >= Y ~ "Increase", volcano_q$"Coefficient (UC JP)" <= -X & log10q_UC >= Y ~ "Decrease", TRUE ~ "None")) -> volcano_q
volcano_q %>% filter(class2 == "Increase" | class2 == "Decrease") -> Hit2
volcano_q$label = if_else(volcano_q$ID %in% Hit2$ID, volcano_q$ID, "")

ggplot(volcano_q, aes(x = `Coefficient (UC JP)`, y = log10q_UC, color=class, size=log10q_UC, label = rownames(ID))) +
  geom_vline(xintercept = c(-X, X), linetype =2, color="Grey") +
  geom_hline(yintercept = 1, linetype =2, color="Grey") +
  geom_point( alpha = .8) + 
  theme_classic() +
  scale_color_manual(values = c("royalblue", "red", "#AFAEAD")) +
  labs(x = "Effect size", y= "-Log (FDR)", color = "Label") +
  annotate("text", x = 2, y = 1.3, label = "FDR = 0.1", size = 5, fontface = "plain") +
  geom_text_repel(aes(label=label), data = subset(volcano_q, `Coefficient (UC JP)` > X),
                  nudge_x = 3 - subset(volcano_q, `Coefficient (UC JP)` > X)$`Coefficient (UC JP)` ,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 1,
                  size = 7,
                  color = "Black"
  ) +
  geom_text_repel(aes(label=label), data  = subset(volcano_q, `Coefficient (UC JP)` < -X),
                  nudge_x = -3 - subset(volcano_q, `Coefficient (UC JP)`< -X)$`Coefficient (UC JP)` ,
                  segment.size = 0.2,
                  segment.color = "grey50",
                  direction = "y",
                  hjust = 1,
                  size = 7,
                  color = "Black"
  ) +
  scale_x_continuous(
    breaks = c(-3, -2, -1, 0, 1, 2, 3),
    limits = c(-3, 3)
  ) +
  theme(text = element_text(size = 25))

#Volcano plot creation for CD
volcano %>% mutate(log10q_IBD = -log10(volcano$"Q-value (IBD JP)"), log10q_UC = -log10(volcano$"Q-value (UC JP)"), log10q_CD = -log10(volcano$"Q-value (CD JP)"))-> volcano_q
Up_Hits = head(arrange(volcano_q, -volcano_q$"Coefficient (CD JP)"),10)
Down_Hits = head(arrange(volcano_q, volcano_q$"Coefficient (CD JP)"),10)
Hits=rbind(Up_Hits, Down_Hits)  
volcano_q$label = if_else(volcano_q$ID %in% Hits$ID, volcano_q$ID, "") #Coef included in Top 10 and Bottom 10 are labeled.

X <- 0.9
Y <- 8 

volcano_q %>% mutate(class = case_when(volcano_q$"Coefficient (CD JP)" >= X & log10q_CD >= 1 ~ "Increase", volcano_q$"Coefficient (CD JP)" <= -X & log10q_CD >= 1 ~ "Decrease", TRUE ~ "None")) -> volcano_q  
volcano_q %>% mutate(class2 = case_when(volcano_q$"Coefficient (CD JP)" >= X & log10q_CD >= 7 ~ "Increase", volcano_q$"Coefficient (CD JP)" <= -X & log10q_CD >= 12 ~ "Decrease", TRUE ~ "None")) -> volcano_q
volcano_q %>% filter(class2 == "Increase" | class2 == "Decrease") -> Hit2
volcano_q$label = if_else(volcano_q$ID %in% Hit2$ID, volcano_q$ID, "")　

ggplot(volcano_q, aes(x = `Coefficient (CD JP)`, y = log10q_CD, color=class, size=log10q_CD, label = rownames(ID))) +
  geom_vline(xintercept = c(-X, X), linetype =2, color="Grey") +
  geom_hline(yintercept = 1, linetype =2, color="Grey") +
  geom_point( alpha = .8) + 
  theme_classic() +
  scale_color_manual(values = c("royalblue", "red", "#AFAEAD")) +
  labs(x = "Effect size", y= "-Log (FDR)", color = "Label") +
  annotate("text", x = 2.5, y = 1.3, label = "FDR = 0.1", size = 5, fontface = "plain") +
  geom_text_repel(aes(label=label), data = subset(volcano_q, `Coefficient (CD JP)` > X),
                  nudge_x = 3 - subset(volcano_q, `Coefficient (CD JP)` > X)$`Coefficient (CD JP)` , #3はラベルの位置
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 1,
                  size = 7,
                  color = "Black"
  ) +
  geom_text_repel(aes(label=label), data  = subset(volcano_q, `Coefficient (CD JP)` < -X),
                  nudge_x = -3 - subset(volcano_q, `Coefficient (CD JP)`< -X)$`Coefficient (CD JP)` ,
                  segment.size = 0.2,
                  segment.color = "grey50",
                  direction = "y",
                  hjust = 1,
                  size = 7,
                  color = "Black"
  ) +
  scale_x_continuous(
    breaks = c(-3, -2, -1, 0, 1, 2, 3),
    limits = c(-3, 3)
  ) +
  theme(text = element_text(size = 25))
