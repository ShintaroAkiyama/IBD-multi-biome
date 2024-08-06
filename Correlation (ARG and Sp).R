#Correlation analysis between resistome (ARGs) and bacterial species
library(dplyr)
library(tidyverse)
library(Maaslin2)
library(ComplexHeatmap)
library(progress)
library(circlize)

#Resistome analysis (MaAsLin2 analysis)  
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/ARG")
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% filter(Meta$UC_Dx == 1| Meta$UC_Dx ==0) -> Meta

SP_input <- read.csv("ARG.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (RPKM)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(Meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("IBD_HC2.x", "Age", "sex", "BMI", "IBD_HC2.1", "UC", "Crohn", "UC_Dx", "CD_Dx", "IBD", "IBD_HC2.y", "Number_of_ARG"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2

rownames(Meta) <- NULL
column_to_rownames(Meta, "ID") -> Meta
Meta %>% select(UC_Dx, Age, sex, BMI) -> JP_meta

keep <- apply(SP_input2, 2, mean) > 1E-4 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-4, prevalence > 0.1
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

SP_input <- read.csv("ARG.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (RPKM)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(Meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("IBD_HC2.x", "Age", "sex", "BMI", "IBD_HC2.1", "UC", "Crohn", "UC_Dx", "CD_Dx", "IBD", "IBD_HC2.y", "Number_of_ARG"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2

rownames(Meta) <- NULL　
column_to_rownames(Meta, "ID") -> Meta
Meta %>% select(CD_Dx, Age, sex, BMI) -> JP_meta

keep <- apply(SP_input2, 2, mean) > 1E-4 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-4, prevalence > 0.1
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

SP_input <- read.csv("ARG.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (RPKM)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(Meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("IBD_HC2.x", "Age", "sex", "BMI", "IBD_HC2.1", "UC", "Crohn", "UC_Dx", "CD_Dx", "IBD", "IBD_HC2.y", "Number_of_ARG"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2

rownames(Meta) <- NULL
column_to_rownames(Meta, "ID") -> Meta
Meta %>% select(IBD, Age, sex, BMI) -> JP_meta

keep <- apply(SP_input2, 2, mean) > 1E-4 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-4, prevalence > 0.1
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
JP_MA %>% mutate(ID = str_sub(JP_MA$feature, start = -8, end = -2)) -> JP_MA
JP_MA %>% mutate(ID = as.numeric(JP_MA$ID)) -> JP_MA

#Bacteriome analysis (MaAsLin2 analysis) 
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/SP_mOTU3_JP") 

#UCvsControl
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% filter(Meta$UC_Dx == 1| Meta$UC_Dx ==0) -> Meta
Meta %>% select(UC_Dx, Age, sex, BMI) -> JP_meta

SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
SP_input %>% filter(SP_input$IBD ==2 |SP_input$UC ==1) -> SP_input
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("ID", "IBD", "UC", "CD"))] 

keep <- apply(SP_input2, 2, mean) > 1E-4 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-4, prevalence > 0.1
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

fit_data$results -> JP_maaslin_SP
JP_maaslin_SP %>% filter(JP_maaslin_SP$"metadata" == "UC_Dx") -> JP_maaslin_UC_SP
JP_maaslin_UC_SP <- JP_maaslin_UC_SP[, -which (colnames(JP_maaslin_UC_SP) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=JP_maaslin_UC_SP, "Coefficient (UC JP)" = "coef") -> JP_maaslin_UC_SP
rename(.data=JP_maaslin_UC_SP, "P-value (UC JP)" = "pval") -> JP_maaslin_UC_SP
rename(.data=JP_maaslin_UC_SP, "Q-value (UC JP)" = "qval") -> JP_maaslin_UC_SP

#CDvsControl
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% filter(Meta$CD_Dx == 1| Meta$CD_Dx ==0) -> Meta
Meta %>% select(CD_Dx, Age, sex, BMI) -> JP_meta

SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
SP_input %>% filter(SP_input$IBD ==2 |SP_input$CD ==1) -> SP_input
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("ID", "IBD", "UC", "CD"))] 

keep <- apply(SP_input2, 2, mean) > 1E-4 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-4, prevalence > 0.1
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

fit_data$results -> JP_maaslin_SP
JP_maaslin_SP %>% filter(JP_maaslin_SP$"metadata" == "CD_Dx") -> JP_maaslin_CD_SP
JP_maaslin_CD_SP <- JP_maaslin_CD_SP[, -which (colnames(JP_maaslin_CD_SP) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=JP_maaslin_CD_SP, "Coefficient (CD JP)" = "coef") -> JP_maaslin_CD_SP
rename(.data=JP_maaslin_CD_SP, "P-value (CD JP)" = "pval") -> JP_maaslin_CD_SP
rename(.data=JP_maaslin_CD_SP, "Q-value (CD JP)" = "qval") -> JP_maaslin_CD_SP

#IBDvsControl
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% select(IBD, Age, sex, BMI) -> JP_meta

SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("ID", "IBD", "UC", "CD"))] 

keep <- apply(SP_input2, 2, mean) > 1E-4 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-4, prevalence > 0.1
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

fit_data$results -> JP_maaslin_SP
JP_maaslin_SP %>% filter(JP_maaslin_SP$"metadata" == "IBD") -> JP_maaslin_IBD_SP
JP_maaslin_IBD_SP <- JP_maaslin_IBD_SP[, -which (colnames(JP_maaslin_IBD_SP) %in% c("metadata", "value", "stderr", "name", "N", "N.not.zero"))] 
rename(.data=JP_maaslin_IBD_SP, "Coefficient (IBD JP)" = "coef") -> JP_maaslin_IBD_SP
rename(.data=JP_maaslin_IBD_SP, "P-value (IBD JP)" = "pval") -> JP_maaslin_IBD_SP
rename(.data=JP_maaslin_IBD_SP, "Q-value (IBD JP)" = "qval") -> JP_maaslin_IBD_SP

#Data integration
JP_MA_SP <- full_join(full_join(JP_maaslin_IBD_SP, JP_maaslin_UC_SP, by = "feature"), JP_maaslin_CD_SP, by = "feature")
JP_MA_SP %>% mutate(ID = str_sub(JP_MA_SP$feature, start = -6, end = -2)) -> JP_MA_SP

#Prepare for ARG abundance
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/ARG")
ARG_input <- read.csv("ARG.csv", header = TRUE, na.strings = c(NA, ''), row.names=1, check.names = FALSE) #abundance data (RPKM)
ARG_input %>% select(Number_of_ARG) -> Number_of_ARG
ARG_input <- ARG_input[, -which (colnames(ARG_input) %in% c("IBD_HC2"))] 
ARG_input <- as.data.frame(t(ARG_input)) 
rownames_to_column(ARG_input, "feature") -> ARG_input 
ARG_input %>% mutate(ID = str_sub(ARG_input$feature, start = -8, end = -2)) -> ARG_input
ARG_input %>% mutate(ID = as.numeric(ARG_input$ID)) -> ARG_input

Ref <- read.csv("Reference.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) 
Ref %>% select("ID", "CLASS2", "Resistance.Mechanism") -> Ref2
left_join(JP_MA, Ref2, by = "ID") -> JP_MA

#Select ARG with coef >1 and FDR <0.1 in UC or CD
JP_MA %>% filter(JP_MA$"Coefficient (CD JP)" > 1 & JP_MA$"Q-value (CD JP)" < 0.1 | JP_MA$"Coefficient (UC JP)" > 1 & JP_MA$"Q-value (UC JP)" < 0.1) -> ARG_coef 
ARG_coef %>% select("feature", "Coefficient (IBD JP)", "P-value (IBD JP)", "Q-value (IBD JP)", "Coefficient (UC JP)", "P-value (UC JP)", "Q-value (UC JP)", 
                    "Coefficient (CD JP)", "P-value (CD JP)", "Q-value (CD JP)", "ID", "CLASS2", "Resistance.Mechanism") -> ARG_anno  #For heatmap annotation

inner_join(ARG_input, ARG_coef, by = "ID") -> ARG_input2
ARG_input2 <- ARG_input2[, -which (colnames(ARG_input2) %in% c("ID", "feature.y", "Coefficient (IBD JP)", "P-value (IBD JP)", "Q-value (IBD JP)", 
                                                               "Coefficient (UC JP)", "P-value (UC JP)", "Q-value (UC JP)", "Coefficient (CD JP)", "P-value (CD JP)", "Q-value (CD JP)", 
                                                               "CLASS2", "Resistance.Mechanism"))] #Remove non-essential variables
column_to_rownames(ARG_input2, "feature.x") -> ARG_input2

ARG_input2 <- as.data.frame(t(ARG_input2))
cbind(ARG_input2, Number_of_ARG)->ARG_input2
target_data  <- ARG_input2[,colSums(ARG_input2) != 0] # all 0 omit
target_list <- colnames(target_data) 

#Prepare for bacterial species abundance data
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/SP_mOTU3(METAFは番号に修正済みIBD症状薬剤解析)")
SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1, check.names = FALSE) 
SP_input %>% select(-c(1:4)) -> SP_input
species_data  <- SP_input[,colSums(SP_input) != 0] # all 0 omit
species_list <- colnames(species_data) 

#SP vs ARG Spearman
#make matrix
matrix_rho <- matrix(NA,length(species_list),length(target_list)) 
rownames(matrix_rho) <- species_list 
colnames(matrix_rho) <- target_list 
matrix_pval <- matrix_rho
matrix_n <- matrix_rho
#spearman iteration
pb <- progress_bar$new(total = length(species_list),
                       format = "[:bar] :percent remain: :eta",
                       clear = TRUE)
for (i in 1:length(species_list)){
  pb$tick()
  ##############################
  for (j in 1:length(target_list)){
    #pair setting
    data_1 <- select(species_data, i) %>% as.matrix() %>% as.numeric()
    data_2 <- select(target_data, j) %>% as.matrix() %>% as.numeric()
    n_all <- nrow(target_data)
    #spearman correlation
    spearman_all <- cor.test(data_1, data_2, method = "spearman")
    rho_all <- spearman_all$estimate
    pval_all <- spearman_all$p.value
    #Store the data in the matrix
    matrix_rho[i,j] <- rho_all
    matrix_pval[i,j] <- pval_all
    matrix_n[i,j] <- n_all
    ##############################
    Sys.sleep(1 / length(species_list))
  }
}

matrix_pval <- as.data.frame(matrix_pval) 

#make matrix for FDR-adjusted Pval
matrix_adp <- matrix(NA,length(species_list),length(target_list)) 
rownames(matrix_adp) <- species_list 
colnames(matrix_adp) <- target_list 
pb <- progress_bar$new(total = length(target_list),
                       format = "[:bar] :percent remain: :eta",
                       clear = TRUE)
for (j in 1:length(target_list)){
  pb$tick()
  data_2 <- select(matrix_pval, j) %>% as.matrix() %>% as.numeric()
  spearman_adp <- p.adjust(data_2, method = "fdr")
  #Store the data in the matrix
  matrix_adp[,j] <- spearman_adp
  ##############################
  Sys.sleep(1 / length(target_list))
}

matrix_rho <- as.data.frame(matrix_rho) 
filter_all(matrix_rho, any_vars(abs(.) >0.4)) -> T_rho #Include any variables with coef >0.4

T_rho <- as.data.frame(T_rho) %>% rownames_to_column(var="feature") 
T_rho %>% mutate(ID = str_sub(T_rho$feature, start = -6, end = -2)) -> T_rho

JP_MA_SP %>% mutate(ID = str_sub(JP_MA_SP$feature, start = -6, end = -2)) -> JP_MA_SP
inner_join(T_rho, JP_MA_SP, by = "ID") ->  T_rho
T_rho %>% arrange(-T_rho$"Coefficient (CD JP)") -> T_rho 
rename(.data=T_rho, "feature" = "feature.x") -> T_rho

T_rho <- T_rho[, -which (colnames(T_rho) %in% c("ID", "feature.y", "Coefficient (IBD JP)", 
                                                "P-value (IBD JP)", "Q-value (IBD JP)", "Coefficient (UC JP)", 
                                                "P-value (UC JP)", "Q-value (UC JP)", "Coefficient (CD JP)", 
                                                "P-value (CD JP)", "Q-value (CD JP)"))] #Remove non-essential variables

T_rho %>% select(feature) -> T_rho_p #List of selected bacterial species

matrix_adp <- as.data.frame(matrix_adp) %>% rownames_to_column(var="feature") 
inner_join(T_rho_p, matrix_adp, by = "feature") -> T_p #FDR-adjusted P
column_to_rownames(T_p, var="feature") -> T_p 
column_to_rownames(T_rho, var="feature") -> T_rho　

T_rho %>% select(Number_of_ARG) -> num_rho　#Number of ARG
colnames(num_rho)[1] <- "Number of ARG"
T_p %>% select(Number_of_ARG) -> num_p

T_rho_p %>% mutate(ID = str_sub(T_rho_p$feature, start = -6, end = -2)) -> T_rho_p
inner_join(T_rho_p, JP_MA_SP, by = "ID") -> T_coef #SP annotation
T_coef <- T_coef[, -which (colnames(T_coef) %in% c("ID", "feature.y"))]
rename(.data=T_coef, "feature" = "feature.x") -> T_coef

# Reorder the rho data according to the ARG coef (IBD)
T_rho <- as.data.frame(t(T_rho))
rownames_to_column(T_rho, "feature") -> T_rho 
T_rho %>% mutate(ID = str_sub(T_rho$feature, start = -8, end = -2)) -> T_rho
T_rho %>% mutate(ID = as.numeric(T_rho$ID)) -> T_rho
inner_join(T_rho, ARG_anno, by = "ID") -> T_rho
T_rho %>% arrange(-T_rho$"Coefficient (IBD JP)") -> T_rho
T_rho <- T_rho[, -which (colnames(T_rho) %in% c("ID", "feature.y", "Coefficient (IBD JP)", 
                                                "P-value (IBD JP)", "Q-value (IBD JP)", "Coefficient (UC JP)", 
                                                "P-value (UC JP)", "Q-value (UC JP)", "Coefficient (CD JP)", 
                                                "P-value (CD JP)", "Q-value (CD JP)", "CLASS2", "Resistance.Mechanism"))] #Remove non-essential variables

T_rho %>% select(feature.x) -> T_rho_KO_p #List of selected ARGs
rownames(T_rho) <- NULL 
column_to_rownames(T_rho, var="feature.x") -> T_rho
T_rho <- as.data.frame(t(T_rho))

#Shorten the two ARG names
rename(.data = T_rho, "Bbif_ileS_MUP [ARO:3003730]" ="Bifidobacterium ileS conferring resistance to mupirocin [ARO:3003730]") -> T_rho
rename(.data = T_rho, "Bado_rpoB_RIF [ARO:3004480]" ="Bifidobacterium adolescentis rpoB mutants conferring resistance to rifampicin [ARO:3004480]") -> T_rho

# Reorder the adjusted P data according to the ARG coef (IBD)
T_p <- as.data.frame(t(T_p))
rownames_to_column(T_p, "feature") -> T_p #KOを列として設定
T_p %>% mutate(ID = str_sub(T_p$feature, start = -8, end = -2)) -> T_p
T_p %>% mutate(ID = as.numeric(T_p$ID)) -> T_p
inner_join(T_p, ARG_anno, by = "ID") -> T_p
T_p %>% arrange(-T_p$"Coefficient (IBD JP)") -> T_p
T_p <- T_p[, -which (colnames(T_p) %in% c("ID", "feature.y", "Coefficient (IBD JP)", 
                                          "P-value (IBD JP)", "Q-value (IBD JP)", "Coefficient (UC JP)", 
                                          "P-value (UC JP)", "Q-value (UC JP)", "Coefficient (CD JP)", 
                                          "P-value (CD JP)", "Q-value (CD JP)", "CLASS2", "Resistance.Mechanism"))] #Remove non-essential variables
rownames(T_p) <- NULL　 
column_to_rownames(T_p, var="feature.x") -> T_p
T_p <- as.data.frame(t(T_p))

T_rho_KO_p %>% mutate(ID = str_sub(T_rho_KO_p$feature, start = -8, end = -2)) -> T_rho_KO_p
T_rho_KO_p %>% mutate(ID = as.numeric(T_rho_KO_p$ID)) -> T_rho_KO_p
inner_join(T_rho_KO_p, ARG_anno, by = "ID") -> K_p
rename(.data = K_p, "Class" = "CLASS2") -> K_p
rename(.data = K_p, "Mechanism" = "Resistance.Mechanism") -> K_p

#Heatmap creation
is_sig = K_p$"Q-value (IBD JP)" < 0.1
pch = rep("*", nrow(K_p))
pch[!is_sig] = NA

is_sig = K_p$"Q-value (UC JP)" < 0.1
pch2 = rep("*", nrow(K_p))
pch2[!is_sig] = NA

is_sig = K_p$"Q-value (CD JP)"< 0.1
pch3 = rep("*", nrow(K_p))
pch3[!is_sig] = NA

m_col_fun = colorRamp2(c(-2.5, 0, 2.5), c("navy", "white", "firebrick3")) 
lgd_pvalue = Legend(title = "MaAsLin coef", col_fun = m_col_fun, at = c(-2.5, 0, 2.5), 
                    labels = c("-2.5", "0", "2.5"))
lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")

KO_ra <- HeatmapAnnotation("UC" = anno_simple(K_p$"Coefficient (UC JP)", col = m_col_fun, pch = pch2, simple_anno_size = unit(3,"mm"), gp = gpar(col = "black")), "CD" = anno_simple(K_p$"Coefficient (CD JP)", col = m_col_fun, pch = pch3, simple_anno_size = unit(3,"mm"), gp = gpar(col = "black")), 
                           "Class"=K_p$Class, "Resistance mechanism"=K_p$Mechanism,   col=list(Class=c("AMINOGLYCOSIDE" = "#ff4b00", "AMPHENICOLS"="#fff100", "ANTIBIOTIC EFFLUX PUMP"="#03af7a", "DRUGS FOR TUBERCULOSIS"="#005aff", "MACROLIDES, LINCOSAMIDES AND STREPTOGRAMINS" = "#4dc4ff", "OTHER ANTIBACTERIALS"="#ffcabf", 
                                                                                                       "OTHER BETA-LACTAM"="#f6aa00", "OTHERS" = "#c8c8cb", "SULFONAMIDES & TRIMETHOPRIM" = "#990099", "TETRACYCLINES" = "#804000", "QUINOLONE" = "#c9ace6"), 
                                                                                               "Resistance mechanism"= c("antibiotic efflux" = "#ff4b00", "antibiotic efflux;antibiotic target alteration"="#fff100", "antibiotic efflux;reduced permeability to antibiotic"="#03af7a", "antibiotic inactivation"="#005aff", 
                                                                                                                         "antibiotic target alteration" = "#4dc4ff", "antibiotic target protection"="#ffcabf", "antibiotic target alteration;antibiotic target replacement"="#f6aa00", "antibiotic target replacement" = "#990099","NoMatch"="#c8c8cb")),  
                           simple_anno_size = unit(3,"mm"), annotation_name_gp= gpar(fontsize = 8), gp = gpar(col = "black"))

n2=pheatmap(as.matrix(num_rho), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(num_p<0.1,"*", ""), nrow(num_p)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-1,0,1), c("#03af7a", "white", "#f6aa00")), name = "Spearman rho",  heatmap_legend_param = list(color_bar = "continuous"))

T_coef %>% select("Q-value (IBD JP)") -> qval_IBD
T_coef %>% select("Q-value (UC JP)", "Q-value (CD JP)") -> qval_UC_CD
qval_IBD[is.na(qval_IBD)] <- 1
qval_UC_CD[is.na(qval_UC_CD)] <- 1

T_coef %>% select("Coefficient (IBD JP)") -> coef_IBD
T_coef %>% select("Coefficient (UC JP)", "Coefficient (CD JP)") ->  coef_UC_CD
data.frame(row.names = T_coef$feature, coef_IBD) -> coef_IBD_2
data.frame(row.names = T_coef$feature, coef_UC_CD) -> coef_UC_CD2

rename(.data = coef_UC_CD2, "UC" = "Coefficient..UC.JP.") -> coef_UC_CD2
rename(.data = coef_UC_CD2, "CD" = "Coefficient..CD.JP.") -> coef_UC_CD2

p2=pheatmap(as.matrix(coef_UC_CD2), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(qval_UC_CD<0.1,"*", ""), nrow(qval_UC_CD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-2.5,0,2.5), c("navy", "white", "firebrick3")), name = "MaAsLin coef", heatmap_legend_param = list(color_bar = "continuous"))
p3=pheatmap(as.matrix(T_rho), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(T_p <0.1,"*", ""), nrow(T_p)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-1,0,1), c("#03af7a", "white", "#f6aa00")), name = "Spearman rho", heatmap_legend_param = list(color_bar = "continuous"), column_title = "ARGs", column_title_gp = gpar(fontsize = 10, fontface = "bold"), top_annotation = KO_ra)

draw(p2+p3+n2, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig))  
