#ARG analysis with CARD annotation (Supplementary Figure 7)
library(dplyr)
library(tidyverse)
library(Maaslin2)
library(ComplexHeatmap)
library(progress)
library(reshape2)

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/ARG")
#UCvsControl
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

Ref <- read.csv("Reference.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) 
#Ref %>% select("ID", "CLASS2", "Resistance.Mechanism") -> Ref2
left_join(JP_MA, Ref, by = "ID") -> JP_MA
JP_MA %>% arrange(-JP_MA$"Coefficient (IBD JP)") -> JP_MA

#Select ARGs with abs coef >1 and FDR <0.1 in UC or CD
JP_MA %>% filter(abs(JP_MA$"Coefficient (CD JP)") > 1 & JP_MA$"Q-value (CD JP)" < 0.1 | abs(JP_MA$"Coefficient (UC JP)") > 1& JP_MA$"Q-value (UC JP)" < 0.1) -> ARG_Heat_qval

species_table_joined <- ARG_Heat_qval[, -which (colnames(ARG_Heat_qval) %in% c("feature", "Coefficient (IBD JP)", "P-value (IBD JP)", "Q-value (IBD JP)", "Coefficient (UC JP)", 
                                                                               "P-value (UC JP)", "Q-value (UC JP)", "Coefficient (CD JP)", "P-value (CD JP)", "Q-value (CD JP)", 
                                                                               "ID", "CLASS", "CLASS2", "Resistance.Mechanism"))] #Remove non-essential variables

column_to_rownames(ARG_Heat_qval, "ARG") -> ARG_Heat_qval
ARG_Heat_qval %>% select("Coefficient (IBD JP)") -> coef_IBD
ARG_Heat_qval %>% select("Coefficient (UC JP)", "Coefficient (CD JP)") -> coef_UC_CD
ARG_Heat_qval %>% select( "Q-value (IBD JP)") -> qval_IBD
ARG_Heat_qval %>% select( "Q-value (UC JP)", "Q-value (CD JP)") -> qval_UC_CD

#Heatmap annotation creation (bar plot for the percentage of significance)
coef_IBD[is.na(coef_IBD)] <- 1 
coef_UC_CD[is.na(coef_UC_CD)] <- 1
qval_IBD[is.na(qval_IBD)] <- 1
qval_UC_CD[is.na(qval_UC_CD)] <- 1
#When calculating positive and negative significant variables, a variable with a NA coef is counted as one.

Num_significant_T_IBD <- as.data.frame(colnames(coef_IBD)) 
colnames(Num_significant_T_IBD) <- "Disease"

Num_significant_T_IBD$positive <- 0
Num_significant_T_IBD$negative <- 0

for (x in 1:nrow(coef_IBD)){
  for (y in 1:ncol(coef_IBD)){
    if(qval_IBD[x,y] < 0.1){
      if(coef_IBD[x,y] > 0){
        Num_significant_T_IBD$positive[y] <- Num_significant_T_IBD$positive[y]+1
      }
      else{Num_significant_T_IBD$negative[y] <- Num_significant_T_IBD$negative[y]+1
      }
    }
  }
}
Num_significant_T_IBD <- column_to_rownames(Num_significant_T_IBD, var="Disease")
Num_significant_T_IBD <- Num_significant_T_IBD/nrow(coef_IBD)*100

Num_significant_T_UC_CD <- as.data.frame(colnames(coef_UC_CD)) 
colnames(Num_significant_T_UC_CD) <- "Disease"

Num_significant_T_UC_CD$positive <- 0
Num_significant_T_UC_CD$negative <- 0

for (x in 1:nrow(coef_UC_CD)){
  for (y in 1:ncol(coef_UC_CD)){
    if(qval_UC_CD[x,y] < 0.1){
      if(coef_UC_CD[x,y] > 0){
        Num_significant_T_UC_CD$positive[y] <- Num_significant_T_UC_CD$positive[y]+1
      }
      else{Num_significant_T_UC_CD$negative[y] <- Num_significant_T_UC_CD$negative[y]+1
      }
    }
  }
}

Num_significant_T_UC_CD <- column_to_rownames(Num_significant_T_UC_CD, var="Disease")
Num_significant_T_UC_CD <- Num_significant_T_UC_CD/nrow(coef_UC_CD)*100

ra_IBD <-  HeatmapAnnotation(Num_sig_IBD = anno_barplot(Num_significant_T_IBD, ylim = c(0, 100), border = FALSE, gp = gpar(fill = c("firebrick3", "blue")), 
                                                        axis_param = list(gp = gpar(fontsize = 6))), empty0 = anno_empty(border = FALSE, height = unit(0.1, "cm")), 
                             empty1 = anno_empty(border = FALSE, height = unit(0.1, "cm")), annotation_name_gp = gpar(fontsize = 0))

ra_UC_CD <- HeatmapAnnotation(Percentage_of_significance = anno_barplot(Num_significant_T_UC_CD, ylim = c(0, 100), border = FALSE, gp = gpar(fill = c("firebrick3", "blue")), 
                                                                        axis_param = list(gp = gpar(fontsize = 0))), empty0 = anno_empty(border = FALSE, height = unit(0.1, "cm")), 
                              empty1 = anno_empty(border = FALSE, height = unit(0.1, "cm")), annotation_name_gp = gpar(fontsize = 7))

lgd_list = list(Legend(labels = c("Positive", "Negative"), title = "Percentage of significance", legend_gp = gpar(fill = c("firebrick3", "blue"))))

#Heatmap annotation of CARD to show species that possess the ARG (Sequence variance)
ARG_Heat_qval %>% select("Coefficient (IBD JP)") -> ARG_coef_IBD
rename(.data= ARG_coef_IBD, "IBD" = "Coefficient (IBD JP)") -> ARG_coef_IBD

ARG_Heat_qval %>% select("Coefficient (UC JP)", "Coefficient (CD JP)") -> ARG_coef_UC_CD
rename(.data= ARG_coef_UC_CD, "UC" = "Coefficient (UC JP)") -> ARG_coef_UC_CD
rename(.data= ARG_coef_UC_CD, "CD" = "Coefficient (CD JP)") -> ARG_coef_UC_CD

ARG_Heat_qval %>% select( "Q-value (IBD JP)") -> qval_IBD
ARG_Heat_qval %>% select( "Q-value (UC JP)", "Q-value (CD JP)") -> qval_UC_CD
qval_IBD[is.na(qval_IBD)] <- 1
qval_UC_CD[is.na(qval_UC_CD)] <- 1 

#Shorten the two ARG names 
rownames(ARG_coef_IBD)[2] <- "Bbif_ileS_MUP [ARO:3003730]" #"Bifidobacterium ileS conferring resistance to mupirocin [ARO:3003730]"
rownames(ARG_coef_IBD)[3] <- "Bado_rpoB_RIF [ARO:3004480]" #"Bifidobacterium adolescentis rpoB mutants conferring resistance to rifampicin [ARO:3004480]"

rownames(ARG_coef_UC_CD)[2] <- "Bbif_ileS_MUP [ARO:3003730]"
rownames(ARG_coef_UC_CD)[3] <- "Bado_rpoB_RIF [ARO:3004480]" 

species_table_joined %>% select(ARG) -> ARG_list

species_count <- species_table_joined %>% column_to_rownames(var="ARG") %>% as.matrix() %>% 
  melt() %>% select(3) %>% table() %>% as.data.frame()

#Species below the threshold (Freq < 10) are rewritten to "Others"
Minorsp_list = filter(species_count, Freq <= 10) %>% select(1) %>% as.matrix()
n <- nrow(ARG_list)
pb <- progress_bar$new(total = n,
                       format = "[:bar] :percent remain: :eta",
                       clear = TRUE)
for (i in 1:n){
  pb$tick()
  for (j in 2:158){
    if (species_table_joined[i, j] %in% Minorsp_list){
      species_table_joined[i, j] <- "Others"
    }
  }
  Sys.sleep(1/n)
}

#Count the number of occurrence of the genus again
species_count <- species_table_joined %>% column_to_rownames(var="ARG") %>% as.matrix() %>% 
  melt() %>% select(3) %>% table() %>% as.data.frame()
species_list <- as.character(species_count[order(species_count$Freq, decreasing = TRUE), 1])
species_list <- c(species_list[-which(species_list %in% "Others")], "Others")

#species_0_1 conversion
annotation_species <- ARG_list %>% as.data.frame()
for (i in 1:length(species_list)){
  for (j in 1:nrow(ARG_list)){
    target <- species_table_joined[j, 2:158] #Location of the CARD data (Sequence variance)
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

ARG_Heat_qval %>% select(CLASS2) -> anottation_2
rename(.data = anottation_2, Class = CLASS2) -> anottation_2

ARG_Heat_qval %>% select(Resistance.Mechanism) -> anottation_3
rename(.data = anottation_3, Mechanism = Resistance.Mechanism) -> anottation_3

#Heatmap creation
p1=pheatmap(as.matrix(ARG_coef_IBD), fontsize = 8, cellwidth = 6, cellheight = 6, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 7, fontsize_row = 6, display_numbers = matrix(ifelse(qval_IBD<0.1,"*", ""), nrow(qval_IBD)), fontsize_number = 6, border_color = "black", col = circlize::colorRamp2(c(-2,0,2), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"), top_annotation = ra_IBD)
p2=pheatmap(as.matrix(ARG_coef_UC_CD), fontsize = 8, cellwidth = 6, cellheight = 6, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 7, fontsize_row = 6, display_numbers = matrix(ifelse(qval_UC_CD<0.1,"*", ""), nrow(qval_UC_CD)), fontsize_number = 6, border_color = "black", col = circlize::colorRamp2(c(-2,0,2), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"), top_annotation = ra_UC_CD)

ra = rowAnnotation(df=anottation_2, col = list(Class=c("AMINOGLYCOSIDE" = "#ff4b00", "AMPHENICOLS"="#fff100", "ANTIBIOTIC EFFLUX PUMP"="#03af7a", "DRUGS FOR TUBERCULOSIS"="#005aff", "MACROLIDES, LINCOSAMIDES AND STREPTOGRAMINS" = "#4dc4ff", "OTHER ANTIBACTERIALS"="#ffcabf", "OTHER BETA-LACTAM"="#f6aa00", "OTHERS" = "#c8c8cb", "SULFONAMIDES & TRIMETHOPRIM" = "#990099", "TETRACYCLINES" = "#804000", "QUINOLONE" = "#c9ace6")), simple_anno_size = unit(2, "mm"), gp = gpar(col = "black"), annotation_name_gp= gpar(fontsize = 7))
ra3 = rowAnnotation(df=anottation_3, col = list(Mechanism=c("antibiotic efflux" = "#ff4b00", "antibiotic efflux;antibiotic target alteration"="#fff100", "antibiotic efflux;reduced permeability to antibiotic"="#03af7a", "antibiotic inactivation"="#005aff", "antibiotic target alteration" = "#4dc4ff", "antibiotic target protection"="#ffcabf", "antibiotic target alteration;antibiotic target replacement"="#f6aa00", "antibiotic target replacement" = "#990099", "NoMatch" = "#c8c8cb")), simple_anno_size = unit(2, "mm"), gp = gpar(col = "black"), annotation_name_gp= gpar(fontsize = 7))

Legend_list <- list(`Escherichia coli` = list(title = "Species that possess the ARG", at = c("infect", "not infect"), labels = c("YES", "NO")))

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")

annotation_species <- annotation_species[, -which (colnames(annotation_species) %in% c("feature"))] 

ra2 = rowAnnotation(df=annotation_species, col = list('Klebsiella pneumoniae' = col_species, `Escherichia coli` = col_species,`Salmonella enterica` = col_species,`Enterobacter hormaechei` = col_species,
                                                      `Shigella flexneri` = col_species, `Shigella sonnei` = col_species, `Escherichia fergusonii` = col_species, `Pseudomonas aeruginosa` = col_species, 
                                                      `Shigella boydii` = col_species, `Acinetobacter baumannii` = col_species, `Shigella dysenteriae` = col_species, `Citrobacter freundii` = col_species, 
                                                      `Klebsiella quasipneumoniae` = col_species, `Serratia marcescens` = col_species, `Escherichia albertii` = col_species, `Klebsiella aerogenes` = col_species, 
                                                      `Citrobacter werkmanii` = col_species, `Staphylococcus aureus` = col_species, `Citrobacter portucalensis` = col_species, `Enterobacter cloacae` = col_species, 
                                                      `Escherichia marmotae` = col_species,`Klebsiella michiganensis` = col_species,`Enterobacter kobei` = col_species,`Klebsiella oxytoca` = col_species,
                                                      `Enterobacter asburiae` = col_species,`Enterobacter roggenkampii` = col_species,`Citrobacter amalonaticus` = col_species,`Citrobacter youngae` = col_species,
                                                      `Enterococcus faecium` = col_species, `Proteus mirabilis` = col_species,`Enterococcus faecalis` = col_species,`Raoultella planticola` = col_species,
                                                      `Streptococcus suis` = col_species,`Citrobacter koseri` = col_species,`Morganella morganii` = col_species,`Enterobacter chengduensis` = col_species,
                                                      `Enterococcus faecium` = col_species,`Leclercia adecarboxylata` = col_species, `Stenotrophomonas maltophilia` = col_species,`Vibrio cholerae` = col_species,
                                                      `Vibrio vulnificus` = col_species, `Bacteroides fragilis` = col_species, `Enterobacter cancerogenus` = col_species,`Bacteroides thetaiotaomicron` = col_species,
                                                      `Salmonella bongori` = col_species,`Mycobacterium avium` = col_species,`Cronobacter sakazakii` = col_species, `Cronobacter dublinensis` = col_species, 
                                                      `Klebsiella huaxiensis` = col_species, `Cronobacter malonaticus` = col_species, `Others` = col_species), gp = gpar(col = "black"), simple_anno_size = unit(2, "mm"),  
                    annotation_name_gp= gpar(fontsize = 7),  show_legend = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,FALSE, FALSE, FALSE, FALSE, FALSE,FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
                                                                             FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE, FALSE,FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
                                                                             FALSE,FALSE, FALSE, FALSE, FALSE), annotation_legend_param = Legend_list)

draw(ra2+ra+ra3+p1+p2,  heatmap_legend_list=lgd_list, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig))
