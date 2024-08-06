# Function analysis of bacteriome (Figure 3c-3h)
library(dplyr)
library(tidyverse)
library(Maaslin2)
library(ComplexHeatmap)

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Function")
#UCvsControl
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta %>% filter(Meta$UC_Dx == 1| Meta$UC_Dx ==0) -> Meta

SP_input <- read.csv("KO.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(Meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("IBD_HC2.x", "Age", "sex", "BMI", "IBD_HC2.1", "UC", "Crohn", "UC_Dx", "CD_Dx", "IBD", "IBD_HC2.y"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2

rownames(Meta) <- NULL　
column_to_rownames(Meta, "ID") -> Meta
Meta %>% select(UC_Dx, Age, sex, BMI) -> JP_meta

keep <- apply(SP_input2, 2, mean) > 1E-8 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-8, prevalence > 0.1
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

SP_input <- read.csv("KO.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(Meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("IBD_HC2.x", "Age", "sex", "BMI", "IBD_HC2.1", "UC", "Crohn", "UC_Dx", "CD_Dx", "IBD", "IBD_HC2.y"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2

rownames(Meta) <- NULL　
column_to_rownames(Meta, "ID") -> Meta
Meta %>% select(CD_Dx, Age, sex, BMI) -> JP_meta

keep <- apply(SP_input2, 2, mean) > 1E-8 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-8, prevalence > 0.1
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

SP_input <- read.csv("KO.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data (normalized)
rownames_to_column(SP_input, "ID") -> SP_input
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
left_join(Meta, SP_input, by = "ID") -> SP_input2
SP_input2 <- SP_input2[, -which (colnames(SP_input2) %in% c("IBD_HC2.x", "Age", "sex", "BMI", "IBD_HC2.1", "UC", "Crohn", "UC_Dx", "CD_Dx", "IBD", "IBD_HC2.y"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2

rownames(Meta) <- NULL　
column_to_rownames(Meta, "ID") -> Meta
Meta %>% select(IBD, Age, sex, BMI) -> JP_meta

keep <- apply(SP_input2, 2, mean) > 1E-8 & apply(SP_input2 > 0, 2, sum) / nrow(SP_input2) > 0.1 # mean abundance > 1E-8, prevalence > 0.1
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

KEGG_input <- read.csv("220430_KEGG_56700.csv", header = TRUE, na.strings = c(NA, '')) #KEGG list
rename(.data=KEGG_input, "feature" = "k") -> KEGG_input
left_join(KEGG_input, JP_MA, by = "feature") -> KO_coef

#Select KOs with abs coef>1 and FDR<0.1 in UC or CD in the specific pathmay (Figure 3c-3d)
KO_coef %>% filter(abs(KO_coef$"Coefficient (CD JP)") > 1 & KO_coef$"Q-value (CD JP)" < 0.1 | abs(KO_coef$"Coefficient (UC JP)") > 1 & KO_coef$"Q-value (UC JP)" < 0.1) -> KO_coef 
KO_coef %>% filter(KO_coef$C=="Terpenoid backbone biosynthesis [PATH:ko00900]") -> KO_coef #"Cationic antimicrobial peptide (CAMP) resistance [PATH:ko01503]", "N-Glycan biosynthesis [PATH:ko00510]", "Glycosaminoglycan degradation [PATH:ko00531]" can be selected.
KO_coef %>% arrange(-KO_coef$"Coefficient (IBD JP)") -> KO_coef

KO_coef %>% distinct(feature, .keep_all=TRUE) -> KO_coef #Remove duplicates
column_to_rownames(KO_coef, "feature") -> KO_coef 

#Heatmap creation
KO_coef %>% select("Q-value (IBD JP)") -> qval_IBD
KO_coef %>% select("Q-value (UC JP)", "Q-value (CD JP)") -> qval_UC_CD
qval_IBD[is.na(qval_IBD)] <- 1
qval_UC_CD[is.na(qval_UC_CD)] <- 1

KO_coef %>% select("Coefficient (IBD JP)") -> coef_IBD
KO_coef %>% select("Coefficient (UC JP)", "Coefficient (CD JP)") ->  coef_UC_CD

colnames(coef_IBD)[1] <- "IBD"
colnames(coef_UC_CD)[1] <- "UC"
colnames(coef_UC_CD)[2] <- "CD"

anno_width = unit(2, "cm")

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")

p1=pheatmap(as.matrix(coef_IBD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(qval_IBD < 0.1,"*", ""), nrow(qval_IBD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-2.5,0,2.5), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))
p2=pheatmap(as.matrix(coef_UC_CD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(qval_UC_CD < 0.1,"*", ""), nrow(qval_UC_CD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-2.5,0,2.5), c("navy", "white", "firebrick3")), column_title = "Terpenoid backbone biosynthesis [PATH:ko00900]", name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))
#"Cationic antimicrobial peptide (CAMP) resistance [PATH:ko01503]", "N-Glycan biosynthesis [PATH:ko00510]", "Glycosaminoglycan degradation [PATH:ko00531]" can be selected in the column_title.

draw(p1+p2, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig))

#PAMPS analysis (Figure 3g)
KEGG_input <- read.csv("220430_KEGG_56700.csv", header = TRUE, na.strings = c(NA, '')) #KEGG list
rename(.data=KEGG_input, "feature" = "k") -> KEGG_input
left_join(KEGG_input, JP_MA, by = "feature") -> KO_coef

KO_coef %>% filter(KO_coef$"Coefficient (CD JP)" > 1 & KO_coef$"Q-value (CD JP)" < 0.1 | KO_coef$"Coefficient (UC JP)" > 1 & KO_coef$"Q-value (UC JP)" < 0.1) -> KO_coef 
KO_coef %>% filter(KO_coef$C=="Flagellar assembly [PATH:ko02040]" | KO_coef$C=="Lipopolysaccharide biosynthesis [PATH:ko00540]"| KO_coef$C=="Peptidoglycan biosynthesis [PATH:ko00550]") -> KO_coef #Select PAMPs
KO_coef %>% arrange(-KO_coef$"Coefficient (IBD JP)") -> KO_coef
KO_coef %>% arrange(KO_coef$C) -> KO_coef

KO_coef %>% distinct(feature, .keep_all=TRUE) -> KO_coef #Remove duplicates
column_to_rownames(KO_coef, "feature") -> KO_coef 

#Heatmap creation
KO_coef %>% select("Q-value (IBD JP)") -> qval_IBD
KO_coef %>% select("Q-value (UC JP)", "Q-value (CD JP)") -> qval_UC_CD
qval_IBD[is.na(qval_IBD)] <- 1
qval_UC_CD[is.na(qval_UC_CD)] <- 1

KO_coef %>% select("Coefficient (IBD JP)") -> coef_IBD
KO_coef %>% select("Coefficient (UC JP)", "Coefficient (CD JP)") ->  coef_UC_CD

colnames(coef_IBD)[1] <- "IBD"
colnames(coef_UC_CD)[1] <- "UC"
colnames(coef_UC_CD)[2] <- "CD"

anno_width = unit(2, "cm")

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")

p1=pheatmap(as.matrix(coef_IBD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(qval_IBD < 0.1,"*", ""), nrow(qval_IBD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-2,0,2), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))
p2=pheatmap(as.matrix(coef_UC_CD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(qval_UC_CD < 0.1,"*", ""), nrow(qval_UC_CD)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-2,0,2), c("navy", "white", "firebrick3")), column_title = "PAMPs", name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))

ra = rowAnnotation("Function"=KO_coef$C, col=list(Function=c("Flagellar assembly [PATH:ko02040]" = "#ff4b00", "Lipopolysaccharide biosynthesis [PATH:ko00540]"="#fff100", "Peptidoglycan biosynthesis [PATH:ko00550]"="#03af7a")), simple_anno_size = unit(3, "mm"), gp = gpar(col = "black"), annotation_name_gp= gpar(fontsize = 8))

draw(ra+p1+p2, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig))

#AIEC virulence factor analysis (Figure 3h)
KEGG_input <- read.csv("220430_KEGG_56700.csv", header = TRUE, na.strings = c(NA, '')) #KEGG list
rename(.data=KEGG_input, "feature" = "k") -> KEGG_input
left_join(KEGG_input, JP_MA, by = "feature") -> KO_coef

AIEC <- read.csv("AIEC.csv", header = TRUE, na.strings = c(NA, '')) 
AIEC  %>% filter(!is.na(AIEC)) -> AIEC
AIEC <- AIEC[, -which (colnames(AIEC) %in% c("A", "B", "C", "D", "rpt"))] 

AIEC%>% distinct(feature, .keep_all=TRUE) -> AIEC #Remove duplicates

left_join(KO_coef, AIEC, by = "feature") -> KO_coef

KO_coef  %>% filter(!is.na(AIEC)) -> KO_coef
KO_coef %>% arrange(AIEC) -> KO_coef
KO_coef %>% distinct(feature, .keep_all=TRUE) -> KO_coef
column_to_rownames(KO_coef, "feature") -> KO_coef
KO_coef$Function[is.na(KO_coef$Function)] <- "Others"

KO_coef %>% select("Q-value (IBD JP)") -> qval_IBD
KO_coef %>% select("Q-value (UC JP)", "Q-value (CD JP)") -> qval_UC_CD
qval_IBD[is.na(qval_IBD)] <- 1
qval_UC_CD[is.na(qval_UC_CD)] <- 1

KO_coef %>% select("Coefficient (IBD JP)") -> coef_IBD
KO_coef %>% select("Coefficient (UC JP)", "Coefficient (CD JP)") ->  coef_UC_CD
KO_coef  %>% select(AIEC) -> AIEC
KO_coef  %>% select(Function) -> Function

colnames(coef_IBD)[1] <- "IBD"
colnames(coef_UC_CD)[1] <- "UC"
colnames(coef_UC_CD)[2] <- "CD"

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")

ra = rowAnnotation("Function"=KO_coef$Function, col=list(Function=c("Flagellar assembly" = "#c9ace6", "Fimbriae type 1"="#ff8082", "Type VI secretion"="#ffff80", "Others"="#d8f255")), simple_anno_size = unit(3, "mm"), gp = gpar(col = "black"), annotation_name_gp= gpar(fontsize = 8))
ra2 = rowAnnotation("AIEC virulence factors"=KO_coef$AIEC, col=list("AIEC virulence factors"=c("Crossing of mucous layer and host defense peptides" = "#fff100", "Interaction with Peyer's patches and macrophages survival"="#03af7a", "Metabolic adaptation and global regulation"="#990099", "Adhesion to and invasion of IECs"="#ff4b00")), simple_anno_size = unit(3, "mm"), gp = gpar(col = "black"), annotation_name_gp= gpar(fontsize = 8))

p1=pheatmap(as.matrix(coef_IBD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(qval_IBD<0.1,"*", ""), nrow(qval_IBD)), fontsize_number = 7, border_color = TRUE, col = circlize::colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")), name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))
p2=pheatmap(as.matrix(coef_UC_CD), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(qval_UC_CD<0.1,"*", ""), nrow(qval_UC_CD)), fontsize_number = 7, border_color = TRUE, col = circlize::colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")), column_title = "AIEC", name = "Effect size", heatmap_legend_param = list(color_bar = "continuous"))

draw(ra2+ra+p1+p2,  heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig))
