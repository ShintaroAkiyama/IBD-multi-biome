library(dplyr)
library(tidyverse)
library(NetCoMi)
library(ComplexHeatmap)
library(progress)

#UC analysis
#Bacteria list
#SP_input2_bacteriome from bacteriome analysis including bacterial species with abs coef >1 and FDR<0.1 in UC or CD
rownames_to_column(SP_input2_bacteriome, "feature") -> SP_input2_bacteriome 

#Repeat from here
SP_input2_bacteriome %>% filter(abs(SP_input2_bacteriome$"Coefficient (UC JP)") > 1 & SP_input2_bacteriome$"Q-value (UC JP)" < 0.1) -> SP_input2_bacteriome2 #For UC
SP_input2_bacteriome2 %>% dplyr::select(feature) -> Bacteriome_No 
#Phage list 
#vOTU_Heat_qval_fig from virome analysis including vOTU with abs coef >1 and FDR<0.1 in UC or CD
vOTU_Heat_qval_fig %>% filter(abs(vOTU_Heat_qval_fig$"Coefficient (UC JP)") > 1 & vOTU_Heat_qval_fig$"Q-value (UC JP)" < 0.1) -> vOTU_Heat_qval_fig2 #For UC
rownames_to_column(vOTU_Heat_qval_fig2, "No") -> vOTU_Heat_qval_fig2
vOTU_Heat_qval_fig2 %>% dplyr::select("No", "feature") -> Virome_No 

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
Ref <- read.csv("Reference2.csv", header = TRUE, na.strings = c(NA, ''))
left_join(Virome_No, Ref, by = "No") ->Virome_No
Virome_No <- Virome_No[, -which (colnames(Virome_No) %in% c("No"))] 
rename(.data=Virome_No, "No" = "vOTU_host") -> Virome_No

#Fungi list
#Fungi_map from mycobiome analysis including fungi with P<0.05 in UC or CD
Fungi_map %>% filter(Fungi_map $"P-value (UC JP)"< 0.05) -> Fungi_map2 #For UC
Fungi_map2 %>% dplyr::select("feature") -> Fungi_No

# Preparation for bacterial species abundance (Japanese 4D cohort, please skip this script to external data if JP 4D is not used)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/SP_mOTU3(METAFは番号に修正済みIBD症状薬剤解析)")
#UCvsControl
#Bacteriome abundance (JP)
SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1, check.names = FALSE) #abundance data
rownames_to_column(SP_input, "ID") ->SP_input
SP_input %>% filter(SP_input$IBD ==2) -> SP_input #For Cont
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
NO <- SP_input%>%dplyr::select(ID) #UC vs Cont is selected.
rename(.data=NO, "NO" = "ID") -> NO
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("METAF", "IBD", "UC", "CD"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2
SP_input2 <- as.data.frame(t(SP_input2)) 
rownames_to_column(SP_input2, "feature") -> SP_input2 
left_join(Bacteriome_No, SP_input2, by = "feature") -> SP_input2
column_to_rownames(SP_input2, "feature") -> Bacteriome
Bacteriome <- as.data.frame(t(Bacteriome))

#vOTU abundance (JP)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
vOTU_input <- read.csv("vOTU.csv", header = TRUE, na.strings = c(NA, '')) #abundance data
left_join(NO, vOTU_input, by = "NO") -> vOTU_input 
vOTU_input <- vOTU_input[, -which (colnames(vOTU_input) %in% c("IBD_HC"))] 
column_to_rownames(vOTU_input, "NO") -> vOTU_input
vOTU_input <- as.data.frame(t(vOTU_input)) 
rownames_to_column(vOTU_input, "feature") -> vOTU_input 
left_join(Virome_No, vOTU_input, by = "feature") -> vOTU_input2
vOTU_input2 <- vOTU_input2[, -which (colnames(vOTU_input2) %in% c("feature"))] #Remove non-essential variables。
column_to_rownames(vOTU_input2, "No") -> Virome
Virome <- as.data.frame(t(Virome))

#Fungi abundance (JP)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Fungi")
F_input <- read.csv("feces_230925.eukaryote.csv", header = TRUE, na.strings = c(NA, ''), row.names=1, check.names = FALSE) 
as.data.frame(t(F_input)) -> F_input
rownames_to_column(F_input, var = "NO") -> F_input
F_input %>% mutate(NO = as.numeric(F_input$NO)) -> F_input
left_join(NO, F_input, by = "NO") -> F_input #Select UC vs Cont
column_to_rownames(F_input, "NO") -> F_input
F_input <- as.data.frame(t(F_input)) 
rownames_to_column(F_input, "feature") -> F_input
Fungi_No$feature <- str_replace_all(Fungi_No$feature, pattern="[.]", replacement=" ")
left_join(Fungi_No, F_input, "feature") -> F_input2
column_to_rownames(F_input2, "feature") -> Mycobiome
Mycobiome <- as.data.frame(t(Mycobiome))

#Network analysis
rownames_to_column(Bacteriome, "ID") -> Bacteriome
rownames_to_column(Virome, "ID") -> Virome
rownames_to_column(Mycobiome, "ID") -> Mycobiome
inner_join(inner_join(Bacteriome, Virome, by = "ID"), Mycobiome, by = "ID") -> SparData
column_to_rownames(.data = SparData, "ID") -> SparData


#Select important strains for the validation
SparData %>% select("vOTU160(Ruminococcus, Eubacterium)", "vOTU161(Ruminococcus, Eubacterium)", "vOTU162(Ruminococcus, Eubacterium)", "vOTU164(Ruminococcus, Eubacterium)", "Ruminococcus bromii [ID:00853]",
                    "vOTU870(Bacteroides, Alistipes)", "vOTU913(Bacteroides)", "Bacteroides thetaiotaomicron [ID:01657]",
                    "vOTU289(Blautia, Fusicatenibacter)", "vOTU291(Blautia, Fusicatenibacter)", "vOTU293(Blautia, Fusicatenibacter)", "vOTU294(Blautia, Fusicatenibacter)", "Fusicatenibacter saccharivorans [ID:02680]",
                    "vOTU497(Clostridium, Anaerotruncus, Ruminococcus, Blautia, Dorea, Coprococcus)", "Dorea longicatena [ID:03693]", "Ruminococcaceae species incertae sedis [ID:12239]", "Blautia massiliensis [ID:03342]", 
                    "Clostridiales species incertae sedis [ID:12476]", "Clostridiales species incertae sedis [ID:12554]",
                    "vOTU331(Faecalibacterium)", "Faecalibacterium prausnitzii [ID:06112]") -> SparData

net <- netConstruct(SparData, #relative abundance
                    filtTax = "none",
                    filtSamp = "none", 
                    filtSampPar = "none", 
                    measure = "spearman",
                    thresh = 0.2,
                    measurePar = list(nlambda=10, 
                                      rep.num=10),　
                    normMethod = "none",
                    sparsMethod = "t-test", #sparsMethod = "threshhold", thresh = 0.01
                    adjust = "fdr",
                    alpha = 0.1, 
                    lfdrThresh = 0.1,
                    sampleSize = 540,  
                    dissFunc = "signed",
                    verbose = 2,
                    seed = 123456)

net$assoMat1 %>% as.data.frame() -> col_fdr #if fdr<0.1, coef included, otherwise 0

net2 <- netConstruct(SparData, #relative abundance
                     filtTax = "none",
                     filtSamp = "none", 
                     filtSampPar = "none", 
                     measure = "spearman",
                     measurePar = list(nlambda=10, 
                                       rep.num=10),　
                     normMethod = "none",
                     sparsMethod = "none", #sparsMethod = "threshhold", thresh = 0.01
                     sampleSize = 540,  
                     dissFunc = "signed",
                     verbose = 2,
                     seed = 123456)

net2$assoMat1 %>% as.data.frame() -> col #fdr not related, pure Spearman correlation 

df <- data.frame(
  Correlation = c("vOTU160 & Ruminococcus bromii", "vOTU161 & Ruminococcus bromii","vOTU162 & Ruminococcus bromii","vOTU164 & Ruminococcus bromii",
                  "vOTU289 & Fusicatenibacter saccharivorans", "vOTU291 & Fusicatenibacter saccharivorans", "vOTU293 & Fusicatenibacter saccharivorans", "vOTU294 & Fusicatenibacter saccharivorans",
                  "vOTU870 & Bacteroides thetaiotaomicron", "vOTU913 & Bacteroides thetaiotaomicron"),
  Coefficient_JP = c(col[5,1], col[5,2], col[5,3], col[5,4],
                  col[13,9], col[13,10], col[13,11], col[13,12],
                  col[8,6], col[8,7])
)

df_fdr <- data.frame(
  Correlation = c("vOTU160 & Ruminococcus bromii", "vOTU161 & Ruminococcus bromii","vOTU162 & Ruminococcus bromii","vOTU164 & Ruminococcus bromii",
                  "vOTU289 & Fusicatenibacter saccharivorans", "vOTU291 & Fusicatenibacter saccharivorans", "vOTU293 & Fusicatenibacter saccharivorans", "vOTU294 & Fusicatenibacter saccharivorans",
                  "vOTU870 & Bacteroides thetaiotaomicron", "vOTU913 & Bacteroides thetaiotaomicron"),
  Coefficient_JP = c(col_fdr[5,1], col_fdr[5,2], col_fdr[5,3], col_fdr[5,4],
                  col_fdr[13,9], col_fdr[13,10], col_fdr[13,11], col_fdr[13,12],
                  col_fdr[8,6], col_fdr[8,7])
)

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")
p_JP =pheatmap(as.matrix(col), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(col_fdr >0 | col_fdr<0,"*", ""), nrow(col_fdr)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.5,0,1), c("navy", "white", "firebrick3")), name = "Spearman rho", heatmap_legend_param = list(color_bar = "continuous"), column_title = "Correlations of UC-related multibiome (Japanese 4D cohort)", column_title_gp = gpar(fontsize = 10, fontface = "bold"))
draw(p_JP, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig)) 

#External data1
SP_input2_bacteriome %>% filter(abs(SP_input2_bacteriome$"Coefficient (UC JP)") > 1 & SP_input2_bacteriome$"Q-value (UC JP)" < 0.1) -> SP_input2_bacteriome2 #For UC
SP_input2_bacteriome2 %>% dplyr::select(feature) -> Bacteriome_No 
#Phage list 
#vOTU_Heat_qval_fig from virome analysis including vOTU with abs coef >1 and FDR<0.1 in UC or CD
vOTU_Heat_qval_fig %>% filter(abs(vOTU_Heat_qval_fig$"Coefficient (UC JP)") > 1 & vOTU_Heat_qval_fig$"Q-value (UC JP)" < 0.1) -> vOTU_Heat_qval_fig2 #For UC
rownames_to_column(vOTU_Heat_qval_fig2, "No") -> vOTU_Heat_qval_fig2
vOTU_Heat_qval_fig2 %>% dplyr::select("No", "feature") -> Virome_No 

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
Ref <- read.csv("Reference2.csv", header = TRUE, na.strings = c(NA, ''))
left_join(Virome_No, Ref, by = "No") ->Virome_No
Virome_No <- Virome_No[, -which (colnames(Virome_No) %in% c("No"))] 
rename(.data=Virome_No, "No" = "vOTU_host") -> Virome_No

#Fungi list
#Fungi_map from mycobiome analysis including fungi with P<0.05 in UC or CD
Fungi_map %>% filter(Fungi_map $"P-value (UC JP)"< 0.05) -> Fungi_map2 #For UC
Fungi_map2 %>% dplyr::select("feature") -> Fungi_No

#External dataset
X<- "Franzosa_2018_IBD" #Select data, "He_2017_Crohn", "ES_MH, MetaHIT"
Y <- "USA" #Select location #"Netherlands", "Spain", "China"

#Bacteriome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/mOTU3_original_data")
SP_WO3<- read.csv("Bork_group.mOUTs3.csv", header = TRUE, na.strings = c(NA, ''))
rename(.data=SP_WO3, "sample_id" = "X",) -> SP_WO3
rename(.data=SP_WO3, "Unassigned species" = "X.1",) -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")

#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
prop.table(data.matrix(US2_SP), 1) -> US2_SP #Bacteriome data is not normalized. This script is needed.
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
US2_SP %>% mutate(ID = str_sub(US2_SP$feature, start = -6, end = -2)) -> US2_SP
Bacteriome_No %>% mutate(ID = str_sub(Bacteriome_No$feature, start = -6, end = -2)) -> Bacteriome_No
left_join(Bacteriome_No, US2_SP, by = "ID") -> US2_SP
column_to_rownames(US2_SP, "feature.x") -> US2_SP
Bacteriome <- US2_SP[, -which (colnames(US2_SP) %in% c("Biome", "ID", "feature.y"))] 
Bacteriome <- as.data.frame(t(Bacteriome))

#Virome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
SP_WO3<- read.csv("phage.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = TRUE)
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")
#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
#prop.table(data.matrix(US2_SP), 1) -> US2_SP (data is RPKM)
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
US2_SP$feature <- str_replace_all(US2_SP$feature, pattern="[:]", replacement=".") #Adjust the phage name before data integration
left_join(Virome_No, US2_SP, by = "feature") -> US2_SP #Please confirm all data is included
column_to_rownames(US2_SP, "No") -> US2_SP
Virome <- US2_SP[, -which (colnames(US2_SP) %in% c("feature"))] 
Virome <- as.data.frame(t(Virome))

#Mycobiome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Fungi")
SP_WO3<- read.csv("eukaryote.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = TRUE)
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")
#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
#prop.table(data.matrix(US2_SP), 1) -> US2_SP (data is RPKM)
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
Fungi_No$feature <- str_replace_all(Fungi_No$feature, pattern="[.]", replacement=" ") 
left_join(Fungi_No, US2_SP, by = "feature") -> US2_SP #Please confirm all data is included
column_to_rownames(US2_SP, "feature") -> Mycobiome
Mycobiome <- as.data.frame(t(Mycobiome))

#Network analysis
rownames_to_column(Bacteriome, "ID") -> Bacteriome
rownames_to_column(Virome, "ID") -> Virome
rownames_to_column(Mycobiome, "ID") -> Mycobiome
inner_join(inner_join(Bacteriome, Virome, by = "ID"), Mycobiome, by = "ID") -> SparData
column_to_rownames(.data = SparData, "ID") -> SparData

#Select important strains for the validation
SparData %>% select("vOTU160(Ruminococcus, Eubacterium)", "vOTU161(Ruminococcus, Eubacterium)", "vOTU162(Ruminococcus, Eubacterium)", "vOTU164(Ruminococcus, Eubacterium)", "Ruminococcus bromii [ID:00853]",
                    "vOTU870(Bacteroides, Alistipes)", "vOTU913(Bacteroides)", "Bacteroides thetaiotaomicron [ID:01657]",
                    "vOTU289(Blautia, Fusicatenibacter)", "vOTU291(Blautia, Fusicatenibacter)", "vOTU293(Blautia, Fusicatenibacter)", "vOTU294(Blautia, Fusicatenibacter)", "Fusicatenibacter saccharivorans [ID:02680]",
                    "vOTU497(Clostridium, Anaerotruncus, Ruminococcus, Blautia, Dorea, Coprococcus)", "Dorea longicatena [ID:03693]", "Ruminococcaceae species incertae sedis [ID:12239]", "Blautia massiliensis [ID:03342]", 
                    "Clostridiales species incertae sedis [ID:12476]", "Clostridiales species incertae sedis [ID:12554]",
                    "vOTU331(Faecalibacterium)", "Faecalibacterium prausnitzii [ID:06112]") -> SparData

net <- netConstruct(SparData, #relative abundance
                    filtTax = "none",
                    filtSamp = "none", 
                    filtSampPar = "none", 
                    measure = "spearman",
                    thresh = 0.2,
                    measurePar = list(nlambda=10, 
                                      rep.num=10),　
                    normMethod = "none",
                    sparsMethod = "t-test", #sparsMethod = "threshhold", thresh = 0.01
                    adjust = "fdr",
                    alpha = 0.1, 
                    lfdrThresh = 0.1,
                    sampleSize = 34,  
                    dissFunc = "signed",
                    verbose = 2,
                    seed = 123456)

net$assoMat1 %>% as.data.frame() -> col_fdr #if fdr<0.1, coef included, otherwise 0

net2 <- netConstruct(SparData, #relative abundance
                     filtTax = "none",
                     filtSamp = "none", 
                     filtSampPar = "none", 
                     measure = "spearman",
                     measurePar = list(nlambda=10, 
                                       rep.num=10),　
                     normMethod = "none",
                     sparsMethod = "none", #sparsMethod = "threshhold", thresh = 0.01
                     sampleSize = 34,  
                     dissFunc = "signed",
                     verbose = 2,
                     seed = 123456)

net2$assoMat1 %>% as.data.frame() -> col #fdr not related, pure Spearman correlation 

df_US <- data.frame(
  Correlation = c("vOTU160 & Ruminococcus bromii", "vOTU161 & Ruminococcus bromii","vOTU162 & Ruminococcus bromii","vOTU164 & Ruminococcus bromii",
                  "vOTU289 & Fusicatenibacter saccharivorans", "vOTU291 & Fusicatenibacter saccharivorans", "vOTU293 & Fusicatenibacter saccharivorans", "vOTU294 & Fusicatenibacter saccharivorans",
                  "vOTU870 & Bacteroides thetaiotaomicron", "vOTU913 & Bacteroides thetaiotaomicron"),
  Coefficient_Franzosa_US = c(col[5,1], col[5,2], col[5,3], col[5,4],
                  col[13,9], col[13,10], col[13,11], col[13,12],
                  col[8,6], col[8,7])
)

df_fdr_US <- data.frame(
  Correlation = c("vOTU160 & Ruminococcus bromii", "vOTU161 & Ruminococcus bromii","vOTU162 & Ruminococcus bromii","vOTU164 & Ruminococcus bromii",
                  "vOTU289 & Fusicatenibacter saccharivorans", "vOTU291 & Fusicatenibacter saccharivorans", "vOTU293 & Fusicatenibacter saccharivorans", "vOTU294 & Fusicatenibacter saccharivorans",
                  "vOTU870 & Bacteroides thetaiotaomicron", "vOTU913 & Bacteroides thetaiotaomicron"),
  Coefficient_Franzosa_US = c(col_fdr[5,1], col_fdr[5,2], col_fdr[5,3], col_fdr[5,4],
                  col_fdr[13,9], col_fdr[13,10], col_fdr[13,11], col_fdr[13,12],
                  col_fdr[8,6], col_fdr[8,7])
)

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")
p_US =pheatmap(as.matrix(col), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(col_fdr >0 | col_fdr<0,"*", ""), nrow(col_fdr)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.5,0,1), c("navy", "white", "firebrick3")), name = "Spearman rho", heatmap_legend_param = list(color_bar = "continuous"), column_title = "Correlations of UC-related multibiome (US cohort)", column_title_gp = gpar(fontsize = 10, fontface = "bold"))
draw(p_US, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig)) 

#External data2
SP_input2_bacteriome %>% filter(abs(SP_input2_bacteriome$"Coefficient (UC JP)") > 1 & SP_input2_bacteriome$"Q-value (UC JP)" < 0.1) -> SP_input2_bacteriome2 #For UC
SP_input2_bacteriome2 %>% dplyr::select(feature) -> Bacteriome_No 
#Phage list 
#vOTU_Heat_qval_fig from virome analysis including vOTU with abs coef >1 and FDR<0.1 in UC or CD
vOTU_Heat_qval_fig %>% filter(abs(vOTU_Heat_qval_fig$"Coefficient (UC JP)") > 1 & vOTU_Heat_qval_fig$"Q-value (UC JP)" < 0.1) -> vOTU_Heat_qval_fig2 #For UC
rownames_to_column(vOTU_Heat_qval_fig2, "No") -> vOTU_Heat_qval_fig2
vOTU_Heat_qval_fig2 %>% dplyr::select("No", "feature") -> Virome_No 

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
Ref <- read.csv("Reference2.csv", header = TRUE, na.strings = c(NA, ''))
left_join(Virome_No, Ref, by = "No") ->Virome_No
Virome_No <- Virome_No[, -which (colnames(Virome_No) %in% c("No"))] 
rename(.data=Virome_No, "No" = "vOTU_host") -> Virome_No

#Fungi list
#Fungi_map from mycobiome analysis including fungi with P<0.05 in UC or CD
Fungi_map %>% filter(Fungi_map $"P-value (UC JP)"< 0.05) -> Fungi_map2 #For UC
Fungi_map2 %>% dplyr::select("feature") -> Fungi_No

#External dataset
X<- "Franzosa_2018_IBD" #Select data, "He_2017_Crohn", "ES_MH, MetaHIT"
Y <- "Netherlands" #Select location #"Netherlands", "Spain", "China"

#Bacteriome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/mOTU3_original_data")
SP_WO3<- read.csv("Bork_group.mOUTs3.csv", header = TRUE, na.strings = c(NA, ''))
rename(.data=SP_WO3, "sample_id" = "X",) -> SP_WO3
rename(.data=SP_WO3, "Unassigned species" = "X.1",) -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")

#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
prop.table(data.matrix(US2_SP), 1) -> US2_SP #Bacteriome data is not normalized. This script is needed.
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
US2_SP %>% mutate(ID = str_sub(US2_SP$feature, start = -6, end = -2)) -> US2_SP
Bacteriome_No %>% mutate(ID = str_sub(Bacteriome_No$feature, start = -6, end = -2)) -> Bacteriome_No
left_join(Bacteriome_No, US2_SP, by = "ID") -> US2_SP
column_to_rownames(US2_SP, "feature.x") -> US2_SP
Bacteriome <- US2_SP[, -which (colnames(US2_SP) %in% c("Biome", "ID", "feature.y"))] 
Bacteriome <- as.data.frame(t(Bacteriome))

#Virome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
SP_WO3<- read.csv("phage.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = TRUE)
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")
#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
#prop.table(data.matrix(US2_SP), 1) -> US2_SP (data is RPKM)
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
US2_SP$feature <- str_replace_all(US2_SP$feature, pattern="[:]", replacement=".") #Adjust the phage name before data integration
left_join(Virome_No, US2_SP, by = "feature") -> US2_SP #Please confirm all data is included
column_to_rownames(US2_SP, "No") -> US2_SP
Virome <- US2_SP[, -which (colnames(US2_SP) %in% c("feature"))] 
Virome <- as.data.frame(t(Virome))

#Mycobiome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Fungi")
SP_WO3<- read.csv("eukaryote.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = TRUE)
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")
#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
#prop.table(data.matrix(US2_SP), 1) -> US2_SP (data is RPKM)
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
Fungi_No$feature <- str_replace_all(Fungi_No$feature, pattern="[.]", replacement=" ") 
left_join(Fungi_No, US2_SP, by = "feature") -> US2_SP #Please confirm all data is included
column_to_rownames(US2_SP, "feature") -> Mycobiome
Mycobiome <- as.data.frame(t(Mycobiome))

#Network analysis
rownames_to_column(Bacteriome, "ID") -> Bacteriome
rownames_to_column(Virome, "ID") -> Virome
rownames_to_column(Mycobiome, "ID") -> Mycobiome
inner_join(inner_join(Bacteriome, Virome, by = "ID"), Mycobiome, by = "ID") -> SparData
column_to_rownames(.data = SparData, "ID") -> SparData

#Select important strains for the validation
SparData %>% select("vOTU160(Ruminococcus, Eubacterium)", "vOTU161(Ruminococcus, Eubacterium)", "vOTU162(Ruminococcus, Eubacterium)", "vOTU164(Ruminococcus, Eubacterium)", "Ruminococcus bromii [ID:00853]",
                    "vOTU870(Bacteroides, Alistipes)", "vOTU913(Bacteroides)", "Bacteroides thetaiotaomicron [ID:01657]",
                    "vOTU289(Blautia, Fusicatenibacter)", "vOTU291(Blautia, Fusicatenibacter)", "vOTU293(Blautia, Fusicatenibacter)", "vOTU294(Blautia, Fusicatenibacter)", "Fusicatenibacter saccharivorans [ID:02680]",
                    "vOTU497(Clostridium, Anaerotruncus, Ruminococcus, Blautia, Dorea, Coprococcus)", "Dorea longicatena [ID:03693]", "Ruminococcaceae species incertae sedis [ID:12239]", "Blautia massiliensis [ID:03342]", 
                    "Clostridiales species incertae sedis [ID:12476]", "Clostridiales species incertae sedis [ID:12554]",
                    "vOTU331(Faecalibacterium)", "Faecalibacterium prausnitzii [ID:06112]") -> SparData

net <- netConstruct(SparData, #relative abundance
                    filtTax = "none",
                    filtSamp = "none", 
                    filtSampPar = "none", 
                    measure = "spearman",
                    thresh = 0.2,
                    measurePar = list(nlambda=10, 
                                      rep.num=10),　
                    normMethod = "none",
                    sparsMethod = "t-test", #sparsMethod = "threshhold", thresh = 0.01
                    adjust = "fdr",
                    alpha = 0.1, 
                    lfdrThresh = 0.1,
                    sampleSize = 22,  
                    dissFunc = "signed",
                    verbose = 2,
                    seed = 123456)

net$assoMat1 %>% as.data.frame() -> col_fdr #if fdr<0.1, coef included, otherwise 0

net2 <- netConstruct(SparData, #relative abundance
                     filtTax = "none",
                     filtSamp = "none", 
                     filtSampPar = "none", 
                     measure = "spearman",
                     measurePar = list(nlambda=10, 
                                       rep.num=10),　
                     normMethod = "none",
                     sparsMethod = "none", #sparsMethod = "threshhold", thresh = 0.01
                     sampleSize = 22,  
                     dissFunc = "signed",
                     verbose = 2,
                     seed = 123456)

net2$assoMat1 %>% as.data.frame() -> col #fdr not related, pure Spearman correlation 

df_Netherlands <- data.frame(
  Correlation = c("vOTU160 & Ruminococcus bromii", "vOTU161 & Ruminococcus bromii","vOTU162 & Ruminococcus bromii","vOTU164 & Ruminococcus bromii",
                  "vOTU289 & Fusicatenibacter saccharivorans", "vOTU291 & Fusicatenibacter saccharivorans", "vOTU293 & Fusicatenibacter saccharivorans", "vOTU294 & Fusicatenibacter saccharivorans",
                  "vOTU870 & Bacteroides thetaiotaomicron", "vOTU913 & Bacteroides thetaiotaomicron"),
  Coefficient_Franzosa_Netherlands = c(col[5,1], col[5,2], col[5,3], col[5,4],
                              col[13,9], col[13,10], col[13,11], col[13,12],
                              col[8,6], col[8,7])
)

df_fdr_Netherlands <- data.frame(
  Correlation = c("vOTU160 & Ruminococcus bromii", "vOTU161 & Ruminococcus bromii","vOTU162 & Ruminococcus bromii","vOTU164 & Ruminococcus bromii",
                  "vOTU289 & Fusicatenibacter saccharivorans", "vOTU291 & Fusicatenibacter saccharivorans", "vOTU293 & Fusicatenibacter saccharivorans", "vOTU294 & Fusicatenibacter saccharivorans",
                  "vOTU870 & Bacteroides thetaiotaomicron", "vOTU913 & Bacteroides thetaiotaomicron"),
  Coefficient_Franzosa_Netherlands = c(col_fdr[5,1], col_fdr[5,2], col_fdr[5,3], col_fdr[5,4],
                              col_fdr[13,9], col_fdr[13,10], col_fdr[13,11], col_fdr[13,12],
                              col_fdr[8,6], col_fdr[8,7])
)

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")
p_NL=pheatmap(as.matrix(col), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(col_fdr >0 | col_fdr<0,"*", ""), nrow(col_fdr)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.5,0,1), c("navy", "white", "firebrick3")), name = "Spearman rho", heatmap_legend_param = list(color_bar = "continuous"), column_title = "Correlations of UC-related multibiome (NL cohort)", column_title_gp = gpar(fontsize = 10, fontface = "bold"))
draw(p_NL, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig)) 

#External dataset3
SP_input2_bacteriome %>% filter(abs(SP_input2_bacteriome$"Coefficient (UC JP)") > 1 & SP_input2_bacteriome$"Q-value (UC JP)" < 0.1) -> SP_input2_bacteriome2 #For UC
SP_input2_bacteriome2 %>% dplyr::select(feature) -> Bacteriome_No 
#Phage list 
#vOTU_Heat_qval_fig from virome analysis including vOTU with abs coef >1 and FDR<0.1 in UC or CD
vOTU_Heat_qval_fig %>% filter(abs(vOTU_Heat_qval_fig$"Coefficient (UC JP)") > 1 & vOTU_Heat_qval_fig$"Q-value (UC JP)" < 0.1) -> vOTU_Heat_qval_fig2 #For UC
rownames_to_column(vOTU_Heat_qval_fig2, "No") -> vOTU_Heat_qval_fig2
vOTU_Heat_qval_fig2 %>% dplyr::select("No", "feature") -> Virome_No 

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
Ref <- read.csv("Reference2.csv", header = TRUE, na.strings = c(NA, ''))
left_join(Virome_No, Ref, by = "No") ->Virome_No
Virome_No <- Virome_No[, -which (colnames(Virome_No) %in% c("No"))] 
rename(.data=Virome_No, "No" = "vOTU_host") -> Virome_No

#Fungi list
#Fungi_map from mycobiome analysis including fungi with P<0.05 in UC or CD
Fungi_map %>% filter(Fungi_map $"P-value (UC JP)"< 0.05) -> Fungi_map2 #For UC
Fungi_map2 %>% dplyr::select("feature") -> Fungi_No

X<- "ES_MH, MetaHIT" #Select data, "He_2017_Crohn", "ES_MH, MetaHIT"
Y <- "Spain" #Select location #"Netherlands", "Spain", "China"

#Bacteriome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/mOTU3_original_data")
SP_WO3<- read.csv("Bork_group.mOUTs3.csv", header = TRUE, na.strings = c(NA, ''))
rename(.data=SP_WO3, "sample_id" = "X",) -> SP_WO3
rename(.data=SP_WO3, "Unassigned species" = "X.1",) -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")

#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
prop.table(data.matrix(US2_SP), 1) -> US2_SP #Bacteriome data is not normalized. This script is needed.
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
US2_SP %>% mutate(ID = str_sub(US2_SP$feature, start = -6, end = -2)) -> US2_SP
Bacteriome_No %>% mutate(ID = str_sub(Bacteriome_No$feature, start = -6, end = -2)) -> Bacteriome_No
left_join(Bacteriome_No, US2_SP, by = "ID") -> US2_SP
column_to_rownames(US2_SP, "feature.x") -> US2_SP
Bacteriome <- US2_SP[, -which (colnames(US2_SP) %in% c("Biome", "ID", "feature.y"))] 
Bacteriome <- as.data.frame(t(Bacteriome))

#Virome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
SP_WO3<- read.csv("phage.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = TRUE)
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")
#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
#prop.table(data.matrix(US2_SP), 1) -> US2_SP (data is RPKM)
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
US2_SP$feature <- str_replace_all(US2_SP$feature, pattern="[:]", replacement=".") #Adjust the phage name before data integration
left_join(Virome_No, US2_SP, by = "feature") -> US2_SP #Please confirm all data is included
column_to_rownames(US2_SP, "No") -> US2_SP
Virome <- US2_SP[, -which (colnames(US2_SP) %in% c("feature"))] 
Virome <- as.data.frame(t(Virome))

#Mycobiome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Fungi")
SP_WO3<- read.csv("eukaryote.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = TRUE)
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")
#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
#prop.table(data.matrix(US2_SP), 1) -> US2_SP (data is RPKM)
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
Fungi_No$feature <- str_replace_all(Fungi_No$feature, pattern="[.]", replacement=" ") 
left_join(Fungi_No, US2_SP, by = "feature") -> US2_SP #Please confirm all data is included
column_to_rownames(US2_SP, "feature") -> Mycobiome
Mycobiome <- as.data.frame(t(Mycobiome))

#Network analysis
rownames_to_column(Bacteriome, "ID") -> Bacteriome
rownames_to_column(Virome, "ID") -> Virome
rownames_to_column(Mycobiome, "ID") -> Mycobiome
inner_join(inner_join(Bacteriome, Virome, by = "ID"), Mycobiome, by = "ID") -> SparData
column_to_rownames(.data = SparData, "ID") -> SparData

SparData <- na.omit(SparData) #Spain data has deficient data of fungi and phage

#Select important strains for the validation
SparData %>% select("vOTU160(Ruminococcus, Eubacterium)", "vOTU161(Ruminococcus, Eubacterium)", "vOTU162(Ruminococcus, Eubacterium)", "vOTU164(Ruminococcus, Eubacterium)", "Ruminococcus bromii [ID:00853]",
                    "vOTU870(Bacteroides, Alistipes)", "vOTU913(Bacteroides)", "Bacteroides thetaiotaomicron [ID:01657]",
                    "vOTU289(Blautia, Fusicatenibacter)", "vOTU291(Blautia, Fusicatenibacter)", "vOTU293(Blautia, Fusicatenibacter)", "vOTU294(Blautia, Fusicatenibacter)", "Fusicatenibacter saccharivorans [ID:02680]",
                    "vOTU497(Clostridium, Anaerotruncus, Ruminococcus, Blautia, Dorea, Coprococcus)", "Dorea longicatena [ID:03693]", "Ruminococcaceae species incertae sedis [ID:12239]", "Blautia massiliensis [ID:03342]", 
                    "Clostridiales species incertae sedis [ID:12476]", "Clostridiales species incertae sedis [ID:12554]",
                    "vOTU331(Faecalibacterium)", "Faecalibacterium prausnitzii [ID:06112]") -> SparData

net <- netConstruct(SparData, #relative abundance
                    filtTax = "none",
                    filtSamp = "none", 
                    filtSampPar = "none", 
                    measure = "spearman",
                    thresh = 0.2,
                    measurePar = list(nlambda=10, 
                                      rep.num=10),　
                    normMethod = "none",
                    sparsMethod = "t-test", #sparsMethod = "threshhold", thresh = 0.01
                    adjust = "fdr",
                    alpha = 0.1, 
                    lfdrThresh = 0.1,
                    sampleSize = 62,
                    dissFunc = "signed",
                    verbose = 2,
                    seed = 123456)

net$assoMat1 %>% as.data.frame() -> col_fdr #if fdr<0.1, coef included, otherwise 0

net2 <- netConstruct(SparData, #relative abundance
                     filtTax = "none",
                     filtSamp = "none", 
                     filtSampPar = "none", 
                     measure = "spearman",
                     measurePar = list(nlambda=10, 
                                       rep.num=10),　
                     normMethod = "none",
                     sparsMethod = "none", #sparsMethod = "threshhold", thresh = 0.01
                     sampleSize = 62,  #要確認
                     dissFunc = "signed",
                     verbose = 2,
                     seed = 123456)

net2$assoMat1 %>% as.data.frame() -> col #fdr not related, pure Spearman correlation 

df_Spain <- data.frame(
  Correlation = c("vOTU160 & Ruminococcus bromii", "vOTU161 & Ruminococcus bromii","vOTU162 & Ruminococcus bromii","vOTU164 & Ruminococcus bromii",
                  "vOTU289 & Fusicatenibacter saccharivorans", "vOTU291 & Fusicatenibacter saccharivorans", "vOTU293 & Fusicatenibacter saccharivorans", "vOTU294 & Fusicatenibacter saccharivorans",
                  "vOTU870 & Bacteroides thetaiotaomicron", "vOTU913 & Bacteroides thetaiotaomicron"),
  Coefficient_Nielsen_Spain = c(col[5,1], col[5,2], col[5,3], col[5,4],
                                       col[13,9], col[13,10], col[13,11], col[13,12],
                                       col[8,6], col[8,7])
)

df_fdr_Spain <- data.frame(
  Correlation = c("vOTU160 & Ruminococcus bromii", "vOTU161 & Ruminococcus bromii","vOTU162 & Ruminococcus bromii","vOTU164 & Ruminococcus bromii",
                  "vOTU289 & Fusicatenibacter saccharivorans", "vOTU291 & Fusicatenibacter saccharivorans", "vOTU293 & Fusicatenibacter saccharivorans", "vOTU294 & Fusicatenibacter saccharivorans",
                  "vOTU870 & Bacteroides thetaiotaomicron", "vOTU913 & Bacteroides thetaiotaomicron"),
  Coefficient_Nielsen_Spain = c(col_fdr[5,1], col_fdr[5,2], col_fdr[5,3], col_fdr[5,4],
                                       col_fdr[13,9], col_fdr[13,10], col_fdr[13,11], col_fdr[13,12],
                                       col_fdr[8,6], col_fdr[8,7])
)

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")
p_ES=pheatmap(as.matrix(col), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(col_fdr >0 | col_fdr<0,"*", ""), nrow(col_fdr)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.5,0,1), c("navy", "white", "firebrick3")), name = "Spearman rho", heatmap_legend_param = list(color_bar = "continuous"), column_title = "Correlations of UC-related multibiome (ES cohort)", column_title_gp = gpar(fontsize = 10, fontface = "bold"))
draw(p_ES, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig)) 


full_join(full_join(full_join(df, df_US, by = "Correlation"), df_Spain, by = "Correlation"), df_Netherlands, by = "Correlation") -> all_df
full_join(full_join(full_join(df_fdr, df_fdr_US, by = "Correlation"), df_fdr_Spain, by = "Correlation"), df_fdr_Netherlands, by = "Correlation") -> all_df_fdr

column_to_rownames(all_df, "Correlation") -> all_df
column_to_rownames(all_df_fdr, "Correlation") -> all_df_fdr

rename(.data = all_df, "Japanese 4D cohort" = "Coefficient_JP") -> all_df
rename(.data = all_df, "US cohort" = "Coefficient_Franzosa_US") -> all_df
rename(.data = all_df, "NL cohort" = "Coefficient_Franzosa_Netherlands") -> all_df
rename(.data = all_df, "ES cohort" = "Coefficient_Nielsen_Spain") -> all_df

rename(.data = all_df_fdr, "Japanese 4D cohort" = "Coefficient_JP") -> all_df_fdr
rename(.data = all_df_fdr, "US cohort" = "Coefficient_Franzosa_US") -> all_df_fdr
rename(.data = all_df_fdr, "NL cohort" = "Coefficient_Franzosa_Netherlands") -> all_df_fdr
rename(.data = all_df_fdr, "ES cohort" = "Coefficient_Nielsen_Spain") -> all_df_fdr

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")
p =pheatmap(as.matrix(all_df), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(all_df_fdr != 0,"*", ""), nrow(all_df_fdr)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.5,0,1), c("navy", "white", "firebrick3")), name = "Spearman rho", heatmap_legend_param = list(color_bar = "continuous"), column_title = "Correlations between UC-related species", column_title_gp = gpar(fontsize = 10, fontface = "bold"))
draw(p, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig)) 


#CD analysis
#Bacteria list
#SP_input2_bacteriome from bacteriome analysis including bacterial species with abs coef >1 and FDR<0.1 in UC or CD
rownames_to_column(SP_input2_bacteriome, "feature") -> SP_input2_bacteriome 

#Repeat from here
SP_input2_bacteriome %>% filter(abs(SP_input2_bacteriome$"Coefficient (CD JP)") > 1 & SP_input2_bacteriome$"Q-value (CD JP)" < 0.1) -> SP_input2_bacteriome2 #For CD 
SP_input2_bacteriome2 %>% dplyr::select(feature) -> Bacteriome_No 
#Phage list 
#vOTU_Heat_qval_fig from virome analysis including vOTU with abs coef >1 and FDR<0.1 in UC or CD
vOTU_Heat_qval_fig %>% filter(abs(vOTU_Heat_qval_fig $"Coefficient (CD JP)") > 1 & vOTU_Heat_qval_fig$"Q-value (CD JP)" < 0.1) -> vOTU_Heat_qval_fig2 #For CD
rownames_to_column(vOTU_Heat_qval_fig2, "No") -> vOTU_Heat_qval_fig2
vOTU_Heat_qval_fig2 %>% dplyr::select("No", "feature") -> Virome_No 

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
Ref <- read.csv("Reference2.csv", header = TRUE, na.strings = c(NA, ''))
left_join(Virome_No, Ref, by = "No") ->Virome_No
Virome_No <- Virome_No[, -which (colnames(Virome_No) %in% c("No"))] 
rename(.data=Virome_No, "No" = "vOTU_host") -> Virome_No

#Fungi list 
#Fungi_map from mycobiome analysis including fungi with P<0.05 in UC or CD
Fungi_map %>% filter(Fungi_map$"P-value (CD JP)"< 0.05) -> Fungi_map2 #For CD
Fungi_map2 %>% dplyr::select("feature") -> Fungi_No

# Preparation for bacterial species abundance (Japanese 4D cohort, please skip this script to external data if JP 4D is not used)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/SP_mOTU3(METAFは番号に修正済みIBD症状薬剤解析)")
#CDvsControl
#Bacteriome abundance (JP)
SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1, check.names = FALSE) #abundance data
rownames_to_column(SP_input, "ID") ->SP_input
SP_input %>% filter(SP_input$IBD ==2) -> SP_input #For Cont
SP_input %>% mutate(ID = as.numeric(SP_input$ID)) -> SP_input
NO <- SP_input%>%dplyr::select(ID) #CD vs Cont is selected.
rename(.data=NO, "NO" = "ID") -> NO
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("METAF", "IBD", "UC", "CD"))] 
column_to_rownames(SP_input2, "ID") -> SP_input2
SP_input2 <- as.data.frame(t(SP_input2)) 
rownames_to_column(SP_input2, "feature") -> SP_input2 
left_join(Bacteriome_No, SP_input2, by = "feature") -> SP_input2
column_to_rownames(SP_input2, "feature") -> Bacteriome
Bacteriome <- as.data.frame(t(Bacteriome))

#vOTU abundance (JP)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
vOTU_input <- read.csv("vOTU.csv", header = TRUE, na.strings = c(NA, '')) #abundance data
left_join(NO, vOTU_input, by = "NO") -> vOTU_input 
vOTU_input <- vOTU_input[, -which (colnames(vOTU_input) %in% c("IBD_HC"))] 
column_to_rownames(vOTU_input, "NO") -> vOTU_input
vOTU_input <- as.data.frame(t(vOTU_input)) 
rownames_to_column(vOTU_input, "feature") -> vOTU_input 
left_join(Virome_No, vOTU_input, by = "feature") -> vOTU_input2
vOTU_input2 <- vOTU_input2[, -which (colnames(vOTU_input2) %in% c("feature"))] #Remove non-essential variables。
column_to_rownames(vOTU_input2, "No") -> Virome
Virome <- as.data.frame(t(Virome))

#Fungi abundance (JP)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Fungi")
F_input <- read.csv("feces_230925.eukaryote.csv", header = TRUE, na.strings = c(NA, ''), row.names=1, check.names = FALSE) 
as.data.frame(t(F_input)) -> F_input
rownames_to_column(F_input, var = "NO") -> F_input
F_input %>% mutate(NO = as.numeric(F_input$NO)) -> F_input
left_join(NO, F_input, by = "NO") -> F_input #Select CD vs Cont
column_to_rownames(F_input, "NO") -> F_input
F_input <- as.data.frame(t(F_input)) 
rownames_to_column(F_input, "feature") -> F_input
Fungi_No$feature <- str_replace_all(Fungi_No$feature, pattern="[.]", replacement=" ")
left_join(Fungi_No, F_input, "feature") -> F_input2
column_to_rownames(F_input2, "feature") -> Mycobiome
Mycobiome <- as.data.frame(t(Mycobiome))

#Network analysis
rownames_to_column(Bacteriome, "ID") -> Bacteriome
rownames_to_column(Virome, "ID") -> Virome
rownames_to_column(Mycobiome, "ID") -> Mycobiome
inner_join(inner_join(Bacteriome, Virome, by = "ID"), Mycobiome, by = "ID") -> SparData
column_to_rownames(.data = SparData, "ID") -> SparData

#Select important strains for validation
SparData %>% select("vOTU913(Bacteroides)", "vOTU914(Bacteroides)", "Bacteroides thetaiotaomicron [ID:01657]", 
                    "vOTU77(Escherichia)", "Escherichia coli [ID:00095]", 
                    "vOTU758(Anaerostipes)", "Anaerostipes hadrus [ID:00857]", "Anaerostipes hadrus [ID:00856]",
                    "vOTU689(Streptococcus)", "Streptococcus sp. [ID:00902]", 
                    "Saccharomyces cerevisiae", "Eubacterium ventriosum [ID:11151]", "Clostridiales species incertae sedis [ID:12288]",
                    "Clostridium sp. [ID:03680]", "Clostridiales species incertae sedis [ID:12563]", "Blautia obeum [ID:05139]", "Clostridiales species incertae sedis [ID:12284]",
                    "vOTU388(Clostridium, Faecalibacterium)", "Faecalibacterium prausnitzii [ID:06110]",
                    "[Clostridium] leptum [ID:03688]", "Clostridiales sp. [ID:03658]",  "Clostridium phoceensis [ID:08865]", "Faecalibacterium prausnitzii [ID:06108]", 
                    "vOTU1328(Unclassified Lachnospiraceae, Unclassified Clostridiales, Flavonifractor)", "Clostridiales species incertae sedis [ID:12473]", 
                    "vOTU756(Lachnoclostridium, Ruminococcus)", "Ruminococcaceae species incertae sedis [ID:12239]",
                    "vOTU497(Clostridium, Anaerotruncus, Ruminococcus, Blautia, Dorea, Coprococcus)", "Dorea longicatena [ID:03693]", "Blautia massiliensis [ID:03342]",
                    "Clostridiales species incertae sedis [ID:12476]","Dorea formicigenerans [ID:03668]", "Clostridiales species incertae sedis [ID:12554]", "Coprococcus comes [ID:03690]", 
                    "vOTU533(Roseburia, Coprococcus, Unclassified Firmicutes, Clostridium, Ruminococcus)", "[Ruminococcus] torques [ID:03659]", "Clostridiales sp. [ID:03658]", "uncultured Flavonifractor sp. [ID:07315]") -> SparData

net <- netConstruct(SparData, #relative abundance
                    filtTax = "none",
                    filtSamp = "none", 
                    filtSampPar = "none", 
                    measure = "spearman",
                    thresh = 0.2,
                    measurePar = list(nlambda=10, 
                                      rep.num=10),　
                    normMethod = "none",
                    sparsMethod = "t-test", #sparsMethod = "threshhold", thresh = 0.01
                    adjust = "fdr",
                    alpha = 0.1, 
                    lfdrThresh = 0.1,
                    sampleSize = 540,  
                    dissFunc = "signed",
                    verbose = 2,
                    seed = 123456)

net$assoMat1 %>% as.data.frame() -> col_fdr #if fdr<0.1, coef included, otherwise 0

net2 <- netConstruct(SparData, #relative abundance
                     filtTax = "none",
                     filtSamp = "none", 
                     filtSampPar = "none", 
                     measure = "spearman",
                     measurePar = list(nlambda=10, 
                                       rep.num=10),　
                     normMethod = "none",
                     sparsMethod = "none", #sparsMethod = "threshhold", thresh = 0.01
                     sampleSize = 540,  
                     dissFunc = "signed",
                     verbose = 2,
                     seed = 123456)

net2$assoMat1 %>% as.data.frame() -> col #fdr not related, pure Spearman correlation 

df <- data.frame(
  Correlation = c("vOTU913 & Bacteroides thetaiotaomicron", "vOTU914 & Bacteroides thetaiotaomicron",
                  "vOTU77 & Escherichia coli",
                  "vOTU758 & Anaerostipes hadrus[ID:00857]", "vOTU758 & Anaerostipes hadrus[ID:00856]", 
                  "vOTU689 & Streptococcus sp.",
                  "Saccharomyces cerevisiae & Eubacterium ventriosum", "Saccharomyces cerevisiae & Clostridiales species incertae sedis",
                  "Escherichia coli & Anaerostipes hadrus[ID:00857]", "Escherichia coli & Anaerostipes hadrus[ID:00856]","Escherichia coli & Blautia obeum", "Escherichia coli & Clostridium sp.","Escherichia coli & Clostridiales species incertae sedis[ID:12563]","Escherichia coli & Clostridiales species incertae sedis[ID:12284]"),
  Coefficient_JP = c(col[3,1], col[3,2], 
                     col[4,5],
                     col[6,7],col[6,8],
                     col[9,10],
                     col[11,12],col[11,13],
                     col[5,7], col[5,8], col[5,16], col[5,14], col[5,15],col[5,17]
                     )
)

df_fdr <- data.frame(
  Correlation = c("vOTU913 & Bacteroides thetaiotaomicron", "vOTU914 & Bacteroides thetaiotaomicron",
                  "vOTU77 & Escherichia coli",
                  "vOTU758 & Anaerostipes hadrus[ID:00857]", "vOTU758 & Anaerostipes hadrus[ID:00856]", 
                  "vOTU689 & Streptococcus sp.",
                  "Saccharomyces cerevisiae & Eubacterium ventriosum", "Saccharomyces cerevisiae & Clostridiales species incertae sedis",
                  "Escherichia coli & Anaerostipes hadrus[ID:00857]", "Escherichia coli & Anaerostipes hadrus[ID:00856]","Escherichia coli & Blautia obeum", "Escherichia coli & Clostridium sp.","Escherichia coli & Clostridiales species incertae sedis[ID:12563]","Escherichia coli & Clostridiales species incertae sedis[ID:12284]"),
  Coefficient_JP = c(col_fdr[3,1], col_fdr[3,2], 
                     col_fdr[4,5],
                     col_fdr[6,7],col_fdr[6,8],
                     col_fdr[9,10],
                     col_fdr[11,12],col_fdr[11,13],
                     col_fdr[5,7], col_fdr[5,8], col_fdr[5,16], col_fdr[5,14], col_fdr[5,15],col_fdr[5,17]
  )
)

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")
p_JPc=pheatmap(as.matrix(col), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(col_fdr >0 | col_fdr<0,"*", ""), nrow(col_fdr)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.5,0,1), c("navy", "white", "firebrick3")), name = "Spearman rho", heatmap_legend_param = list(color_bar = "continuous"), column_title = "Correlations of CD-related multibiome (Japanese 4D cohort)", column_title_gp = gpar(fontsize = 10, fontface = "bold"))
draw(p_JPc, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig)) 

#External data1
#Repeat from here
SP_input2_bacteriome %>% filter(abs(SP_input2_bacteriome$"Coefficient (CD JP)") > 1 & SP_input2_bacteriome$"Q-value (CD JP)" < 0.1) -> SP_input2_bacteriome2 #For CD 
SP_input2_bacteriome2 %>% dplyr::select(feature) -> Bacteriome_No 
#Phage list 
#vOTU_Heat_qval_fig from virome analysis including vOTU with abs coef >1 and FDR<0.1 in UC or CD
vOTU_Heat_qval_fig %>% filter(abs(vOTU_Heat_qval_fig $"Coefficient (CD JP)") > 1 & vOTU_Heat_qval_fig$"Q-value (CD JP)" < 0.1) -> vOTU_Heat_qval_fig2 #For CD
rownames_to_column(vOTU_Heat_qval_fig2, "No") -> vOTU_Heat_qval_fig2
vOTU_Heat_qval_fig2 %>% dplyr::select("No", "feature") -> Virome_No 

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
Ref <- read.csv("Reference2.csv", header = TRUE, na.strings = c(NA, ''))
left_join(Virome_No, Ref, by = "No") ->Virome_No
Virome_No <- Virome_No[, -which (colnames(Virome_No) %in% c("No"))] 
rename(.data=Virome_No, "No" = "vOTU_host") -> Virome_No

#Fungi list 
#Fungi_map from mycobiome analysis including fungi with P<0.05 in UC or CD
Fungi_map %>% filter(Fungi_map$"P-value (CD JP)"< 0.05) -> Fungi_map2 #For CD
Fungi_map2 %>% dplyr::select("feature") -> Fungi_No

#External dataset
X<- "Franzosa_2018_IBD" #Select data, "He_2017_Crohn", "ES_MH, MetaHIT"
Y <- "USA" #Select location #"Netherlands", "Spain", "China"

#Bacteriome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/mOTU3_original_data")
SP_WO3<- read.csv("Bork_group.mOUTs3.csv", header = TRUE, na.strings = c(NA, ''))
rename(.data=SP_WO3, "sample_id" = "X",) -> SP_WO3
rename(.data=SP_WO3, "Unassigned species" = "X.1",) -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")

#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
prop.table(data.matrix(US2_SP), 1) -> US2_SP #Bacteriome data is not normalized. This script is needed.
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
US2_SP %>% mutate(ID = str_sub(US2_SP$feature, start = -6, end = -2)) -> US2_SP
Bacteriome_No %>% mutate(ID = str_sub(Bacteriome_No$feature, start = -6, end = -2)) -> Bacteriome_No
left_join(Bacteriome_No, US2_SP, by = "ID") -> US2_SP
column_to_rownames(US2_SP, "feature.x") -> US2_SP
Bacteriome <- US2_SP[, -which (colnames(US2_SP) %in% c("Biome", "ID", "feature.y"))] 
Bacteriome <- as.data.frame(t(Bacteriome))

#Virome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
SP_WO3<- read.csv("phage.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = TRUE)
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")
#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
#prop.table(data.matrix(US2_SP), 1) -> US2_SP (data is RPKM)
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
US2_SP$feature <- str_replace_all(US2_SP$feature, pattern="[:]", replacement=".") #Adjust the phage name before data integration
left_join(Virome_No, US2_SP, by = "feature") -> US2_SP #Please confirm all data is included
column_to_rownames(US2_SP, "No") -> US2_SP
Virome <- US2_SP[, -which (colnames(US2_SP) %in% c("feature"))] 
Virome <- as.data.frame(t(Virome))

#Mycobiome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Fungi")
SP_WO3<- read.csv("eukaryote.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = TRUE)
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")
#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
#prop.table(data.matrix(US2_SP), 1) -> US2_SP (data is RPKM)
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
Fungi_No$feature <- str_replace_all(Fungi_No$feature, pattern="[.]", replacement=" ") 
left_join(Fungi_No, US2_SP, by = "feature") -> US2_SP #Please confirm all data is included
column_to_rownames(US2_SP, "feature") -> Mycobiome
Mycobiome <- as.data.frame(t(Mycobiome))

#Network analysis
rownames_to_column(Bacteriome, "ID") -> Bacteriome
rownames_to_column(Virome, "ID") -> Virome
rownames_to_column(Mycobiome, "ID") -> Mycobiome
inner_join(inner_join(Bacteriome, Virome, by = "ID"), Mycobiome, by = "ID") -> SparData
column_to_rownames(.data = SparData, "ID") -> SparData

#Select important strains for validation
SparData %>% select("vOTU913(Bacteroides)", "vOTU914(Bacteroides)", "Bacteroides thetaiotaomicron [ID:01657]", 
                    "vOTU77(Escherichia)", "Escherichia coli [ID:00095]", 
                    "vOTU758(Anaerostipes)", "Anaerostipes hadrus [ID:00857]", "Anaerostipes hadrus [ID:00856]",
                    "vOTU689(Streptococcus)", "Streptococcus sp. [ID:00902]", 
                    "Saccharomyces cerevisiae", "Eubacterium ventriosum [ID:11151]", "Clostridiales species incertae sedis [ID:12288]",
                    "Clostridium sp. [ID:03680]", "Clostridiales species incertae sedis [ID:12563]", "Blautia obeum [ID:05139]", "Clostridiales species incertae sedis [ID:12284]",
                    "vOTU388(Clostridium, Faecalibacterium)", "Faecalibacterium prausnitzii [ID:06110]",
                    "[Clostridium] leptum [ID:03688]", "Clostridiales sp. [ID:03658]",  "Clostridium phoceensis [ID:08865]", "Faecalibacterium prausnitzii [ID:06108]", 
                    "vOTU1328(Unclassified Lachnospiraceae, Unclassified Clostridiales, Flavonifractor)", "Clostridiales species incertae sedis [ID:12473]", 
                    "vOTU756(Lachnoclostridium, Ruminococcus)", "Ruminococcaceae species incertae sedis [ID:12239]",
                    "vOTU497(Clostridium, Anaerotruncus, Ruminococcus, Blautia, Dorea, Coprococcus)", "Dorea longicatena [ID:03693]", "Blautia massiliensis [ID:03342]",
                    "Clostridiales species incertae sedis [ID:12476]","Dorea formicigenerans [ID:03668]", "Clostridiales species incertae sedis [ID:12554]", "Coprococcus comes [ID:03690]", 
                    "vOTU533(Roseburia, Coprococcus, Unclassified Firmicutes, Clostridium, Ruminococcus)", "[Ruminococcus] torques [ID:03659]", "Clostridiales sp. [ID:03658]", "uncultured Flavonifractor sp. [ID:07315]") -> SparData

net <- netConstruct(SparData, #relative abundance
                    filtTax = "none",
                    filtSamp = "none", 
                    filtSampPar = "none", 
                    measure = "spearman",
                    thresh = 0.2,
                    measurePar = list(nlambda=10, 
                                      rep.num=10),　
                    normMethod = "none",
                    sparsMethod = "t-test", #sparsMethod = "threshhold", thresh = 0.01
                    adjust = "fdr",
                    alpha = 0.1, 
                    lfdrThresh = 0.1,
                    sampleSize = 34,  
                    dissFunc = "signed",
                    verbose = 2,
                    seed = 123456)

net$assoMat1 %>% as.data.frame() -> col_fdr #if fdr<0.1, coef included, otherwise 0

net2 <- netConstruct(SparData, #relative abundance
                     filtTax = "none",
                     filtSamp = "none", 
                     filtSampPar = "none", 
                     measure = "spearman",
                     measurePar = list(nlambda=10, 
                                       rep.num=10),　
                     normMethod = "none",
                     sparsMethod = "none", #sparsMethod = "threshhold", thresh = 0.01
                     sampleSize = 34,  
                     dissFunc = "signed",
                     verbose = 2,
                     seed = 123456)

net2$assoMat1 %>% as.data.frame() -> col #fdr not related, pure Spearman correlation 

df_US <- data.frame(
  Correlation = c("vOTU913 & Bacteroides thetaiotaomicron", "vOTU914 & Bacteroides thetaiotaomicron",
                  "vOTU77 & Escherichia coli",
                  "vOTU758 & Anaerostipes hadrus[ID:00857]", "vOTU758 & Anaerostipes hadrus[ID:00856]", 
                  "vOTU689 & Streptococcus sp.",
                  "Saccharomyces cerevisiae & Eubacterium ventriosum", "Saccharomyces cerevisiae & Clostridiales species incertae sedis",
                  "Escherichia coli & Anaerostipes hadrus[ID:00857]", "Escherichia coli & Anaerostipes hadrus[ID:00856]","Escherichia coli & Blautia obeum", "Escherichia coli & Clostridium sp.","Escherichia coli & Clostridiales species incertae sedis[ID:12563]","Escherichia coli & Clostridiales species incertae sedis[ID:12284]"),
  Coefficient_Franzosa_US = c(col[3,1], col[3,2], 
                     col[4,5],
                     col[6,7],col[6,8],
                     col[9,10],
                     col[11,12],col[11,13],
                     col[5,7], col[5,8], col[5,16], col[5,14], col[5,15],col[5,17]
  )
)

df_fdr_US <- data.frame(
  Correlation = c("vOTU913 & Bacteroides thetaiotaomicron", "vOTU914 & Bacteroides thetaiotaomicron",
                  "vOTU77 & Escherichia coli",
                  "vOTU758 & Anaerostipes hadrus[ID:00857]", "vOTU758 & Anaerostipes hadrus[ID:00856]", 
                  "vOTU689 & Streptococcus sp.",
                  "Saccharomyces cerevisiae & Eubacterium ventriosum", "Saccharomyces cerevisiae & Clostridiales species incertae sedis",
                  "Escherichia coli & Anaerostipes hadrus[ID:00857]", "Escherichia coli & Anaerostipes hadrus[ID:00856]","Escherichia coli & Blautia obeum", "Escherichia coli & Clostridium sp.","Escherichia coli & Clostridiales species incertae sedis[ID:12563]","Escherichia coli & Clostridiales species incertae sedis[ID:12284]"),
  Coefficient_Franzosa_US = c(col_fdr[3,1], col_fdr[3,2], 
                     col_fdr[4,5],
                     col_fdr[6,7],col_fdr[6,8],
                     col_fdr[9,10],
                     col_fdr[11,12],col_fdr[11,13],
                     col_fdr[5,7], col_fdr[5,8], col_fdr[5,16], col_fdr[5,14], col_fdr[5,15],col_fdr[5,17]
  )
)

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")
p_USc=pheatmap(as.matrix(col), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(col_fdr >0 | col_fdr<0,"*", ""), nrow(col_fdr)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.5,0,1), c("navy", "white", "firebrick3")), name = "Spearman rho", heatmap_legend_param = list(color_bar = "continuous"), column_title = "Correlations of CD-related multibiome (US cohort)", column_title_gp = gpar(fontsize = 10, fontface = "bold"))
draw(p_USc, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig)) 

#External data2
#Repeat from here
SP_input2_bacteriome %>% filter(abs(SP_input2_bacteriome$"Coefficient (CD JP)") > 1 & SP_input2_bacteriome$"Q-value (CD JP)" < 0.1) -> SP_input2_bacteriome2 #For CD 
SP_input2_bacteriome2 %>% dplyr::select(feature) -> Bacteriome_No 
#Phage list 
#vOTU_Heat_qval_fig from virome analysis including vOTU with abs coef >1 and FDR<0.1 in UC or CD
vOTU_Heat_qval_fig %>% filter(abs(vOTU_Heat_qval_fig $"Coefficient (CD JP)") > 1 & vOTU_Heat_qval_fig$"Q-value (CD JP)" < 0.1) -> vOTU_Heat_qval_fig2 #For CD
rownames_to_column(vOTU_Heat_qval_fig2, "No") -> vOTU_Heat_qval_fig2
vOTU_Heat_qval_fig2 %>% dplyr::select("No", "feature") -> Virome_No 

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
Ref <- read.csv("Reference2.csv", header = TRUE, na.strings = c(NA, ''))
left_join(Virome_No, Ref, by = "No") ->Virome_No
Virome_No <- Virome_No[, -which (colnames(Virome_No) %in% c("No"))] 
rename(.data=Virome_No, "No" = "vOTU_host") -> Virome_No

#Fungi list 
#Fungi_map from mycobiome analysis including fungi with P<0.05 in UC or CD
Fungi_map %>% filter(Fungi_map$"P-value (CD JP)"< 0.05) -> Fungi_map2 #For CD
Fungi_map2 %>% dplyr::select("feature") -> Fungi_No

#External dataset
X<- "Franzosa_2018_IBD" #Select data, "He_2017_Crohn", "ES_MH, MetaHIT"
Y <- "Netherlands" #Select location #"Netherlands", "Spain", "China"

#Bacteriome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/mOTU3_original_data")
SP_WO3<- read.csv("Bork_group.mOUTs3.csv", header = TRUE, na.strings = c(NA, ''))
rename(.data=SP_WO3, "sample_id" = "X",) -> SP_WO3
rename(.data=SP_WO3, "Unassigned species" = "X.1",) -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")

#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
prop.table(data.matrix(US2_SP), 1) -> US2_SP #Bacteriome data is not normalized. This script is needed.
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
US2_SP %>% mutate(ID = str_sub(US2_SP$feature, start = -6, end = -2)) -> US2_SP
Bacteriome_No %>% mutate(ID = str_sub(Bacteriome_No$feature, start = -6, end = -2)) -> Bacteriome_No
left_join(Bacteriome_No, US2_SP, by = "ID") -> US2_SP
column_to_rownames(US2_SP, "feature.x") -> US2_SP
Bacteriome <- US2_SP[, -which (colnames(US2_SP) %in% c("Biome", "ID", "feature.y"))] 
Bacteriome <- as.data.frame(t(Bacteriome))

#Virome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
SP_WO3<- read.csv("phage.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = TRUE)
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")
#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
#prop.table(data.matrix(US2_SP), 1) -> US2_SP (data is RPKM)
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
US2_SP$feature <- str_replace_all(US2_SP$feature, pattern="[:]", replacement=".") #Adjust the phage name before data integration
left_join(Virome_No, US2_SP, by = "feature") -> US2_SP #Please confirm all data is included
column_to_rownames(US2_SP, "No") -> US2_SP
Virome <- US2_SP[, -which (colnames(US2_SP) %in% c("feature"))] 
Virome <- as.data.frame(t(Virome))

#Mycobiome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Fungi")
SP_WO3<- read.csv("eukaryote.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = TRUE)
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")
#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
#prop.table(data.matrix(US2_SP), 1) -> US2_SP (data is RPKM)
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
Fungi_No$feature <- str_replace_all(Fungi_No$feature, pattern="[.]", replacement=" ") 
left_join(Fungi_No, US2_SP, by = "feature") -> US2_SP #Please confirm all data is included
column_to_rownames(US2_SP, "feature") -> Mycobiome
Mycobiome <- as.data.frame(t(Mycobiome))

#Network analysis
rownames_to_column(Bacteriome, "ID") -> Bacteriome
rownames_to_column(Virome, "ID") -> Virome
rownames_to_column(Mycobiome, "ID") -> Mycobiome
inner_join(inner_join(Bacteriome, Virome, by = "ID"), Mycobiome, by = "ID") -> SparData
column_to_rownames(.data = SparData, "ID") -> SparData

#Select important strains for validation
SparData %>% select("vOTU913(Bacteroides)", "vOTU914(Bacteroides)", "Bacteroides thetaiotaomicron [ID:01657]", 
                    "vOTU77(Escherichia)", "Escherichia coli [ID:00095]", 
                    "vOTU758(Anaerostipes)", "Anaerostipes hadrus [ID:00857]", "Anaerostipes hadrus [ID:00856]",
                    "vOTU689(Streptococcus)", "Streptococcus sp. [ID:00902]", 
                    "Saccharomyces cerevisiae", "Eubacterium ventriosum [ID:11151]", "Clostridiales species incertae sedis [ID:12288]",
                    "Clostridium sp. [ID:03680]", "Clostridiales species incertae sedis [ID:12563]", "Blautia obeum [ID:05139]", "Clostridiales species incertae sedis [ID:12284]",
                    "vOTU388(Clostridium, Faecalibacterium)", "Faecalibacterium prausnitzii [ID:06110]",
                    "[Clostridium] leptum [ID:03688]", "Clostridiales sp. [ID:03658]",  "Clostridium phoceensis [ID:08865]", "Faecalibacterium prausnitzii [ID:06108]", 
                    "vOTU1328(Unclassified Lachnospiraceae, Unclassified Clostridiales, Flavonifractor)", "Clostridiales species incertae sedis [ID:12473]", 
                    "vOTU756(Lachnoclostridium, Ruminococcus)", "Ruminococcaceae species incertae sedis [ID:12239]",
                    "vOTU497(Clostridium, Anaerotruncus, Ruminococcus, Blautia, Dorea, Coprococcus)", "Dorea longicatena [ID:03693]", "Blautia massiliensis [ID:03342]",
                    "Clostridiales species incertae sedis [ID:12476]","Dorea formicigenerans [ID:03668]", "Clostridiales species incertae sedis [ID:12554]", "Coprococcus comes [ID:03690]", 
                    "vOTU533(Roseburia, Coprococcus, Unclassified Firmicutes, Clostridium, Ruminococcus)", "[Ruminococcus] torques [ID:03659]", "Clostridiales sp. [ID:03658]", "uncultured Flavonifractor sp. [ID:07315]") -> SparData

net <- netConstruct(SparData, #relative abundance
                    filtTax = "none",
                    filtSamp = "none", 
                    filtSampPar = "none", 
                    measure = "spearman",
                    thresh = 0.2,
                    measurePar = list(nlambda=10, 
                                      rep.num=10),　
                    normMethod = "none",
                    sparsMethod = "t-test", #sparsMethod = "threshhold", thresh = 0.01
                    adjust = "fdr",
                    alpha = 0.1, 
                    lfdrThresh = 0.1,
                    sampleSize = 22,  
                    dissFunc = "signed",
                    verbose = 2,
                    seed = 123456)

net$assoMat1 %>% as.data.frame() -> col_fdr #if fdr<0.1, coef included, otherwise 0

net2 <- netConstruct(SparData, #relative abundance
                     filtTax = "none",
                     filtSamp = "none", 
                     filtSampPar = "none", 
                     measure = "spearman",
                     measurePar = list(nlambda=10, 
                                       rep.num=10),　
                     normMethod = "none",
                     sparsMethod = "none", #sparsMethod = "threshhold", thresh = 0.01
                     sampleSize = 22,  
                     dissFunc = "signed",
                     verbose = 2,
                     seed = 123456)

net2$assoMat1 %>% as.data.frame() -> col #fdr not related, pure Spearman correlation 

df_Netherlands <- data.frame(
  Correlation = c("vOTU913 & Bacteroides thetaiotaomicron", "vOTU914 & Bacteroides thetaiotaomicron",
                  "vOTU77 & Escherichia coli",
                  "vOTU758 & Anaerostipes hadrus[ID:00857]", "vOTU758 & Anaerostipes hadrus[ID:00856]", 
                  "vOTU689 & Streptococcus sp.",
                  "Saccharomyces cerevisiae & Eubacterium ventriosum", "Saccharomyces cerevisiae & Clostridiales species incertae sedis",
                  "Escherichia coli & Anaerostipes hadrus[ID:00857]", "Escherichia coli & Anaerostipes hadrus[ID:00856]","Escherichia coli & Blautia obeum", "Escherichia coli & Clostridium sp.","Escherichia coli & Clostridiales species incertae sedis[ID:12563]","Escherichia coli & Clostridiales species incertae sedis[ID:12284]"),
  Coefficient_Franzosa_Netherlands = c(col[3,1], col[3,2], 
                              col[4,5],
                              col[6,7],col[6,8],
                              col[9,10],
                              col[11,12],col[11,13],
                              col[5,7], col[5,8], col[5,16], col[5,14], col[5,15],col[5,17]
  )
)

df_fdr_Netherlands <- data.frame(
  Correlation = c("vOTU913 & Bacteroides thetaiotaomicron", "vOTU914 & Bacteroides thetaiotaomicron",
                  "vOTU77 & Escherichia coli",
                  "vOTU758 & Anaerostipes hadrus[ID:00857]", "vOTU758 & Anaerostipes hadrus[ID:00856]", 
                  "vOTU689 & Streptococcus sp.",
                  "Saccharomyces cerevisiae & Eubacterium ventriosum", "Saccharomyces cerevisiae & Clostridiales species incertae sedis",
                  "Escherichia coli & Anaerostipes hadrus[ID:00857]", "Escherichia coli & Anaerostipes hadrus[ID:00856]","Escherichia coli & Blautia obeum", "Escherichia coli & Clostridium sp.","Escherichia coli & Clostridiales species incertae sedis[ID:12563]","Escherichia coli & Clostridiales species incertae sedis[ID:12284]"),
  Coefficient_Franzosa_Netherlands = c(col_fdr[3,1], col_fdr[3,2], 
                              col_fdr[4,5],
                              col_fdr[6,7],col_fdr[6,8],
                              col_fdr[9,10],
                              col_fdr[11,12],col_fdr[11,13],
                              col_fdr[5,7], col_fdr[5,8], col_fdr[5,16], col_fdr[5,14], col_fdr[5,15],col_fdr[5,17]
  )
)

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")
p_NLc=pheatmap(as.matrix(col), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(col_fdr >0 | col_fdr<0,"*", ""), nrow(col_fdr)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.5,0,1), c("navy", "white", "firebrick3")), name = "Spearman rho", heatmap_legend_param = list(color_bar = "continuous"), column_title = "Correlations of CD-related multibiome (NL cohort)", column_title_gp = gpar(fontsize = 10, fontface = "bold"))
draw(p_NLc, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig)) 

#External dataset3
#Repeat from here
SP_input2_bacteriome %>% filter(abs(SP_input2_bacteriome$"Coefficient (CD JP)") > 1 & SP_input2_bacteriome$"Q-value (CD JP)" < 0.1) -> SP_input2_bacteriome2 #For CD 
SP_input2_bacteriome2 %>% dplyr::select(feature) -> Bacteriome_No 
#Phage list 
#vOTU_Heat_qval_fig from virome analysis including vOTU with abs coef >1 and FDR<0.1 in UC or CD
vOTU_Heat_qval_fig %>% filter(abs(vOTU_Heat_qval_fig $"Coefficient (CD JP)") > 1 & vOTU_Heat_qval_fig$"Q-value (CD JP)" < 0.1) -> vOTU_Heat_qval_fig2 #For CD
rownames_to_column(vOTU_Heat_qval_fig2, "No") -> vOTU_Heat_qval_fig2
vOTU_Heat_qval_fig2 %>% dplyr::select("No", "feature") -> Virome_No 

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
Ref <- read.csv("Reference2.csv", header = TRUE, na.strings = c(NA, ''))
left_join(Virome_No, Ref, by = "No") ->Virome_No
Virome_No <- Virome_No[, -which (colnames(Virome_No) %in% c("No"))] 
rename(.data=Virome_No, "No" = "vOTU_host") -> Virome_No

#Fungi list 
#Fungi_map from mycobiome analysis including fungi with P<0.05 in UC or CD
Fungi_map %>% filter(Fungi_map$"P-value (CD JP)"< 0.05) -> Fungi_map2 #For CD
Fungi_map2 %>% dplyr::select("feature") -> Fungi_No

X<- "ES_MH, MetaHIT" #Select data, "He_2017_Crohn", "ES_MH, MetaHIT"
Y <- "Spain" #Select location #"Netherlands", "Spain", "China"

#Bacteriome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/mOTU3_original_data")
SP_WO3<- read.csv("Bork_group.mOUTs3.csv", header = TRUE, na.strings = c(NA, ''))
rename(.data=SP_WO3, "sample_id" = "X",) -> SP_WO3
rename(.data=SP_WO3, "Unassigned species" = "X.1",) -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")

#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
prop.table(data.matrix(US2_SP), 1) -> US2_SP #Bacteriome data is not normalized. This script is needed.
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
US2_SP %>% mutate(ID = str_sub(US2_SP$feature, start = -6, end = -2)) -> US2_SP
Bacteriome_No %>% mutate(ID = str_sub(Bacteriome_No$feature, start = -6, end = -2)) -> Bacteriome_No
left_join(Bacteriome_No, US2_SP, by = "ID") -> US2_SP
column_to_rownames(US2_SP, "feature.x") -> US2_SP
Bacteriome <- US2_SP[, -which (colnames(US2_SP) %in% c("Biome", "ID", "feature.y"))] 
Bacteriome <- as.data.frame(t(Bacteriome))

#Virome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
SP_WO3<- read.csv("phage.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = TRUE)
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")
#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
#prop.table(data.matrix(US2_SP), 1) -> US2_SP (data is RPKM)
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
US2_SP$feature <- str_replace_all(US2_SP$feature, pattern="[:]", replacement=".") #Adjust the phage name before data integration
left_join(Virome_No, US2_SP, by = "feature") -> US2_SP #Please confirm all data is included
column_to_rownames(US2_SP, "No") -> US2_SP
Virome <- US2_SP[, -which (colnames(US2_SP) %in% c("feature"))] 
Virome <- as.data.frame(t(Virome))

#Mycobiome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Fungi")
SP_WO3<- read.csv("eukaryote.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = TRUE)
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")
#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
#prop.table(data.matrix(US2_SP), 1) -> US2_SP (data is RPKM)
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
Fungi_No$feature <- str_replace_all(Fungi_No$feature, pattern="[.]", replacement=" ") 
left_join(Fungi_No, US2_SP, by = "feature") -> US2_SP #Please confirm all data is included
column_to_rownames(US2_SP, "feature") -> Mycobiome
Mycobiome <- as.data.frame(t(Mycobiome))

#Network analysis
rownames_to_column(Bacteriome, "ID") -> Bacteriome
rownames_to_column(Virome, "ID") -> Virome
rownames_to_column(Mycobiome, "ID") -> Mycobiome
inner_join(inner_join(Bacteriome, Virome, by = "ID"), Mycobiome, by = "ID") -> SparData
column_to_rownames(.data = SparData, "ID") -> SparData

SparData <- na.omit(SparData) #Spain data has deficient data of fungi and phage

#Select important strains for validation
SparData %>% select("vOTU913(Bacteroides)", "vOTU914(Bacteroides)", "Bacteroides thetaiotaomicron [ID:01657]", 
                    "vOTU77(Escherichia)", "Escherichia coli [ID:00095]", 
                    "vOTU758(Anaerostipes)", "Anaerostipes hadrus [ID:00857]", "Anaerostipes hadrus [ID:00856]",
                    "vOTU689(Streptococcus)", "Streptococcus sp. [ID:00902]", 
                    "Saccharomyces cerevisiae", "Eubacterium ventriosum [ID:11151]", "Clostridiales species incertae sedis [ID:12288]",
                    "Clostridium sp. [ID:03680]", "Clostridiales species incertae sedis [ID:12563]", "Blautia obeum [ID:05139]", "Clostridiales species incertae sedis [ID:12284]",
                    "vOTU388(Clostridium, Faecalibacterium)", "Faecalibacterium prausnitzii [ID:06110]",
                    "[Clostridium] leptum [ID:03688]", "Clostridiales sp. [ID:03658]",  "Clostridium phoceensis [ID:08865]", "Faecalibacterium prausnitzii [ID:06108]", 
                    "vOTU1328(Unclassified Lachnospiraceae, Unclassified Clostridiales, Flavonifractor)", "Clostridiales species incertae sedis [ID:12473]", 
                    "vOTU756(Lachnoclostridium, Ruminococcus)", "Ruminococcaceae species incertae sedis [ID:12239]",
                    "vOTU497(Clostridium, Anaerotruncus, Ruminococcus, Blautia, Dorea, Coprococcus)", "Dorea longicatena [ID:03693]", "Blautia massiliensis [ID:03342]",
                    "Clostridiales species incertae sedis [ID:12476]","Dorea formicigenerans [ID:03668]", "Clostridiales species incertae sedis [ID:12554]", "Coprococcus comes [ID:03690]", 
                    "vOTU533(Roseburia, Coprococcus, Unclassified Firmicutes, Clostridium, Ruminococcus)", "[Ruminococcus] torques [ID:03659]", "Clostridiales sp. [ID:03658]", "uncultured Flavonifractor sp. [ID:07315]") -> SparData

net <- netConstruct(SparData, #relative abundance
                    filtTax = "none",
                    filtSamp = "none", 
                    filtSampPar = "none", 
                    measure = "spearman",
                    thresh = 0.2,
                    measurePar = list(nlambda=10, 
                                      rep.num=10),　
                    normMethod = "none",
                    sparsMethod = "t-test", #sparsMethod = "threshhold", thresh = 0.01
                    adjust = "fdr",
                    alpha = 0.1, 
                    lfdrThresh = 0.1,
                    sampleSize = 62,  
                    dissFunc = "signed",
                    verbose = 2,
                    seed = 123456)

net$assoMat1 %>% as.data.frame() -> col_fdr #if fdr<0.1, coef included, otherwise 0

net2 <- netConstruct(SparData, #relative abundance
                     filtTax = "none",
                     filtSamp = "none", 
                     filtSampPar = "none", 
                     measure = "spearman",
                     measurePar = list(nlambda=10, 
                                       rep.num=10),　
                     normMethod = "none",
                     sparsMethod = "none", #sparsMethod = "threshhold", thresh = 0.01
                     sampleSize = 62,  
                     dissFunc = "signed",
                     verbose = 2,
                     seed = 123456)

net2$assoMat1 %>% as.data.frame() -> col #fdr not related, pure Spearman correlation 

df_Spain <- data.frame(
  Correlation = c("vOTU913 & Bacteroides thetaiotaomicron", "vOTU914 & Bacteroides thetaiotaomicron",
                  "vOTU77 & Escherichia coli",
                  "vOTU758 & Anaerostipes hadrus[ID:00857]", "vOTU758 & Anaerostipes hadrus[ID:00856]", 
                  "vOTU689 & Streptococcus sp.",
                  "Saccharomyces cerevisiae & Eubacterium ventriosum", "Saccharomyces cerevisiae & Clostridiales species incertae sedis",
                  "Escherichia coli & Anaerostipes hadrus[ID:00857]", "Escherichia coli & Anaerostipes hadrus[ID:00856]","Escherichia coli & Blautia obeum", "Escherichia coli & Clostridium sp.","Escherichia coli & Clostridiales species incertae sedis[ID:12563]","Escherichia coli & Clostridiales species incertae sedis[ID:12284]"),
  Coefficient_Nielsen_Spain = c(col[3,1], col[3,2], 
                                       col[4,5],
                                       col[6,7],col[6,8],
                                       col[9,10],
                                       col[11,12],col[11,13],
                                       col[5,7], col[5,8], col[5,16], col[5,14], col[5,15],col[5,17]
  )
)

df_fdr_Spain <- data.frame(
  Correlation = c("vOTU913 & Bacteroides thetaiotaomicron", "vOTU914 & Bacteroides thetaiotaomicron",
                  "vOTU77 & Escherichia coli",
                  "vOTU758 & Anaerostipes hadrus[ID:00857]", "vOTU758 & Anaerostipes hadrus[ID:00856]", 
                  "vOTU689 & Streptococcus sp.",
                  "Saccharomyces cerevisiae & Eubacterium ventriosum", "Saccharomyces cerevisiae & Clostridiales species incertae sedis",
                  "Escherichia coli & Anaerostipes hadrus[ID:00857]", "Escherichia coli & Anaerostipes hadrus[ID:00856]","Escherichia coli & Blautia obeum", "Escherichia coli & Clostridium sp.","Escherichia coli & Clostridiales species incertae sedis[ID:12563]","Escherichia coli & Clostridiales species incertae sedis[ID:12284]"),
  Coefficient_Nielsen_Spain = c(col_fdr[3,1], col_fdr[3,2], 
                                       col_fdr[4,5],
                                       col_fdr[6,7],col_fdr[6,8],
                                       col_fdr[9,10],
                                       col_fdr[11,12],col_fdr[11,13],
                                       col_fdr[5,7], col_fdr[5,8], col_fdr[5,16], col_fdr[5,14], col_fdr[5,15],col_fdr[5,17]
  )
)

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")
p_ESc=pheatmap(as.matrix(col), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(col_fdr >0 | col_fdr<0,"*", ""), nrow(col_fdr)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.5,0,1), c("navy", "white", "firebrick3")), name = "Spearman rho", heatmap_legend_param = list(color_bar = "continuous"), column_title = "Correlations of CD-related multibiome (ES cohort)", column_title_gp = gpar(fontsize = 10, fontface = "bold"))
draw(p_ESc, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig)) 


#External dataset4
#Repeat from here
SP_input2_bacteriome %>% filter(abs(SP_input2_bacteriome$"Coefficient (CD JP)") > 1 & SP_input2_bacteriome$"Q-value (CD JP)" < 0.1) -> SP_input2_bacteriome2 #For CD 
SP_input2_bacteriome2 %>% dplyr::select(feature) -> Bacteriome_No 
#Phage list 
#vOTU_Heat_qval_fig from virome analysis including vOTU with abs coef >1 and FDR<0.1 in UC or CD
vOTU_Heat_qval_fig %>% filter(abs(vOTU_Heat_qval_fig $"Coefficient (CD JP)") > 1 & vOTU_Heat_qval_fig$"Q-value (CD JP)" < 0.1) -> vOTU_Heat_qval_fig2 #For CD
rownames_to_column(vOTU_Heat_qval_fig2, "No") -> vOTU_Heat_qval_fig2
vOTU_Heat_qval_fig2 %>% dplyr::select("No", "feature") -> Virome_No 

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
Ref <- read.csv("Reference2.csv", header = TRUE, na.strings = c(NA, ''))
left_join(Virome_No, Ref, by = "No") ->Virome_No
Virome_No <- Virome_No[, -which (colnames(Virome_No) %in% c("No"))] 
rename(.data=Virome_No, "No" = "vOTU_host") -> Virome_No

#Fungi list 
#Fungi_map from mycobiome analysis including fungi with P<0.05 in UC or CD
Fungi_map %>% filter(Fungi_map$"P-value (CD JP)"< 0.05) -> Fungi_map2 #For CD
Fungi_map2 %>% dplyr::select("feature") -> Fungi_No

X<- "He_2017_Crohn" #Select data, "He_2017_Crohn", "ES_MH, MetaHIT"
Y <- "China" #Select location #"Netherlands", "Spain", "China"

#Bacteriome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/mOTU3_original_data")
SP_WO3<- read.csv("Bork_group.mOUTs3.csv", header = TRUE, na.strings = c(NA, ''))
rename(.data=SP_WO3, "sample_id" = "X",) -> SP_WO3
rename(.data=SP_WO3, "Unassigned species" = "X.1",) -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")

#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
prop.table(data.matrix(US2_SP), 1) -> US2_SP #Bacteriome data is not normalized. This script is needed.
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
US2_SP %>% mutate(ID = str_sub(US2_SP$feature, start = -6, end = -2)) -> US2_SP
Bacteriome_No %>% mutate(ID = str_sub(Bacteriome_No$feature, start = -6, end = -2)) -> Bacteriome_No
left_join(Bacteriome_No, US2_SP, by = "ID") -> US2_SP
column_to_rownames(US2_SP, "feature.x") -> US2_SP
Bacteriome <- US2_SP[, -which (colnames(US2_SP) %in% c("Biome", "ID", "feature.y"))] 
Bacteriome <- as.data.frame(t(Bacteriome))

#Virome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
SP_WO3<- read.csv("phage.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = TRUE)
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")
#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
#prop.table(data.matrix(US2_SP), 1) -> US2_SP (data is RPKM)
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
US2_SP$feature <- str_replace_all(US2_SP$feature, pattern="[:]", replacement=".") #Adjust the phage name before data integration
left_join(Virome_No, US2_SP, by = "feature") -> US2_SP #Please confirm all data is included
column_to_rownames(US2_SP, "No") -> US2_SP
Virome <- US2_SP[, -which (colnames(US2_SP) %in% c("feature"))] 
Virome <- as.data.frame(t(Virome))

#Mycobiome (External)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Fungi")
SP_WO3<- read.csv("eukaryote.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1, check.names = TRUE)
SP_WO3 %>% t() %>% as.data.frame() -> SP_WO3
rownames_to_column(SP_WO3, "sample_id") -> SP_WO3
Meta_WO3<- read.csv("Bork_group.metadata.csv", header = TRUE, na.strings = c(NA, ''), row.names = 1)
WO3 <- full_join(SP_WO3, Meta_WO3, by ="sample_id")
#Select Study
WO3 %>% filter(WO3$study == X & WO3$geographic_location == Y) -> US2 #timepoint=0がall 
US2 %>% distinct(subject_id, .keep_all=TRUE) -> US2 #remove duplicates
US2 %>% filter(US2$disease == "Control") -> US2 #Select Cont from US data
US2_SP <- US2[, -which (colnames(US2) %in% c("subject_id", "environment_material", "timepoint", 
                                             "Sex", "Age", "geographic_location", "disease", "study", "publications", "environment_feature", "collection_date", "intervention", "weight_kg", 
                                             "height_cm", "bmi", "diet", "smoker", "antibiotic", "bristol_stool_scale"))] 
column_to_rownames(US2_SP, "sample_id") -> US2_SP
#prop.table(data.matrix(US2_SP), 1) -> US2_SP (data is RPKM)
#Data_integration 
US2_SP %>% t() %>% as.data.frame -> US2_SP
rownames_to_column(US2_SP, "feature") -> US2_SP
Fungi_No$feature <- str_replace_all(Fungi_No$feature, pattern="[.]", replacement=" ") 
left_join(Fungi_No, US2_SP, by = "feature") -> US2_SP #Please confirm all data is included
column_to_rownames(US2_SP, "feature") -> Mycobiome
Mycobiome <- as.data.frame(t(Mycobiome))

#Network analysis
rownames_to_column(Bacteriome, "ID") -> Bacteriome
rownames_to_column(Virome, "ID") -> Virome
rownames_to_column(Mycobiome, "ID") -> Mycobiome
inner_join(inner_join(Bacteriome, Virome, by = "ID"), Mycobiome, by = "ID") -> SparData
column_to_rownames(.data = SparData, "ID") -> SparData

#Select important strains for validation
SparData %>% select("vOTU913(Bacteroides)", "vOTU914(Bacteroides)", "Bacteroides thetaiotaomicron [ID:01657]", 
                    "vOTU77(Escherichia)", "Escherichia coli [ID:00095]", 
                    "vOTU758(Anaerostipes)", "Anaerostipes hadrus [ID:00857]", "Anaerostipes hadrus [ID:00856]",
                    "vOTU689(Streptococcus)", "Streptococcus sp. [ID:00902]", 
                    "Saccharomyces cerevisiae", "Eubacterium ventriosum [ID:11151]", "Clostridiales species incertae sedis [ID:12288]",
                    "Clostridium sp. [ID:03680]", "Clostridiales species incertae sedis [ID:12563]", "Blautia obeum [ID:05139]", "Clostridiales species incertae sedis [ID:12284]",
                    "vOTU388(Clostridium, Faecalibacterium)", "Faecalibacterium prausnitzii [ID:06110]",
                    "[Clostridium] leptum [ID:03688]", "Clostridiales sp. [ID:03658]",  "Clostridium phoceensis [ID:08865]", "Faecalibacterium prausnitzii [ID:06108]", 
                    "vOTU1328(Unclassified Lachnospiraceae, Unclassified Clostridiales, Flavonifractor)", "Clostridiales species incertae sedis [ID:12473]", 
                    "vOTU756(Lachnoclostridium, Ruminococcus)", "Ruminococcaceae species incertae sedis [ID:12239]",
                    "vOTU497(Clostridium, Anaerotruncus, Ruminococcus, Blautia, Dorea, Coprococcus)", "Dorea longicatena [ID:03693]", "Blautia massiliensis [ID:03342]",
                    "Clostridiales species incertae sedis [ID:12476]","Dorea formicigenerans [ID:03668]", "Clostridiales species incertae sedis [ID:12554]", "Coprococcus comes [ID:03690]", 
                    "vOTU533(Roseburia, Coprococcus, Unclassified Firmicutes, Clostridium, Ruminococcus)", "[Ruminococcus] torques [ID:03659]", "Clostridiales sp. [ID:03658]", "uncultured Flavonifractor sp. [ID:07315]") -> SparData

net <- netConstruct(SparData, #relative abundance
                    filtTax = "none",
                    filtSamp = "none", 
                    filtSampPar = "none", 
                    measure = "spearman",
                    thresh = 0.2,
                    measurePar = list(nlambda=10, 
                                      rep.num=10),　
                    normMethod = "none",
                    sparsMethod = "t-test", #sparsMethod = "threshhold", thresh = 0.01
                    adjust = "fdr",
                    alpha = 0.1, 
                    lfdrThresh = 0.1,
                    sampleSize = 60, 
                    dissFunc = "signed",
                    verbose = 2,
                    seed = 123456)

net$assoMat1 %>% as.data.frame() -> col_fdr #if fdr<0.1, coef included, otherwise 0

net2 <- netConstruct(SparData, #relative abundance
                     filtTax = "none",
                     filtSamp = "none", 
                     filtSampPar = "none", 
                     measure = "spearman",
                     measurePar = list(nlambda=10, 
                                       rep.num=10),　
                     normMethod = "none",
                     sparsMethod = "none", #sparsMethod = "threshhold", thresh = 0.01
                     sampleSize = 60, 
                     dissFunc = "signed",
                     verbose = 2,
                     seed = 123456)

net2$assoMat1 %>% as.data.frame() -> col #fdr not related, pure Spearman correlation 

df_China <- data.frame(
  Correlation = c("vOTU913 & Bacteroides thetaiotaomicron", "vOTU914 & Bacteroides thetaiotaomicron",
                  "vOTU77 & Escherichia coli",
                  "vOTU758 & Anaerostipes hadrus[ID:00857]", "vOTU758 & Anaerostipes hadrus[ID:00856]", 
                  "vOTU689 & Streptococcus sp.",
                  "Saccharomyces cerevisiae & Eubacterium ventriosum", "Saccharomyces cerevisiae & Clostridiales species incertae sedis",
                  "Escherichia coli & Anaerostipes hadrus[ID:00857]", "Escherichia coli & Anaerostipes hadrus[ID:00856]","Escherichia coli & Blautia obeum", "Escherichia coli & Clostridium sp.","Escherichia coli & Clostridiales species incertae sedis[ID:12563]","Escherichia coli & Clostridiales species incertae sedis[ID:12284]"),
  Coefficient_He_China = c(col[3,1], col[3,2], 
                                col[4,5],
                                col[6,7],col[6,8],
                                col[9,10],
                                col[11,12],col[11,13],
                                col[5,7], col[5,8], col[5,16], col[5,14], col[5,15],col[5,17]
  )
)

df_fdr_China <- data.frame(
  Correlation = c("vOTU913 & Bacteroides thetaiotaomicron", "vOTU914 & Bacteroides thetaiotaomicron",
                  "vOTU77 & Escherichia coli",
                  "vOTU758 & Anaerostipes hadrus[ID:00857]", "vOTU758 & Anaerostipes hadrus[ID:00856]", 
                  "vOTU689 & Streptococcus sp.",
                  "Saccharomyces cerevisiae & Eubacterium ventriosum", "Saccharomyces cerevisiae & Clostridiales species incertae sedis",
                  "Escherichia coli & Anaerostipes hadrus[ID:00857]", "Escherichia coli & Anaerostipes hadrus[ID:00856]","Escherichia coli & Blautia obeum", "Escherichia coli & Clostridium sp.","Escherichia coli & Clostridiales species incertae sedis[ID:12563]","Escherichia coli & Clostridiales species incertae sedis[ID:12284]"),
  Coefficient_He_China = c(col_fdr[3,1], col_fdr[3,2], 
                                col_fdr[4,5],
                                col_fdr[6,7],col_fdr[6,8],
                                col_fdr[9,10],
                                col_fdr[11,12],col_fdr[11,13],
                                col_fdr[5,7], col_fdr[5,8], col_fdr[5,16], col_fdr[5,14], col_fdr[5,15],col_fdr[5,17]
  )
)

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")
p_CNc=pheatmap(as.matrix(col), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(col_fdr >0 | col_fdr<0,"*", ""), nrow(col_fdr)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.5,0,1), c("navy", "white", "firebrick3")), name = "Spearman rho", heatmap_legend_param = list(color_bar = "continuous"), column_title = "Correlations of CD-related multibiome (CN cohort)", column_title_gp = gpar(fontsize = 10, fontface = "bold"))
draw(p_CNc, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig)) 


full_join(full_join(full_join(full_join(df, df_US, by = "Correlation"), df_Spain, by = "Correlation"), df_Netherlands, by = "Correlation"), df_China, by ="Correlation") -> all_df
full_join(full_join(full_join(full_join(df_fdr, df_fdr_US, by = "Correlation"), df_fdr_Spain, by = "Correlation"), df_fdr_Netherlands, by = "Correlation"), df_fdr_China, by = "Correlation") -> all_df_fdr

column_to_rownames(all_df, "Correlation") -> all_df
column_to_rownames(all_df_fdr, "Correlation") -> all_df_fdr

rename(.data = all_df, "Japanese 4D cohort" = "Coefficient_JP") -> all_df
rename(.data = all_df, "US cohort" = "Coefficient_Franzosa_US") -> all_df
rename(.data = all_df, "NL cohort" = "Coefficient_Franzosa_Netherlands") -> all_df
rename(.data = all_df, "ES cohort" = "Coefficient_Nielsen_Spain") -> all_df
rename(.data = all_df, "CN cohort" = "Coefficient_He_China") -> all_df

rename(.data = all_df_fdr, "Japanese 4D cohort" = "Coefficient_JP") -> all_df_fdr
rename(.data = all_df_fdr, "US cohort" = "Coefficient_Franzosa_US") -> all_df_fdr
rename(.data = all_df_fdr, "NL cohort" = "Coefficient_Franzosa_Netherlands") -> all_df_fdr
rename(.data = all_df_fdr, "ES cohort" = "Coefficient_Nielsen_Spain") -> all_df_fdr
rename(.data = all_df_fdr, "CN cohort" = "Coefficient_He_China") -> all_df_fdr

lgd_sig = Legend(pch = "*", type = "points", labels = "FDR < 0.1")
p =pheatmap(as.matrix(all_df), fontsize = 5, cellwidth = 8, cellheight = 8, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_col = 8, fontsize_row = 8, display_numbers = matrix(ifelse(all_df_fdr != 0,"*", ""), nrow(all_df_fdr)), fontsize_number = 7, border_color = "black", col = circlize::colorRamp2(c(-0.5,0,1), c("navy", "white", "firebrick3")), name = "Spearman rho", heatmap_legend_param = list(color_bar = "continuous"), column_title = "Correlations of CD-related species", column_title_gp = gpar(fontsize = 10, fontface = "bold"))
draw(p, heatmap_legend_side = "left", annotation_legend_side = "left", annotation_legend_list = list(lgd_sig)) 

#Script double check done (7/28/2024)
