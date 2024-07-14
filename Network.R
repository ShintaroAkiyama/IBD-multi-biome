library(dplyr)
library(tidyverse)
library(NetCoMi)

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

#External dataset
X<- "Franzosa_2018_IBD" #Select data
Y <- "USA" #Select location

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

SP_input2_bacteriome %>% filter(abs(SP_input2_bacteriome$"Coefficient (UC JP)") > 1 & SP_input2_bacteriome$"Q-value (UC JP)" < 0.1) -> SP_input2_bacteriome2 #For UC
SP_input2_bacteriome2 %>% dplyr::select(feature) -> Bacteriome_No 
Bacteriome_No %>% mutate("Biome" = "Bacteriome") -> Bacteriome_No

Virome_No %>% mutate("Biome" = "Virome") -> Virome_No 
Virome_No <- Virome_No[, -which (colnames(Virome_No) %in% c("feature"))] 
rename(.data = Virome_No, "feature" = "No") -> Virome_No
Fungi_No %>% mutate("Biome" = "Mycobiome") -> Fungi_No 
rbind(Bacteriome_No, Virome_No, Fungi_No) -> Biome
Bio <- as.factor(Biome[, "Biome"])
names(Bio) <- Biome[, "feature"]

biomecol <- c("#ff4b00", "#77d9a8", "dodgerblue3")

net <- netConstruct(SparData, #relative abundance
                    filtTax = "none",
                    filtSamp = "none", 
                    filtSampPar = "none", 
                    measure = "spearman",
                    measurePar = list(nlambda=10, 
                                      rep.num=10),　
                    normMethod = "none",
                    sparsMethod = "t-test", #sparsMethod = "threshhold", thresh = 0.01
                    adjust = "fdr",
                    alpha = 0.1, 
                    sampleSize = 540,  #注意して。UC111, CDでは31, Contは540 
                    dissFunc = "signed",
                    verbose = 2,
                    seed = 123456)

props <- netAnalyze(net, 
                    centrLCC = TRUE,
                    clustMethod = "cluster_fast_greedy",
                    hubPar = "eigenvector",
                    weightDeg = FALSE, normDeg = FALSE)

summary(props, numbNodes = 1L)

#https://rdrr.io/github/stefpeschel/NetCoMi/man/plot.microNetProps.html

p <- plot(props, 
          #nodeFilter = "highestConnect",
          #nodeFilterPar = 200,
          edgeInvisFilter = "threshold",
          edgeInvisPar = 0.2, 
          nodeSize = "eigenvector",
          repulsion =  0.7, 
          rmSingles = "all", 
          labelScale = FALSE,
          cexLabels = 0.8,  
          nodeSizeSpread = 3,
          cexNodes = 2,
          nodeColor = "feature", 
          featVecCol = Bio, 
          colorVec =  biomecol,
          posCol = "darkturquoise", 
          negCol = "orange",
          edgeTranspLow = 0,
          edgeTranspHigh = 50,
          hubBorderCol = "darkgray",
          title1 = "Network of UC multi-biome in controls", #CD (n = 31) Control (n = 540) 
          showTitle = TRUE,
          cexTitle = 1)

biomecol_transp <- colToTransp(biomecol, 60)

legend(-1, 1.1, cex = 1, pt.cex = 2.5, title = "Biome:", 
       legend=levels(Bio), col = biomecol_transp, bty = "n", pch = 16) 

legend(0.7, 1.1, cex = 1, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("darkturquoise","orange"), 
       bty = "n", horiz = TRUE)

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
SP_input %>% filter(SP_input$IBD ==2) -> SP_input  #Select Cont
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
left_join(NO, vOTU_input, by = "NO") -> vOTU_input #Select CD vs Cont
vOTU_input <- vOTU_input[, -which (colnames(vOTU_input) %in% c("IBD_HC"))] 
column_to_rownames(vOTU_input, "NO") -> vOTU_input
vOTU_input <- as.data.frame(t(vOTU_input)) 
rownames_to_column(vOTU_input, "feature") -> vOTU_input 
left_join(Virome_No, vOTU_input, by = "feature") -> vOTU_input2
vOTU_input2 <- vOTU_input2[, -which (colnames(vOTU_input2) %in% c("feature"))] #Remove non-essential variables
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

#External dataset
X<- "Franzosa_2018_IBD" #Select data
Y <- "USA" #Select location

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

SP_input2_bacteriome %>% filter(abs(SP_input2_bacteriome$"Coefficient (CD JP)") > 1 & SP_input2_bacteriome$"Q-value (CD JP)" < 0.1) -> SP_input2_bacteriome2 #For CD 
SP_input2_bacteriome2 %>% dplyr::select(feature) -> Bacteriome_No 
Bacteriome_No %>% mutate("Biome" = "Bacteriome") -> Bacteriome_No 

Virome_No %>% mutate("Biome" = "Virome") -> Virome_No 
Virome_No <- Virome_No[, -which (colnames(Virome_No) %in% c("feature"))] 
rename(.data = Virome_No, "feature" = "No") -> Virome_No
Fungi_No %>% mutate("Biome" = "Mycobiome") -> Fungi_No 
rbind(Bacteriome_No, Virome_No, Fungi_No) -> Biome
Bio <- as.factor(Biome[, "Biome"])
names(Bio) <- Biome[, "feature"]

biomecol <- c("#ff4b00", "#77d9a8", "dodgerblue3")

net <- netConstruct(SparData, #relative abundance
                    filtTax = "none",
                    filtSamp = "none", 
                    filtSampPar = "none", 
                    measure = "spearman",
                    measurePar = list(nlambda=10, 
                                      rep.num=10),　
                    normMethod = "none",
                    sparsMethod = "t-test", #sparsMethod = "threshhold", thresh = 0.01
                    adjust = "fdr",
                    alpha = 0.1, 
                    sampleSize = 540,  
                    dissFunc = "signed",
                    verbose = 2,
                    seed = 123456)

props <- netAnalyze(net, 
                    centrLCC = TRUE,
                    clustMethod = "cluster_fast_greedy",
                    hubPar = "eigenvector",
                    weightDeg = FALSE, normDeg = FALSE)

summary(props, numbNodes = 5L)

#https://rdrr.io/github/stefpeschel/NetCoMi/man/plot.microNetProps.html

p <- plot(props, 
          edgeInvisFilter = "threshold",
          edgeInvisPar = 0.2, 
          nodeSize = "eigenvector",
          repulsion =  0.8, #CDだと0.8
          rmSingles = "all",
          labelScale = FALSE,
          cexLabels = 0.6,  #CDだと0.6
          nodeSizeSpread = 3,
          cexNodes = 2,
          nodeColor = "feature", 
          featVecCol = Bio, 
          colorVec =  biomecol,
          posCol = "darkturquoise", 
          negCol = "orange",
          edgeTranspLow = 0,
          edgeTranspHigh = 50,
          hubBorderCol = "darkgray",
          title1 = "Network of CD multi-biome in controls",
          showTitle = TRUE,
          cexTitle = 1)

biomecol_transp <- colToTransp(biomecol, 60)

legend(-1, 1.1, cex = 1, pt.cex = 2.5, title = "Biome:", 
       legend=levels(Bio), col = biomecol_transp, bty = "n", pch = 16) 

legend(0.7, 1.1, cex = 1, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("darkturquoise","orange"), 
       bty = "n", horiz = TRUE)