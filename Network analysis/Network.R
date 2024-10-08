library(dplyr)
library(tidyverse)
library(NetCoMi)
library(ComplexHeatmap)
library(progress)

#UC analysis (Figure 5b)
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

#Can skip here
#Select important strains for the validation
SparData %>% select("vOTU160(Ruminococcus, Eubacterium)", "vOTU161(Ruminococcus, Eubacterium)", "vOTU162(Ruminococcus, Eubacterium)", "vOTU164(Ruminococcus, Eubacterium)", "Ruminococcus bromii [ID:00853]",
                    "vOTU870(Bacteroides, Alistipes)", "vOTU913(Bacteroides)", "Bacteroides thetaiotaomicron [ID:01657]",
                    "vOTU289(Blautia, Fusicatenibacter)", "vOTU291(Blautia, Fusicatenibacter)", "vOTU293(Blautia, Fusicatenibacter)", "vOTU294(Blautia, Fusicatenibacter)", "Fusicatenibacter saccharivorans [ID:02680]",
                    "vOTU497(Clostridium, Anaerotruncus, Ruminococcus, Blautia, Dorea, Coprococcus)", "Dorea longicatena [ID:03693]", "Ruminococcaceae species incertae sedis [ID:12239]", "Blautia massiliensis [ID:03342]", 
                    "Clostridiales species incertae sedis [ID:12476]", "Clostridiales species incertae sedis [ID:12554]",
                    "vOTU331(Faecalibacterium)", "Faecalibacterium prausnitzii [ID:06112]") -> SparData

#Start again
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
                    thresh = 0.2,#not related
                    measurePar = list(nlambda=10, 
                                      rep.num=10),　
                    normMethod = "none",
                    sparsMethod = "t-test", #sparsMethod = "threshhold", thresh = 0.01
                    adjust = "fdr",
                    alpha = 0.1, #fdr
                    lfdrThresh = 0.1, #not related
                    sampleSize = 540,  #注意して。UC111, CDでは31, Contは540 
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
          #nodeFilter = "highestConnect",
          #nodeFilterPar = 100,
          edgeInvisFilter = "threshold",
          edgeInvisPar = 0.2, 
          nodeSize = "eigenvector",
          repulsion =  0.8, #Fullは0.6
          rmSingles = "all", 
          labelScale = FALSE,
          cexLabels = 1.3,  
          nodeSizeSpread = 3,
          cexNodes = 2.5,
          nodeColor = "feature", 
          featVecCol = Bio, 
          colorVec =  biomecol,
          posCol = "darkturquoise", 
          negCol = "orange",
          edgeTranspLow = 0,
          edgeTranspHigh = 50,
          hubBorderCol = "darkgray",
          title1 = "Network of UC-related species in controls", #CD (n = 31) Control (n = 540) 
          showTitle = TRUE,
          cexTitle = 1.5)

biomecol_transp <- colToTransp(biomecol, 60)

legend(-1.1, 1.1, cex = 0.8, pt.cex = 2.5, title = "Biome:", 
       legend=levels(Bio), col = biomecol_transp, bty = "n", pch = 16) 

legend(0.7, 1.1, cex = 1, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("darkturquoise","orange"), 
       bty = "n", horiz = TRUE)


#CD analysis (Figure 5c)
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
left_join(NO, vOTU_input, by = "NO") -> vOTU_input #Select Cont
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
left_join(NO, F_input, by = "NO") -> F_input #Select Cont
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

#Can skip here
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
                    "vOTU533(Roseburia, Coprococcus, Unclassified Firmicutes, Clostridium, Ruminococcus)", "[Ruminococcus] torques [ID:03659]", "Clostridiales sp. [ID:03658]","uncultured Flavonifractor sp. [ID:07315]") -> SparData

#Restart here
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

# SparData <- na.omit(SparData) #for China data

biomecol <- c("#ff4b00", "#77d9a8", "dodgerblue3")

net <- netConstruct(SparData, #relative abundance
                    #filtTax = "numbSamp",
                    #filtTaxPar = list("numbSamp" = 180), 
                    measure = "spearman",
                    thresh = 0.2,# not related
                    measurePar = list(nlambda=10, 
                                      rep.num=10),　
                    normMethod = "none",
                    sparsMethod = "t-test", #sparsMethod = "threshhold", thresh = 0.01
                    adjust = "fdr",
                    alpha = 0.1, #fdr
                    lfdrThresh = 0.1, #not related
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
          #nodeFilter = "highestConnect",
          #nodeFilterPar = 100,
          edgeInvisFilter = "threshold",
          edgeInvisPar = 0.2, 
          nodeSize = "eigenvector",
          repulsion =  0.6, #CD fullも0.6
          rmSingles = "all",
          labelScale = FALSE,
          cexLabels = 1.1,  #CDだと0.6
          nodeSizeSpread = 3,
          cexNodes = 2.5,
          nodeColor = "feature", 
          featVecCol = Bio, 
          colorVec =  biomecol,
          posCol = "darkturquoise", 
          negCol = "orange",
          edgeTranspLow = 0,
          edgeTranspHigh = 70,
          hubBorderCol = "darkgray",
          title1 = "Network of CD-related species in controls",
          showTitle = TRUE,
          cexTitle = 1.5)

biomecol_transp <- colToTransp(biomecol, 60)

legend(-1, 1.1, cex = 1, pt.cex = 2.5, title = "Biome:", 
       legend=levels(Bio), col = biomecol_transp, bty = "n", pch = 16) 

legend(0.7, 1.1, cex = 1, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("darkturquoise","orange"), 
       bty = "n", horiz = TRUE)
