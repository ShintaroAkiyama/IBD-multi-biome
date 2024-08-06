library(dplyr)
library(tidyverse)
library(ggplot2)
library("phyloseq")
library(vegan)
library(ggstatsplot)

#Diversity analysis
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/SP_mOTU3_JP") 
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta <- rename(Meta, Group = UC_CD)
#rownames_to_column(Meta, "METAF") -> Meta #Only for alpha diversity, not for beta diversity
METADATA=sample_data(Meta, errorIfNULL = TRUE)

#alpa-diversity for bacterial species (Figure 1a, Supplementary Fig 1a-1b)
#Shannon
SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("ID", "IBD", "UC", "CD"))] 
DATA=otu_table(SP_input2, taxa_are_rows = FALSE)

Shannon_SP<-diversity(SP_input2,index="shannon",base=2) %>% as.data.frame()  #Select simpson, invsimpson
colnames(Shannon_SP)[1] <- "Shannon"
rownames_to_column(Shannon_SP, "METAF") -> Shannon_SP

Shannon_SP <- full_join(Shannon_SP, Meta, by = "METAF")

ggplot(Shannon_SP, aes(Group, Shannon, fill = Group)) + 
  geom_violin(aes(group = Group), width = 0.5) +
  scale_fill_manual(values = c("royalblue","#84919e","red"))+
  geom_boxplot(aes(group = Group), fill = "white", width = 0.1) + 
  theme_classic() +
  labs(x = "", y = "Shannon for species") +
  scale_y_continuous(limits = c(0, 8.5)) +
  scale_x_discrete(limit=c("Controls","UC","CD"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 22),
        text = element_text(size = 22)) +
  annotate("text", x = 1.5, y = 8, label = "P < 0.001", size = 6, fontface = "plain") +
  geom_segment(x = 1, xend = 2, y = 7.8, yend = 7.8) +
  annotate("text", x = 2, y = 8.5, label = "P < 0.001", size = 6, fontface = "plain")+
  geom_segment(x = 1, xend = 3, y = 8.3, yend = 8.3)+
  annotate("text", x = 2.5, y = 7.7, label = "P = 0.007", size = 6, fontface = "plain")+
  geom_segment(x = 2, xend = 3, y = 7.5, yend = 7.5)+
  guides(fill=guide_legend(title=NULL))

#Simpson
SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("ID", "IBD", "UC", "CD"))] 
DATA=otu_table(SP_input2, taxa_are_rows = FALSE)

Simpson_SP<-diversity(SP_input2,index="simpson",base=2) %>% as.data.frame()  #Select simpson, invsimpson
colnames(Simpson_SP)[1] <- "Simpson"
rownames_to_column(Simpson_SP, "METAF") -> Simpson_SP

Simpson_SP <- full_join(Simpson_SP, Meta, by = "METAF")

ggplot(Simpson_SP, aes(Group, Simpson, fill = Group)) + 
  geom_violin(aes(group = Group), width = 0.5) +
  scale_fill_manual(values = c("royalblue","#84919e","red"))+
  geom_boxplot(aes(group = Group), fill = "white", width = 0.1) + 
  theme_classic() +
  labs(x = "", y = "Simpson for species") +
  scale_y_continuous(limits = c(0, 1.3)) +
  scale_x_discrete(limit=c("Controls","UC","CD"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 18),
        text = element_text(size = 18)) +
  annotate("text", x = 1.5, y = 1.15, label = "P < 0.001", size = 6, fontface = "plain") +
  geom_segment(x = 1, xend = 2, y = 1.1, yend = 1.1) +
  annotate("text", x = 2, y = 1.25, label = "P < 0.001", size = 6, fontface = "plain")+
  geom_segment(x = 1, xend = 3, y = 1.2, yend = 1.2)+
  annotate("text", x = 2.5, y = 1.1, label = "P = 0.004", size = 6, fontface = "plain")+
  geom_segment(x = 2, xend = 3, y = 1.05, yend = 1.05)+
  guides(fill=guide_legend(title=NULL))

#Inverse Simpson
SP_input <- read.csv("SP_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("ID", "IBD", "UC", "CD"))] 
DATA=otu_table(SP_input2, taxa_are_rows = FALSE)

Simpson_SP<-diversity(SP_input2,index="invsimpson",base=2) %>% as.data.frame()  #Select simpson, invsimpson
colnames(Simpson_SP)[1] <- "Simpson"
rownames_to_column(Simpson_SP, "METAF") -> Simpson_SP

Simpson_SP <- full_join(Simpson_SP, Meta, by = "METAF")

ggplot(Simpson_SP, aes(Group, Simpson, fill = Group)) + 
  geom_violin(aes(group = Group), width = 0.5) +
  scale_fill_manual(values = c("royalblue","#84919e","red"))+
  geom_boxplot(aes(group = Group), fill = "white", width = 0.1) + 
  theme_classic() +
  labs(x = "", y = "Inverse Simpson for species") +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_discrete(limit=c("Controls","UC","CD"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 18),
        text = element_text(size = 18)) +
  annotate("text", x = 1.5, y = 79, label = "P < 0.001", size = 6, fontface = "plain") +
  geom_segment(x = 1, xend = 2, y = 76, yend = 76) +
  annotate("text", x = 2, y = 85, label = "P < 0.001", size = 6, fontface = "plain")+
  geom_segment(x = 1, xend = 3, y = 82, yend = 82)+
  annotate("text", x = 2.5, y = 75, label = "P = 0.004", size = 6, fontface = "plain")+
  geom_segment(x = 2, xend = 3, y = 72, yend = 72)+
  guides(fill=guide_legend(title=NULL))

#beta-diversity (Figure 1b)
adonis2(SP_input2~Group, data=Meta, permutation=9999, method="bray") #obtain p value for PERMANOVA

data_phylo<-phyloseq(DATA, METADATA)
bc=ordinate(data_phylo, method = "MDS", distane = "bray") #NMDSなど選べる
plot_ordination(data_phylo, bc, color= "Group")+
  geom_point(size=3)+
  theme_classic() +
  stat_ellipse() +
  labs(x = "MDS1", y = "MDS2") +
  scale_color_manual(values = c("royalblue","#84919e","red"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 24),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 24),
        axis.line = element_line(colour = "black", linewidth= 0.7, linetype = "solid"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 28))+
  guides(fill=guide_legend(title=NULL))+
  annotate("text", x = 0.4, y = -0.42, label = "P = 1e-04", size = 8, fontface = "plain")

#Richness (Supplementary Fig 1c)
SP_input <- read.csv("SP_IBD_num.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) 
SP_input2 <- SP_input[, -which (colnames(SP_input) %in% c("ID", "IBD", "UC", "CD"))] 
rownames_to_column(SP_input2, "METAF") -> SP_input2

Richness_SP <- full_join(SP_input2, Meta, by = "METAF")

ggplot(Richness_SP, aes(Group, Richness, fill = Group)) + 
  geom_violin(aes(group = Group), width = 0.5) +
  scale_fill_manual(values = c("royalblue","#84919e","red"))+
  geom_boxplot(aes(group = Group), fill = "white", width = 0.1) + 
  theme_classic() +
  labs(x = "", y = "Richness for species") +
  scale_y_continuous(limits = c(0,  600)) +
  scale_x_discrete(limit=c("Controls","UC","CD"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 18),
        text = element_text(size = 18)) +
  annotate("text", x = 1.5, y = 570, label = "P < 0.001", size = 6, fontface = "plain") +
  geom_segment(x = 1, xend = 2, y = 555, yend = 555) +
  annotate("text", x = 2, y = 600, label = "P < 0.001", size = 6, fontface = "plain")+
  geom_segment(x = 1, xend = 3, y = 585, yend = 585)+
  annotate("text", x = 2.5, y = 540, label = "P = 0.007", size = 6, fontface = "plain")+
  geom_segment(x = 2, xend = 3, y = 525, yend = 525)+
  guides(fill=guide_legend(title=NULL))

#Diversity for ARG (Figure 2a, Supplementary Fig 1d-1e)
#Shannon
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Figure用解析/Bacteriome_all_results_for_pval/Diversity/Microbiome_diversity")
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta <- rename(Meta, Group = UC_CD)
METADATA=sample_data(Meta, errorIfNULL = TRUE)

ARG_input <- read.csv("ARG_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data
ARG_input2 <- ARG_input[, -which (colnames(ARG_input) %in% c("IBD_HC2", "IBD_1", "UC", "Crohn"))] 
DATA=otu_table(ARG_input2, taxa_are_rows = FALSE)

Shannon_ARG<-diversity(ARG_input2,index="shannon",base=2) %>% as.data.frame() 
colnames(Shannon_ARG)[1] <- "Shannon"
Shannon_ARG <- cbind(Shannon_ARG, Meta)

ggplot(Shannon_ARG, aes(Group, Shannon, fill = Group)) + 
  geom_violin(aes(group = Group), width = 0.5) +
  scale_fill_manual(values = c("royalblue","#84919e","red"))+
  geom_boxplot(aes(group = Group), fill = "white", width = 0.1) + 
  theme_classic() +
  labs(x = "", y = "Shannon for ARG") +
  scale_y_continuous(limits = c(0, 8)) +
  scale_x_discrete(limit=c("Controls","UC","CD"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 22),
        text = element_text(size = 22)) +
  guides(fill=guide_legend(title=NULL))

#simpson
ARG_input <- read.csv("ARG_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data
ARG_input2 <- ARG_input[, -which (colnames(ARG_input) %in% c("IBD_HC2", "IBD_1", "UC", "Crohn"))] 
DATA=otu_table(ARG_input2, taxa_are_rows = FALSE)

Simpson_ARG<-diversity(ARG_input2,index="simpson",base=2) %>% as.data.frame() 
colnames(Simpson_ARG)[1] <- "Simpson"
Simpson_ARG <- cbind(Simpson_ARG, Meta)

ggplot(Simpson_ARG, aes(Group, Simpson, fill = Group)) + 
  geom_violin(aes(group = Group), width = 0.5) +
  scale_fill_manual(values = c("royalblue","#84919e","red"))+
  geom_boxplot(aes(group = Group), fill = "white", width = 0.1) + 
  theme_classic() +
  labs(x = "", y = "Simpson for ARG") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(limit=c("Controls","UC","CD"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 18),
        text = element_text(size = 18)) +
  guides(fill=guide_legend(title=NULL))

#Inverse simpson
ARG_input <- read.csv("ARG_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data
ARG_input2 <- ARG_input[, -which (colnames(ARG_input) %in% c("IBD_HC2", "IBD_1", "UC", "Crohn"))] 
DATA=otu_table(ARG_input2, taxa_are_rows = FALSE)

Simpson_ARG<-diversity(ARG_input2,index="invsimpson",base=2) %>% as.data.frame() 
colnames(Simpson_ARG)[1] <- "Simpson"
Simpson_ARG <- cbind(Simpson_ARG, Meta)

ggplot(Simpson_ARG, aes(Group, Simpson, fill = Group)) + 
  geom_violin(aes(group = Group), width = 0.5) +
  scale_fill_manual(values = c("royalblue","#84919e","red"))+
  geom_boxplot(aes(group = Group), fill = "white", width = 0.1) + 
  theme_classic() +
  labs(x = "", y = "Inverse Simpson for ARG") +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_discrete(limit=c("Controls","UC","CD"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 18),
        text = element_text(size = 18)) +
  guides(fill=guide_legend(title=NULL))

#beta-diversity (Figure 2b)
adonis2(ARG_input2~Group, data=Meta, permutation=9999, method="bray") #obtain p value for PERMANOVA 

data_phylo<-phyloseq(DATA, METADATA)
bc=ordinate(data_phylo, method = "MDS", distane = "bray") 
plot_ordination(data_phylo, bc, color= "Group")+
  geom_point(size=3)+
  theme_classic() +
  stat_ellipse() +
  labs(x = "MDS1", y = "MDS2") +
  scale_color_manual(values = c("royalblue","#84919e","red"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 19),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 19),
        axis.line = element_line(colour = "black", linewidth= 0.7, linetype = "solid"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 28))+
  guides(fill=guide_legend(title=NULL))+
  annotate("text", x = 0.4, y = -0.35, label = "P = 1e-04", size = 8, fontface = "plain")

#Richness(Supplementary Fig 1f)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Figure用解析/Bacteriome_all_results_for_pval/Diversity/Microbiome_diversity")
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta <- rename(Meta, Group = UC_CD)
rownames_to_column(Meta, "ID") -> Meta

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/ARG")
SP_input <- read.csv("ARG.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) 
rownames_to_column(SP_input, "ID") -> SP_input2

Richness_SP <- full_join(SP_input2, Meta, by = "ID")

ggplot(Richness_SP, aes(Group, Number_of_ARG, fill = Group)) + 
  geom_violin(aes(group = Group), width = 0.5) +
  scale_fill_manual(values = c("royalblue","#84919e","red"))+
  geom_boxplot(aes(group = Group), fill = "white", width = 0.1) + 
  theme_classic() +
  labs(x = "", y = "Richness for ARG") +
  scale_y_continuous(limits = c(0,  200)) +
  scale_x_discrete(limit=c("Controls","UC","CD"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 18),
        text = element_text(size = 18)) +
  #annotate("text", x = 1.5, y = 570, label = "P < 0.001", size = 6, fontface = "plain") +
  #geom_segment(x = 1, xend = 2, y = 555, yend = 555) +
  #annotate("text", x = 2, y = 600, label = "P < 0.001", size = 6, fontface = "plain")+
  #geom_segment(x = 1, xend = 3, y = 585, yend = 585)+
  annotate("text", x = 2.5, y = 175, label = "P = 0.02", size = 6, fontface = "plain")+
  geom_segment(x = 2, xend = 3, y = 165, yend = 165)+
  guides(fill=guide_legend(title=NULL))

#Diversity for vOTU (Figure 3a, Supplementary Fig 1g-1h)
#Shannon
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Figure用解析/Bacteriome_all_results_for_pval/Diversity/Microbiome_diversity")
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta <- rename(Meta, Group = UC_CD)
METADATA=sample_data(Meta, errorIfNULL = TRUE)

vOTU_input <- read.csv("vOTU_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data
vOTU_input2 <- vOTU_input[, -which (colnames(vOTU_input) %in% c("IBD_HC2", "UC", "Crohn"))] 
DATA=otu_table(vOTU_input2, taxa_are_rows = FALSE)

Shannon_vOTU<-diversity(vOTU_input2,index="shannon",base=2) %>% as.data.frame() 
colnames(Shannon_vOTU)[1] <- "Shannon"
Shannon_vOTU <- cbind(Shannon_vOTU, Meta)

ggplot(Shannon_vOTU, aes(Group, Shannon, fill = Group)) + 
  geom_violin(aes(group = Group), width = 0.5) +
  scale_fill_manual(values = c("royalblue","#84919e","red"))+
  geom_boxplot(aes(group = Group), fill = "white", width = 0.1) + 
  theme_classic() +
  labs(x = "", y = "Shannon for vOTU") +
  scale_y_continuous(limits = c(0, 8.5)) +
  scale_x_discrete(limit=c("Controls","UC","CD"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 22),
        text = element_text(size = 22)) +
  annotate("text", x = 1.5, y = 7.7, label = "P < 0.001", size = 6, fontface = "plain") +
  geom_segment(x = 1, xend = 2, y = 7.5, yend = 7.5) +
  annotate("text", x = 2, y = 8.2, label = "P < 0.001", size = 6, fontface = "plain")+
  geom_segment(x = 1, xend = 3, y = 8, yend = 8)+
  annotate("text", x = 2.5, y = 7.4, label = "P = 0.02", size = 6, fontface = "plain")+
  geom_segment(x = 2, xend = 3, y = 7.2, yend = 7.2)+
  guides(fill=guide_legend(title=NULL))

#simpson
vOTU_input <- read.csv("vOTU_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data
vOTU_input2 <- vOTU_input[, -which (colnames(vOTU_input) %in% c("IBD_HC2", "UC", "Crohn"))] 
DATA=otu_table(vOTU_input2, taxa_are_rows = FALSE)

Simpson_vOTU<-diversity(vOTU_input2,index="simpson",base=2) %>% as.data.frame() 
colnames(Simpson_vOTU)[1] <- "Simpson"
Simpson_vOTU <- cbind(Simpson_vOTU, Meta)

ggplot(Simpson_vOTU, aes(Group, Simpson, fill = Group)) + 
  geom_violin(aes(group = Group), width = 0.5) +
  scale_fill_manual(values = c("royalblue","#84919e","red"))+
  geom_boxplot(aes(group = Group), fill = "white", width = 0.1) + 
  theme_classic() +
  labs(x = "", y = "Simpson for vOTU") +
  scale_y_continuous(limits = c(0, 1.3)) +
  scale_x_discrete(limit=c("Controls","UC","CD"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 18),
        text = element_text(size = 18)) +
  annotate("text", x = 1.5, y = 1.15, label = "P < 0.001", size = 6, fontface = "plain") +
  geom_segment(x = 1, xend = 2, y = 1.1, yend = 1.1) +
  annotate("text", x = 2, y = 1.25, label = "P < 0.001", size = 6, fontface = "plain")+
  geom_segment(x = 1, xend = 3, y = 1.2, yend = 1.2)+
  annotate("text", x = 2.5, y = 1.1, label = "P = 0.02", size = 6, fontface = "plain")+
  geom_segment(x = 2, xend = 3, y = 1.05, yend = 1.05)+
  guides(fill=guide_legend(title=NULL))

#Inverse Simpson
vOTU_input <- read.csv("vOTU_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data
vOTU_input2 <- vOTU_input[, -which (colnames(vOTU_input) %in% c("IBD_HC2", "UC", "Crohn"))] 
DATA=otu_table(vOTU_input2, taxa_are_rows = FALSE)

Simpson_vOTU<-diversity(vOTU_input2,index="invsimpson",base=2) %>% as.data.frame() 
colnames(Simpson_vOTU)[1] <- "Simpson"
Simpson_vOTU <- cbind(Simpson_vOTU, Meta)

ggplot(Simpson_vOTU, aes(Group, Simpson, fill = Group)) + 
  geom_violin(aes(group = Group), width = 0.5) +
  scale_fill_manual(values = c("royalblue","#84919e","red"))+
  geom_boxplot(aes(group = Group), fill = "white", width = 0.1) + 
  theme_classic() +
  labs(x = "", y = "Inverse Simpson for vOTU") +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_discrete(limit=c("Controls","UC","CD"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 18),
        text = element_text(size = 18)) +
  annotate("text", x = 1.5, y = 79, label = "P < 0.001", size = 6, fontface = "plain") +
  geom_segment(x = 1, xend = 2, y = 76, yend = 76) +
  annotate("text", x = 2, y = 85, label = "P < 0.001", size = 6, fontface = "plain")+
  geom_segment(x = 1, xend = 3, y = 82, yend = 82)+
  annotate("text", x = 2.5, y = 75, label = "P = 0.02", size = 6, fontface = "plain")+
  geom_segment(x = 2, xend = 3, y = 72, yend = 72)+
  guides(fill=guide_legend(title=NULL))

#beta-diversity (Figure 4b)
adonis2(vOTU_input2~Group, data=Meta, permutation=9999, method="bray") 

data_phylo<-phyloseq(DATA, METADATA)
bc=ordinate(data_phylo, method = "MDS", distane = "bray") 
plot_ordination(data_phylo, bc, color= "Group")+
  geom_point(size=3)+
  theme_classic() +
  stat_ellipse() +
  labs(x = "MDS1", y = "MDS2") +
  scale_color_manual(values = c("royalblue","#84919e","red"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.line = element_line(colour = "black", linewidth= 0.7, linetype = "solid"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 28))+
  guides(fill=guide_legend(title=NULL))+
  annotate("text", x = 0.4, y = -0.35, label = "P = 1e-04", size = 8, fontface = "plain")

#Richness (Supplementary Fig 1i)
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Figure用解析/Bacteriome_all_results_for_pval/Diversity/Microbiome_diversity")
Meta <- read.csv("Meta.csv", header = TRUE, na.strings = c(NA, ''), row.names=1)
Meta <- rename(Meta, Group = UC_CD)
rownames_to_column(Meta, "ID") -> Meta

setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Manuscript/Nat Com/Revised/MaAsLin/Phage")
SP_input <- read.csv("vOTU_num.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) 
rownames_to_column(SP_input, "ID") -> SP_input2

Richness_SP <- full_join(SP_input2, Meta, by = "ID")

ggplot(Richness_SP, aes(Group, Richness, fill = Group)) + 
  geom_violin(aes(group = Group), width = 0.5) +
  scale_fill_manual(values = c("royalblue","#84919e","red"))+
  geom_boxplot(aes(group = Group), fill = "white", width = 0.1) + 
  theme_classic() +
  labs(x = "", y = "Richness for vOTU") +
  scale_y_continuous(limits = c(0,  280)) +
  scale_x_discrete(limit=c("Controls","UC","CD"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 18),
        text = element_text(size = 18)) +
  annotate("text", x = 1.5, y = 230, label = "P < 0.001", size = 6, fontface = "plain") +
  geom_segment(x = 1, xend = 2, y = 220, yend = 220) +
  annotate("text", x = 2, y = 260, label = "P < 0.001", size = 6, fontface = "plain")+
  geom_segment(x = 1, xend = 3, y = 250, yend = 250)+
  #annotate("text", x = 2.5, y = 180, label = "P = 0.02", size = 6, fontface = "plain")+
  #geom_segment(x = 2, xend = 3, y = 165, yend = 165)+
  guides(fill=guide_legend(title=NULL))
