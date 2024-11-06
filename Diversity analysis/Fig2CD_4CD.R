library(ggstatsplot)
library(dplyr)
library(ggplot2)

#Figure 2c, 2d
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Figure用解析/Bacteriome_all_results_for_pval/ARG_Pval/Speceis vs ARG spearman")
#ARG number ans sum 
ARG_input <- read.csv("ARG_IBD.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data 
ARG_input  <- as.data.frame(t(ARG_input)) 
ARG_input %>% filter(ARG_input$IBD_1 == 0 | ARG_input$UC == 1 | ARG_input$Crohn == 1) -> ARG_compare

ARG_compare %>% mutate(Group = case_when(ARG_compare$"UC"==1  ~ "UC",  ARG_compare$"Crohn" == 1  ~ "CD", ARG_compare$"IBD_1" == 0 ~ "Controls")) -> ARG_compare 

#rename(.data = ARG_compare, "Number of ARG" = "Number_of_ARG") -> data
#rename(.data = data, "Sum of ARG abundance (RPKM)" = "Sum_of_ARG") -> data

#Number of ARG
ggplot(ARG_compare, aes(Group, Number_of_ARG, fill = Group)) + 
  geom_violin(aes(group = Group), width = 0.5) +
  scale_fill_manual(values = c("royalblue","#84919e","red"))+
  geom_boxplot(aes(group = Group), fill = "white", width = 0.1) + 
  theme_classic() +
  labs(x = "", y = "Number of ARG") +
  scale_y_continuous(limits = c(0, 200)) +
  scale_x_discrete(limit=c("Controls","UC","CD"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 22),
        text = element_text(size = 22)) +
  annotate("text", x = 2, y = 185, label = "P = 0.05", size = 6, fontface = "plain") +
  geom_segment(x = 1, xend = 3, y = 180, yend = 180) +
  annotate("text", x = 2.5, y = 175, label = "P = 0.02", size = 6, fontface = "plain")+
  geom_segment(x = 2, xend = 3, y = 170, yend = 170)+
  guides(fill=guide_legend(title=NULL))

#Sum of ARG
ggplot(ARG_compare, aes(Group, Sum_of_ARG, fill = Group)) + 
  geom_violin(aes(group = Group), width = 0.5) +
  scale_fill_manual(values = c("royalblue","#84919e","red"))+
  geom_boxplot(aes(group = Group), fill = "white", width = 0.1) + 
  theme_classic() +
  labs(x = "", y = "Sum of ARG abundance (RPKM)") +
  scale_y_continuous(limits = c(0, 8000)) +
  scale_x_discrete(limit=c("Controls","UC","CD"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 22),
        text = element_text(size = 22)) +
  annotate("text", x = 2, y = 7950, label = "P = 0.0049", size = 6, fontface = "plain") +
  geom_segment(x = 1, xend = 3, y = 7750, yend = 7750) +
  guides(fill=guide_legend(title=NULL))

#Check Pval
P_num <- ggbetweenstats(plot.type = "box", data = ARG_compare, x = Group, y = Number_of_ARG, type = "nonparametirc", title = "", 
                        centrality.plotting = FALSE, ggplot.component = list(theme(text=element_text(size=20)), scale_x_discrete(limits=c("Controls", "UC", "CD"))), results.subtitle = FALSE, centrality.point.args = list(size = 4, color = "darkred"), centrality.label.args = list(size = 4, nudge_x = 0.5, segment.linetype = 2, min.segment.length = 0), ggsignif.args = list(textsize = 6, tip_length = 0.01))
P_sum <- ggbetweenstats(plot.type = "box", data = ARG_compare, x = Group, y = Sum_of_ARG, type = "nonparametirc", title = "", 
                        centrality.plotting = FALSE, ggplot.component = list(theme(text=element_text(size=20)), scale_x_discrete(limits=c("Controls", "UC", "CD"))), results.subtitle = FALSE, centrality.point.args = list(size = 4, color = "darkred"), centrality.label.args = list(size = 4, nudge_x = 0.5, segment.linetype = 2, min.segment.length = 0), ggsignif.args = list(textsize = 6, tip_length = 0.01))

combine_plots(
  list(P_num, P_sum),
  plotgrid.args = list(nrow = 1),
  annotation.args = list(
    title = "",
    caption = ""
  )
)

#Figure 4c, 4d
setwd("/Users/akiyama/Documents/筑波大学/筑波大学研究/プロジェクト/Microbiome共同研究/Figure用解析/vOTU_heatmap (pval)/vOTU_heatmap/解析用CSVファイル（上書き禁止）")
vOTU_input <- read.csv("vOTU_heatmap_update_abundance_lyso.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data 
colnames(vOTU_input)[5] <- "Sum_of_abundance_Lysogenic"

#vOTU_input <- read.csv("Generalist.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data 
#colnames(vOTU_input)[5] <- "Sum_of_abundance_Generalist"

#Sum of abundance (Lysogenic)
ggplot(vOTU_input, aes(Group, Sum_of_abundance_Lysogenic, fill = Group)) + 
  geom_violin(aes(group = Group), width = 0.5) +
  scale_fill_manual(values = c("royalblue","#84919e","red"))+
  geom_boxplot(aes(group = Group), fill = "white", width = 0.1) + 
  theme_classic() +
  labs(x = "", y = "Sum of abundance (Lysogenic)") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(limit=c("Controls","UC","CD"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 22),
        text = element_text(size = 22)) +
        guides(fill=guide_legend(title=NULL))

#Sum of abundance (Specialist)
vOTU_input <- read.csv("Specialist.csv", header = TRUE, na.strings = c(NA, ''), row.names=1) #abundance data 
colnames(vOTU_input)[5] <- "Sum_of_abundance_Specialist"

ggplot(vOTU_input, aes(Group, Sum_of_abundance_Specialist, fill = Group)) + 
  geom_violin(aes(group = Group), width = 0.5) +
  scale_fill_manual(values = c("royalblue","#84919e","red"))+
  geom_boxplot(aes(group = Group), fill = "white", width = 0.1) + 
  theme_classic() +
  labs(x = "", y = "Sum of abundance (Specialist)") +
  scale_y_continuous(limits = c(0, 1.13)) +
  scale_x_discrete(limit=c("Controls","UC","CD"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 22),
        text = element_text(size = 22)) +
  annotate("text", x = 2, y = 1.12, label = "P = 0.03", size = 6, fontface = "plain") +
  geom_segment(x = 1, xend = 3, y = 1.08, yend = 1.08) +
  annotate("text", x = 1.5, y = 1.04, label = "P = 0.02", size = 6, fontface = "plain")+
  geom_segment(x = 1, xend = 2, y = 1, yend = 1)+
  guides(fill=guide_legend(title=NULL))

Lyso <- ggbetweenstats(plot.type = "box", data = vOTU_input, x = Group, y = Sum_of_abundance_Lysogenic, type = "nonparametirc", title = "", centrality.plotting = FALSE, ggplot.component = list(theme(text=element_text(size=20)), scale_x_discrete(limits=c("Controls", "UC", "CD"))), results.subtitle = FALSE, centrality.point.args = list(size = 4, color = "darkred"), centrality.label.args = list(size = 4, nudge_x = 0.5, segment.linetype = 2, min.segment.length = 0), ggsignif.args = list(textsize = 6, tip_length = 0.01))
#Gen<-ggbetweenstats(plot.type = "box", data = vOTU_input, x = Group, y = Sum_of_abundance_Generalist, type = "nonparametirc", title = "", centrality.plotting = FALSE, ggplot.component = list(theme(text=element_text(size=20)), scale_x_discrete(limits=c("Controls", "UC", "CD"))), results.subtitle = FALSE, centrality.point.args = list(size = 4, color = "darkred"), centrality.label.args = list(size = 4, nudge_x = 0.5, segment.linetype = 2, min.segment.length = 0), ggsignif.args = list(textsize = 6, tip_length = 0.01))
Spe<-ggbetweenstats(plot.type = "box", data = vOTU_input, x = Group, y = Sum_of_abundance_Specialist, type = "nonparametirc", title = "", centrality.plotting = FALSE, ggplot.component = list(theme(text=element_text(size=20)), scale_x_discrete(limits=c("Controls", "UC", "CD"))), results.subtitle = FALSE, centrality.point.args = list(size = 4, color = "darkred"), centrality.label.args = list(size = 4, nudge_x = 0.5, segment.linetype = 2, min.segment.length = 0), ggsignif.args = list(textsize = 6, tip_length = 0.01))

combine_plots(
  list(Lyso, Spe),
  plotgrid.args = list(nrow = 1),
  annotation.args = list(
    title = "",
    caption = ""
  )
)
