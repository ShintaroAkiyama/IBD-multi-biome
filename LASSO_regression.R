## script to construct LASSO model to distinguish between IBD and controls
## Suguru Nishijima (nishjima.suguru@gmail.com)

## read library
library(tidyverse)
library(data.table)
library(pROC)
library(caret)


## setting
train_control <- trainControl(method = "repeatedcv", savePredictions = T, repeats = 5, classProbs = T, summaryFunction = twoClassSummary)
tune_grid.lasso <- expand.grid(alpha = 1, lambda = 10 ^ (1:10 * -1))


check_roc <- function(model, disease){
  best <- model$bestTune
  
  ## LASSO
  print("LASSO")    
  pred <- model$pred %>% filter(alpha == best$alpha & lambda == best$lambda)
  pred.score <- tapply(pred[, 5], pred$rowIndex, mean)        
  
  roc <- roc(disease, pred.score)
  return(roc)
}



## read data
d <- read.delim("data/feces_230402.motus3.sp.txt", header = T, row.names = 1, check.names = F)
colnames(d) <- colnames(d) %>% str_pad(width = 4, pad = "0")

md <- read.delim("data/Japanese_4D_IBD.txt", header = T)
md$ID <- md$ID %>% str_pad(width = 4, pad = "0")
md <- md[order(md$ID), ]

keep <- colnames(d) %in% md$ID
d <- d[, keep]

md$disease <- ifelse(md$IBD == 2, "Control", ifelse(md$CD == 1, "Crohn's disease", "Ulcerative colitis"))
d <- d %>% t() %>% data.frame(check.names = F)
colnames(d) <- colnames(d) %>% str_remove("ID\\:") %>% str_replace("Unassigned species", "-1")

## remove minor species
keep <- apply(d, 2, mean) > 1E-5 & apply(d>0, 2, mean) > 0.1
d_all <- d
d <- d[, keep]


## construct LASSO model
## LASSO (control vs IBD)
label <- ifelse(md$disease == "Control", "Control", "IBD")
d.min <- min(d[d != 0])/2

d.scale <- scale(log10(d + d.min))
  
model <- caret::train(d.scale, label, trControl = train_control, method = "glmnet", metric = "ROC", tuneGrid = tune_grid.lasso) 
  
## save model
saveRDS(model, file = "out/rds/LASSO_IBD.rds")
  
## plot AUC
res <- check_roc(model, label)  
df.roc <- data.frame(fpr = 1 - res$specificities, tpr = res$sensitivities)

auc <- res$auc %>% round(digits = 2)
auc %>% print()
  
roc <- ggplot(df.roc, aes(x = fpr, y = tpr)) + 
  theme_bw() +
  geom_path() +
  ggtitle("IBD") +
  xlab("False positive rate") +
  ylab("True positive rate") +
  #scale_color_manual(values = tidy_col[5:6]) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom") +
  annotate("text", x = 0.3, y = 0, label = paste0("AUC=", auc), hjust = 0, vjust = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "gray")
roc
  
  
## LASSO (control vs CD)
keep <- md$disease %>% str_detect("Ulcerative") 
d2 <- d[!keep, ]
md2 <- md[!keep, ]    

label <- ifelse(md2$disease == "Control", "Control", "CD")
d2.min <- min(d2[d2 != 0])/2

keep <- apply(d2 == 0, 2, sum) == nrow(d2)
d2 <- d2[, !keep]

d.scale <- scale(log10(d2 + d2.min))
model <- caret::train(d.scale, label, trControl = train_control, method = "glmnet", metric = "ROC", tuneGrid = tune_grid.lasso) 

## save model
saveRDS(model, file = "out/rds/LASSO_CD.rds")

## plot AUC
res <- check_roc(model, label)  
df.roc <- data.frame(fpr = 1 - res$specificities, tpr = res$sensitivities)

auc <- res$auc %>% round(digits = 2)
auc %>% print()

roc <- ggplot(df.roc, aes(x = fpr, y = tpr)) + 
  theme_bw() +
  geom_path() +
  ggtitle("CD") +
  xlab("False positive rate") +
  ylab("True positive rate") +
  #scale_color_manual(values = tidy_col[5:6]) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom") +
  annotate("text", x = 0.3, y = 0, label = paste0("AUC=", auc), hjust = 0, vjust = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "gray")
roc
    
  
## LASSO (control vs UC)
keep <- md$disease %>% str_detect("Crohn") 
d2 <- d[!keep, ]
md2 <- md[!keep, ]    
  
label <- ifelse(md2$disease == "Control", "Control", "UC")
d2.min <- min(d2[d2 != 0])/2
    
keep <- apply(d2 == 0, 2, sum) == nrow(d2)
d2 <- d2[, !keep]
    
d.scale <- scale(log10(d2 + d2.min))
model <- caret::train(d.scale, label, trControl = train_control, method = "glmnet", metric = "ROC", tuneGrid = tune_grid.lasso) 
    
## save model
saveRDS(model, file = "out/rds/LASSO_UC.rds")
    
## plot AUC
res <- check_roc(model, label)  
df.roc <- data.frame(fpr = 1 - res$specificities, tpr = res$sensitivities)
    
auc <- res$auc %>% round(digits = 2)
auc %>% print()
    
roc <- ggplot(df.roc, aes(x = fpr, y = tpr)) + 
  theme_bw() +
  geom_path() +
  ggtitle("UC") +
  xlab("False positive rate") +
  ylab("True positive rate") +
  #scale_color_manual(values = tidy_col[5:6]) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "bottom") +
  annotate("text", x = 0.3, y = 0, label = paste0("AUC=", auc), hjust = 0, vjust = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "gray")
roc
    




## feature importance
files <- list.files("out/rds/", pattern = "LASSO.*rds") %>% paste0("out/rds/", .)
files

i <- files[1]
for(i in files){
  print(i)
  out <- i %>% str_replace("\\.rds", ".feature_importance.tsv")
  
  model <- read_rds(i)
  
  beta <- stats::predict(model$finalModel, s = model$bestTune$lambda, type = "coef")
  res <- data.frame(Overall = beta[, 1]) %>% data.frame() 
  res$Overall <- res$Overall %>% abs()
  res <- res %>% arrange(-Overall) %>% rownames_to_column() %>% filter(rowname != "(Intercept)")

  out <- i %>% str_replace("\\.rds", ".feature_importance.tsv")
  write_tsv(res, file = out)
}
