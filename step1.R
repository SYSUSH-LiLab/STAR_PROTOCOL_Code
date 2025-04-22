## ========= 1. Load the required libraries ==============
library(tidyverse)
library(ggplot2)
library(survival)
library(broom)
library(survminer)
library(randomForest)
library(RTN)
library(glmnet)

## ======= 2. Aging-associated gene with LASSO ================
load("Datasets/TCGA-STAD.rda")
tmp <- gc.expr[aging.gene$Symbol,] %>% t() %>% as.data.frame() %>% bind_cols(gc.subtype)
age.os.cox <- lapply(aging.gene$Symbol, function(x){
  formu <- sprintf("coxph(Surv(OS.time, OS) ~ %s, data = tmp)", x)
  res <- eval(parse(text = formu)) %>% tidy(exponentiate = T)
}) %>% do.call(rbind, .)
train.gene <- age.os.cox$term[age.os.cox$p.value < 0.1]


## ==== 3. LASSO to select most important genes ========
lasso.gene <- lapply(1:100, function(x){
  fit <- cv.glmnet(gc.expr[train.gene,] %>% t(), 
                           as.matrix(Surv(gc.subtype$OS.time,gc.subtype$OS)),
                           family = "cox", alpha=1, nfolds = 5)
  efs <- coef(fit, s=fit$lambda.min) %>% as.matrix()
  gene.lasso <- rownames(efs)[efs[,1] != 0]
  return(train.gene %in% gene.lasso)
}) %>% do.call(cbind, .) %>% magrittr::set_rownames(train.gene) %>% rowSums()
final.gene <- names(lasso.gene)[lasso.gene > 85]

## ====== 4 Aging-associated Index (AAI) =============
asi.expr <- gc.expr[final.gene,] %>% t() %>% as.data.frame() %>% bind_cols(gc.subtype[,c('OS.time','OS')])
asi.fit <- coxph(Surv(OS.time, OS) ~ ., data = asi.expr)
gc.subtype$asi.score <- sapply(gc.expr[final.gene,], function(x){
  sum(x * coef(asi.fit)) %>% return()
})
coxph(Surv(OS.time, OS) ~ asi.score, data = gc.subtype) %>% broom::tidy(exponentiate = T)
coxph(Surv(DFI.time, DFI) ~ asi.score, data = gc.subtype) %>% broom::tidy(exponentiate = T)

gc.subtype$ASI.HL <- ifelse(gc.subtype$asi.score > median(gc.subtype$asi.score), 
                            'ASI-H', 'ASI-L')

survfit(Surv(OS.time, OS) ~ ASI.HL, data=gc.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='OS (%)',
    legend.title = '',legend.lab=c('ASI-H', 'ASI-L'),
    ggtheme = theme_classic2(),
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12
  )

survfit(Surv(DFI.time, DFI) ~ ASI.HL, data=gc.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='DFS (%)',
    legend.title = '',legend.lab=c('ASI-H', 'ASI-L'),
    ggtheme = theme_classic2(),
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12
  )

## ====== 5. AAI validation =============
load("Datasets/GSE62254.rda")
identical(colnames(GSE62254.expr), GSE62254.subtype$GEO_ID)
GSE62254.subtype$asi.score <- sapply(GSE62254.expr[names(asi.fit$coefficients),], function(x){
  sum(x * coef(asi.fit))
})
coxph(Surv(OS.m, Death) ~ asi.score, data = GSE62254.subtype) %>% broom::tidy(exponentiate = T)
coxph(Surv(DFS.m, Recur) ~ asi.score, data = GSE62254.subtype) %>% broom::tidy(exponentiate = T)

GSE62254.subtype$ASI.HL <- ifelse(GSE62254.subtype$asi.score > median(GSE62254.subtype$asi.score), 'ASI-H','ASI-L')
survfit(Surv(OS.m, Death) ~ ASI.HL, data=GSE62254.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='OS (%)',
    legend.title = '',legend.lab=c('ASI-H', 'ASI-L'),
    ggtheme = theme_classic2(),
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12
  )

survfit(Surv(DFS.m, Recur) ~ ASI.HL, data=GSE62254.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='DFS (%)',
    legend.title = '',legend.lab=c('ASI-H', 'ASI-L'),
    ggtheme = theme_classic2(),
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12
  )

