# ======= Step 1 =============
## ========= 1.1 Load the required libraries ==============
library(tidyverse)
library(ggplot2)
library(survival)
library(broom)
library(survminer)
library(randomForest)
library(RTN)
library(glmnet)

## ======= 1.2 Aging-associated gene with LASSO ================
load("Datasets/TCGA-STAD.rda")
tmp <- gc.expr[aging.gene$Symbol,] %>% t() %>% as.data.frame() %>% bind_cols(gc.subtype)
age.os.cox <- lapply(aging.gene$Symbol, function(x){
  formu <- sprintf("coxph(Surv(OS.time, OS) ~ %s, data = tmp)", x)
  res <- eval(parse(text = formu)) %>% tidy(exponentiate = T)
}) %>% do.call(rbind, .)
train.gene <- age.os.cox$term[age.os.cox$p.value < 0.1]

## ==== 1.3 LASSO to select most important genes ========
lasso.gene <- lapply(1:100, function(x){
  fit <- cv.glmnet(gc.expr[train.gene,] %>% t(), 
                   as.matrix(Surv(gc.subtype$OS.time,gc.subtype$OS)),
                   family = "cox", alpha=1, nfolds = 5)
  efs <- coef(fit, s=fit$lambda.min) %>% as.matrix()
  gene.lasso <- rownames(efs)[efs[,1] != 0]
  return(train.gene %in% gene.lasso)
}) %>% do.call(cbind, .) %>% magrittr::set_rownames(train.gene) %>% rowSums()
final.gene <- names(lasso.gene)[lasso.gene > 85]

f1<- lasso.gene %>% as.data.frame() %>% rownames_to_column('term') %>%
  magrittr::set_colnames(c('term','freq')) %>% filter(freq > 0) %>% 
  left_join(age.os.cox[,c('term','estimate','p.value')]) %>% 
  dplyr::mutate(HR = ifelse(estimate > 1, 'Risk', 'Protective'), estimate=round(estimate, 1)) %>% 
  ggdotchart(x = "term", y = "freq", color='HR', dot.size = 8,
             sorting = "descending", add = "segments", palette = 'Dark1',                      
             ggtheme = theme_pubr(),label = 'estimate',
             font.label = list(color = "white", size = 9, 
                               vjust = 0.5))+
  geom_hline(yintercept=85, linetype='dashed', color='#FC4E07')+
  labs(x='', y='Occurance', color='Hazard ratio')
ggsave("Figure1.pdf", plot = f1, width = 10, height = 8, units = "in")

## ====== 1.4 Aging-associated Index (AAI) =============
asi.expr <- gc.expr[final.gene,] %>% t() %>% as.data.frame() %>% bind_cols(gc.subtype[,c('OS.time','OS')])
asi.fit <- coxph(Surv(OS.time, OS) ~ ., data = asi.expr)
gc.subtype$asi.score <- sapply(gc.expr[final.gene,], function(x){
  sum(x * coef(asi.fit)) %>% return()
})
coxph(Surv(OS.time, OS) ~ asi.score, data = gc.subtype) %>% broom::tidy(exponentiate = T)
coxph(Surv(DFI.time, DFI) ~ asi.score, data = gc.subtype) %>% broom::tidy(exponentiate = T)

gc.subtype$ASI.HL <- ifelse(gc.subtype$asi.score > median(gc.subtype$asi.score), 
                            'ASI-H', 'ASI-L')

f2 <- list()
f2[[1]]<-survfit(Surv(OS.time, OS) ~ ASI.HL, data=gc.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='OS (%)',
    legend.title = '',legend.lab=c('AAI-H', 'AAI-L'),
    ggtheme = theme_classic2(), palette = c('#FC4E07','#E7B800'), title = 'TCGA-STAD',
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12
  )

f2[[3]]<-survfit(Surv(DFI.time, DFI) ~ ASI.HL, data=gc.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='DFS (%)',
    legend.title = '',legend.lab=c('AAI-H', 'AAI-L'),
    ggtheme = theme_classic2(),palette = c('#FC4E07','#E7B800'), title = 'TCGA-STAD',
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12
  )


## ====== 1.5 AAI validation =============
load("Datasets/GSE62254.rda")
identical(colnames(GSE62254.expr), GSE62254.subtype$GEO_ID)
GSE62254.subtype$asi.score <- sapply(GSE62254.expr[names(asi.fit$coefficients),], function(x){
  sum(x * coef(asi.fit))
})
coxph(Surv(OS.m, Death) ~ asi.score, data = GSE62254.subtype) %>% broom::tidy(exponentiate = T)
coxph(Surv(DFS.m, Recur) ~ asi.score, data = GSE62254.subtype) %>% broom::tidy(exponentiate = T)

GSE62254.subtype$ASI.HL <- ifelse(GSE62254.subtype$asi.score > median(GSE62254.subtype$asi.score), 'ASI-H','ASI-L')
f2[[2]]<-survfit(Surv(OS.m, Death) ~ ASI.HL, data=GSE62254.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='OS (%)',
    legend.title = '',legend.lab=c('AAI-H', 'AAI-L'),
    ggtheme = theme_classic2(),palette = c('#FC4E07','#E7B800'), title = 'GSE62254',
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12
  )


f2[[4]]<-survfit(Surv(DFS.m, Recur) ~ ASI.HL, data=GSE62254.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='DFS (%)',
    legend.title = '',legend.lab=c('AAI-H', 'AAI-L'),
    ggtheme = theme_classic2(), palette = c('#FC4E07','#E7B800'), title = 'GSE62254',
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12
  )
ggsave(
  "Figure2.pdf",
  arrange_ggsurvplots(f2, print = F,ncol = 2, 
                      nrow = 2, risk.table.height = 0.3),
  width = 10, height = 12
)

# ========= Step 2 =============
## ========= 2.1 Determine the optimal number of clusters =============
best.cluster <- NbClust::NbClust(
  data = t(gc.expr[train.gene,]), method = 'kmeans',
  min.nc = 2, max.nc = 5
)

best_nc <- as.data.frame(t(best.cluster$Best.nc), stringsAsFactors = TRUE)
f3a <- table(best_nc$Number_clusters) %>% as.data.frame() %>% filter(Var1 != 0) %>% 
  ggbarplot(x='Var1', y='Freq', fill='darkgrey', color='darkgrey')+
  labs(x='Number of clusters', y='Frequency among all indices') + 
  ggtitle('Optimal number of clusters')
ggsave(f3a, filename = "Figure3A.pdf", width = 8, height = 6)

gc.cc.fit <- ConsensusClusterPlus::ConsensusClusterPlus(
  d = as.matrix(gc.expr[train.gene,]), maxK = 5, clusterAlg = 'km',
  plot='pdf'
)

pdf("Figure3B.pdf", width = 10, height = 8)
factoextra::fviz_silhouette(
  sil.obj = cluster::silhouette(
    gc.cc.fit[[2]]$consensusClass,
    dist(t(gc.expr[train.gene,]))
  ), palette = c('#609EA2','#C27664')
)
dev.off()

gc.subtype$Consensus.subtype <- paste0('Cluster', gc.cc.fit[[2]]$consensusClass)

f3cd <- list()
f3cd[[1]]<-survfit(Surv(OS.time, OS) ~ Consensus.subtype, data=gc.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='OS (%)', palette = c('#609EA2','#C27664'),
    legend.title = '',legend.lab=paste0('Cluster',1:2),
    ggtheme = theme_classic2(), title = 'TCGA-STAD',
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12
  )

f3cd[[3]]<-survfit(Surv(DFI.time, DFI) ~ Consensus.subtype, data=gc.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='DFS (%)',
    legend.title = '',legend.lab=paste0('Cluster',1:2),title = 'TCGA-STAD',
    ggtheme = theme_classic2(), palette = c('#609EA2','#C27664'),
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12
  )

## ========= 2. Build a molecular subtyping model =============
library(randomForest)
gc.consensus.cluster.rf <- randomForest::randomForest(
  x = gc.expr[names(sort(apply(gc.expr, 1, mad), decreasing = T)[1:1000]),] %>% t() %>% scale(), 
  y = as.factor(gc.subtype$Consensus.subtype)
)

## ======== 3. Validate the model performace =============
GSE62254.subtype$Consensus.subtype <- predict(
  gc.consensus.cluster.rf,
  GSE62254.expr[rownames(importance(gc.consensus.cluster.rf)),] %>% t() %>% scale()
) %>% as.character()

f3cd[[2]]<-survfit(Surv(OS.m, Death) ~ Consensus.subtype, data=GSE62254.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='OS (%)',
    legend.title = '',legend.lab=paste0('Cluster',1:2),
    ggtheme = theme_classic2(), title = 'GSE62254',
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12,palette = c('#609EA2','#C27664')
  )

f3cd[[4]]<-survfit(Surv(DFS.m, Recur) ~ Consensus.subtype, data=GSE62254.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='DFS (%)',
    legend.title = '',legend.lab=paste0('Cluster',1:2),
    ggtheme = theme_classic2(),title = 'GSE62254',
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12,palette = c('#609EA2','#C27664')
  )
ggsave(
  "Figure3C-D.pdf",
  arrange_ggsurvplots(f3cd, print = F,ncol = 2, 
                      nrow = 2, risk.table.height = 0.3),
  width = 10, height = 12
)


# ========= Step 3 =============
## ======== 3.1 Perform the differential analysis ===========
identical(colnames(gc.expr), gc.subtype$sample)
consen.subtype.limma <- limma::lmFit(object = gc.expr, design = model.matrix(~ifelse(gc.subtype$Consensus.subtype == 'Cluster2',1,0))) %>% 
  limma::eBayes() %>% limma::topTable(coef = 2, number = Inf, adjust.method = 'BH')

## ======== 3.2 Prioritize the potential regulons and target genes ===========
load('Datasets/TF_EMT_signature.rda')
mRNA.limma <- consen.subtype.limma[setdiff(rownames(gc.expr), tf.emt.sgt$TF$Symbol),]
tf.limma <- consen.subtype.limma[tf.emt.sgt$TF$Symbol,]
data4enrich <- mRNA.limma$logFC
names(data4enrich) <- rownames(mRNA.limma)

TF.regulon <- rownames(tf.limma %>% filter(adj.P.Val < 0.001 & abs(logFC) > 0.5))
mRNA.target <- rownames(mRNA.limma %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1))

## ======== 3.3 Establish the integrative network inference ===========
rtni <- RTN::tni.constructor(expData = gc.expr[TF.regulon,] %>% t() %>% scale() %>% t() %>% as.data.frame() %>% 
                               bind_rows(
                                 gc.expr[mRNA.target,] %>% t() %>% scale() %>% t() %>% as.data.frame()
                               ) %>% as.matrix(), regulatoryElements = TF.regulon) %>% 
  RTN::tni.permutation(nPermutations = 1000) %>% RTN::tni.bootstrap() %>% RTN::tni.dpi.filter()

## ======== 3.4 Perform the master regulatory analysis ===========
rtna <- RTN::tni2tna.preprocess(object = rtni, phenotype = data4enrich, hits = tf.emt.sgt$EMT.gene$GeneSymbol) %>% 
  RTN::tna.mra() %>% RTN::tna.get(what="mra", ntop = -1)

## ======== 3.5 Prioritize the master regulatory TFs ===========
regulon <- rtna$Regulon[rtna$Pvalue < 0.05]
tmp <- gc.expr[regulon,] %>% t() %>% as.data.frame() %>% bind_cols(gc.subtype)
regulon <- lapply(regulon, function(x){
  res.cox<-base::eval(parse(text = sprintf('coxph(Surv(OS.time, OS) ~ %s, data=tmp)',x)))
  return(res.cox %>% broom::tidy())
}) %>% do.call(rbind, .) %>% filter(p.value < 0.05) %>% pull(term)

# ========= Step 4 =============
## ======== 4.1 Perform the drug sensitivity analysis ===========
load('Datasets/CTRP_GC.rda')
genes<-c('ELK3','SOX7')
drug.gene.df<-expand.grid(genes, rownames(ctrp.drug),stringsAsFactors = F) %>% as.data.frame()
colnames(drug.gene.df)<-c('gene','drug')
drug.gene.df$anova.p<-NA
drug.gene.df$auc<-NA

for(i in 1:nrow(drug.gene.df)){
  df<-data.frame(
    gene=ctrp.expr[as.character(drug.gene.df$gene[i]),] %>% unlist(),
    drug=ctrp.drug[as.character(drug.gene.df$drug[i]),] %>% unlist()
  )
  drug.gene.df$anova.p[i]<-anova(lm(drug~gene, df))[1,5]
  drug.gene.df$auc[i]<-pROC::roc(
    ifelse(df$drug > median(df$drug),1,0),
    df$gene,quiet=T,direction='<'
  )$auc
}
rm(i,df)


f4a<-drug.gene.df[drug.gene.df$gene == 'ELK3',] %>% left_join(drug.info[,c('ID','CLINICAL.STATUS')],by=c('drug'='ID')) %>% 
  ggplot(aes(x=auc,y=-log10(anova.p)))+
  geom_point(aes(color=I(ifelse((auc >= 0.75 | auc <= 0.25) & anova.p < 0.05 & CLINICAL.STATUS =='FDA approved', 'red','lightgrey'))))+
  theme_classic()+labs(x='AUC',y='ANOVA -log10 (P)',title='ELK3')+
  geom_hline(yintercept = -log10(0.05),lty='dashed',color='lightgrey')+
  geom_vline(xintercept = c(0.25,0.75),lty='dashed',color='lightgrey')+
  ggrepel::geom_label_repel(aes(label=ifelse((auc >= 0.75 | auc <= 0.25) & anova.p < 0.05 & CLINICAL.STATUS =='FDA approved', drug,'')),label.size = NA,max.overlaps=Inf)+
  theme(
    legend.position = 'none',
    axis.text= element_text(size=15,family ="sans"),
    axis.title= element_text(size = 15,family ="sans")
  )
ggsave('Figure4A.pdf', f4a, width = 6, height = 6, units = 'in')

f4b<-drug.gene.df[drug.gene.df$gene == 'SOX7',] %>% left_join(drug.info[,c('ID','CLINICAL.STATUS')],by=c('drug'='ID')) %>% 
  ggplot(aes(x=auc,y=-log10(anova.p)))+
  geom_point(aes(color=I(ifelse((auc >= 0.75 | auc <= 0.25) & anova.p < 0.05, 'red','lightgrey'))))+
  theme_classic()+labs(x='AUC',y='ANOVA -log10 (P)',title='SOX7')+
  geom_hline(yintercept = -log10(0.05),lty='dashed',color='lightgrey')+
  geom_vline(xintercept = c(0.25,0.75),lty='dashed',color='lightgrey')+
  ggrepel::geom_label_repel(aes(label=ifelse((auc >= 0.75 | auc <= 0.25) & anova.p < 0.05, drug,'')),label.size = NA,max.overlaps=Inf)+
  theme(
    legend.position = 'none',
    axis.text= element_text(size=15,family ="sans"),
    axis.title= element_text(size = 15,family ="sans")
  )
ggsave('Figure4B.pdf', f4b, width = 6, height = 6, units = 'in')
