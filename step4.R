## 1. ================Data pre-processing of CTRP dataset==================
genes<-c('ELK3','SOX7')
ctrp.drug<-data.table::fread('raw_data/data_CTRP-Broad-MIT_act.txt',data.table = F)
drug.info<-ctrp.drug[,1:4]
ctrp.drug<-ctrp.drug[,5:ncol(ctrp.drug)]
rownames(ctrp.drug)<-drug.info$ID
ctrp.expr<-data.table::fread('raw_data/data_CTRP-Broad-MIT_exp.txt',data.table = F)
rownames(ctrp.expr)<-ctrp.expr$ID
ctrp.expr$ID<-NULL

ctrp.sample.info<-data.table::fread('raw_data/Cell_line_annotation_ctrp.txt',data.table = F)
ctrp.sample.info<-ctrp.sample.info[str_detect(ctrp.sample.info$TissueType, 'stomach'),]
cell.line<-intersect(ctrp.sample.info$Name, colnames(ctrp.expr)) %>% intersect(colnames(ctrp.drug))
ctrp.drug<-ctrp.drug[,cell.line]
ctrp.expr<-ctrp.expr[,cell.line]
ctrp.drug<-ctrp.drug[apply(ctrp.drug,1,function(x){sum(is.na(x)) < (length(cell.line)*0.5)}),]
ctrp.drug<-impute::impute.knn(ctrp.drug %>% as.matrix())[['data']] %>% as.data.frame()
ctrp.drug<- (1-ctrp.drug/30)

## ======== 2. Perform the drug sensitivity analysis ===========
load('Datasets/CTRP_GC.rda')
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

lapply(genes, function(x){
  df<-drug.gene.df[drug.gene.df$gene == x,] %>% left_join(drug.info[,c('ID','CLINICAL.STATUS')],by=c('drug'='ID'))
  ggplot(df,aes(x=auc,y=-log10(anova.p)))+geom_point(aes(color=I(ifelse((auc >= 0.75 | auc <= 0.25) & anova.p < 0.05, 'red','lightgrey'))))+
    theme_classic()+labs(x='AUC',y='ANOVA -log10 (P)',title=x)+
    geom_hline(yintercept = -log10(0.05),lty='dashed',color='lightgrey')+
    geom_vline(xintercept = c(0.25,0.75),lty='dashed',color='lightgrey')+
    ggrepel::geom_label_repel(aes(label=ifelse((auc >= 0.75 | auc <= 0.25) & anova.p < 0.05, drug,'')),label.size = NA,max.overlaps=Inf)+
    theme(
      legend.position = 'none',
      axis.text= element_text(size=15,family ="sans"),
      axis.title= element_text(size = 15,family ="sans")
    )
}) %>% cowplot::plot_grid(plotlist = ., ncol=2)
