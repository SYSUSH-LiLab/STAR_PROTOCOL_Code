#Since the data is too large to upload, we only provide the code and processed data.
## ======== 1. Perform the drug sensitivity analysis ===========
### Load the processed CTRP drug sensitivity and gene expression dataset
load('Datasets/CTRP_GC.rda')

### Define the genes of interest for drug association analysis
genes<-c('ELK3','SOX7')

### Create a data frame with all possible combinations of gene-drug pairs
drug.gene.df<-expand.grid(genes, rownames(ctrp.drug),stringsAsFactors = F) %>% as.data.frame()
colnames(drug.gene.df)<-c('gene','drug')
drug.gene.df$anova.p<-NA
drug.gene.df$auc<-NA

### Calculate the AUC and ANOVA p-value for each gene-drug pair ====
for(i in 1:nrow(drug.gene.df)){
  ### Extract the expression level of the current gene and the drug sensitivity value
  df<-data.frame(
    gene=ctrp.expr[as.character(drug.gene.df$gene[i]),] %>% unlist(),
    drug=ctrp.drug[as.character(drug.gene.df$drug[i]),] %>% unlist()
  )
  
  ### Perform ANOVA test to assess linear association between gene expression and drug response
  drug.gene.df$anova.p[i]<-anova(lm(drug~gene, df))[1,5]
  
  ### Compute AUC from ROC curve to evaluate classification power of gene expression
  drug.gene.df$auc[i]<-pROC::roc(
    ifelse(df$drug > median(df$drug),1,0), ### define binary drug response
    df$gene,quiet=T,direction='<' ### test whether higher gene expression predicts sensitivity
  )$auc
}
rm(i,df)

### ===== Visualize the results for each gene-drug pair ====
drug.gene.df[drug.gene.df$gene == 'ELK3',] %>% left_join(drug.info[,c('ID','CLINICAL.STATUS')],by=c('drug'='ID')) %>% 
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

drug.gene.df[drug.gene.df$gene == 'SOX7',] %>% left_join(drug.info[,c('ID','CLINICAL.STATUS')],by=c('drug'='ID')) %>% 
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
