## ========= 1. Determine the optimal number of clusters =============
### Use the NbClust package to determine the optimal number of clusters based on multiple indices.
### Input: gene expression matrix (transposed), limited to selected prognostic genes.
### K-means clustering is applied, testing cluster numbers from 2 to 5.
best.cluster <- NbClust::NbClust(
  data = t(gc.expr[train.gene,]), method = 'kmeans',
  min.nc = 2, max.nc = 5
)

### Extract the best number of clusters selected by each index
best_nc <- as.data.frame(t(best.cluster$Best.nc), stringsAsFactors = TRUE)
### Summarize and plot the frequency of cluster numbers chosen across all indices
table(best_nc$Number_clusters) %>% as.data.frame() %>% filter(Var1 != 0) %>% 
  ggbarplot(x='Var1', y='Freq', fill='darkgrey', color='darkgrey')+
  labs(x='Number of clusters', y='Frequency among all indices') + ggtitle('Optimal number of clusters')

### Perform consensus clustering using k-means (km) on the training dataset
### The maximum number of clusters considered is 5
### Consensus matrices and tracking plots will be saved as PDF files in the working directory
gc.cc.fit <- ConsensusClusterPlus::ConsensusClusterPlus(
  d = as.matrix(gc.expr[train.gene,]), maxK = 5, clusterAlg = 'km',
  plot='pdf'
)

### Compute and visualize the silhouette width for each sample to assess clustering quality
factoextra::fviz_silhouette(
  sil.obj = cluster::silhouette(
    gc.cc.fit[[2]]$consensusClass,
    dist(t(gc.expr[train.gene,]))
  ), palette = c('#609EA2','#C27664')
)

### Only the results for 2 clusters (gc.cc.fit[[2]]) are visualized here
### Assign patients to one of the two consensus clusters and label as 'Cluster1' or 'Cluster2'
gc.subtype$Consensus.subtype <- paste0('Cluster', gc.cc.fit[[2]]$consensusClass)

### Generate a Kaplan-Meier plot for overall survival (OS) stratified by consensus subtype in TCGA-STAD.
survfit(Surv(OS.time, OS) ~ Consensus.subtype, data=gc.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='OS (%)', palette = c('#609EA2','#C27664'),
    legend.title = '',legend.lab=paste0('Cluster',1:2),
    ggtheme = theme_classic2(), title = 'TCGA-STAD',
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12
  )

### Generate a Kaplan-Meier plot for disease-free survival (DFS) by consensus subtype.
survfit(Surv(DFI.time, DFI) ~ Consensus.subtype, data=gc.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='DFS (%)',
    legend.title = '',legend.lab=paste0('Cluster',1:2),title = 'TCGA-STAD',
    ggtheme = theme_classic2(), palette = c('#609EA2','#C27664'),
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12
  )

## ========= 2. Build a molecular subtyping model =============
### Train a random forest classifier using the top 1000 most variable genes (based on MAD)
### Inputs are z-score scaled expression values
### Output variable is the consensus subtype label (Cluster1 or Cluster2)
gc.consensus.cluster.rf <- randomForest::randomForest(
  x = gc.expr[names(sort(apply(gc.expr, 1, mad), decreasing = T)[1:1000]),] %>% t() %>% scale(), 
  y = as.factor(gc.subtype$Consensus.subtype)
)

## ======== 3. Validate the model performance =============
### Predict consensus subtypes in the validation dataset (GSE62254) using the trained random forest classifier
### Only genes used in the classifier are retained, and input expression is standardized
GSE62254.subtype$Consensus.subtype <- predict(
  gc.consensus.cluster.rf,
  GSE62254.expr[rownames(importance(gc.consensus.cluster.rf)),] %>% t() %>% scale()
) %>% as.character()

### Perform the survival analysis for GSE62254 consensus subtypes in OS and DFS ====
survfit(Surv(OS.m, Death) ~ Consensus.subtype, data=GSE62254.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='OS (%)',
    legend.title = '',legend.lab=paste0('Cluster',1:2),
    ggtheme = theme_classic2(), title = 'GSE62254',
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12,palette = c('#609EA2','#C27664')
  )

survfit(Surv(DFS.m, Recur) ~ Consensus.subtype, data=GSE62254.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='DFS (%)',
    legend.title = '',legend.lab=paste0('Cluster',1:2),
    ggtheme = theme_classic2(),title = 'GSE62254',
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12,palette = c('#609EA2','#C27664')
  )
