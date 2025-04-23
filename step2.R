## ========= 1. Determine the optimal number of clusters =============
best.cluster <- NbClust::NbClust(
  data = t(gc.expr[train.gene,]), method = 'kmeans',
  min.nc = 2, max.nc = 5
)

gc.cc.fit <- ConsensusClusterPlus::ConsensusClusterPlus(
  d = as.matrix(gc.expr[train.gene,]), maxK = 5, clusterAlg = 'km',
  plot='pdf'
)

factoextra::fviz_silhouette(
  sil.obj = cluster::silhouette(
    gc.cc.fit[[2]]$consensusClass,
    dist(t(gc.expr[train.gene,]))
  ), palette = c('#609EA2','#C27664')
)

gc.subtype$Consensus.subtype <- paste0('Cluster', gc.cc.fit[[2]]$consensusClass)
survfit(Surv(OS.time, OS) ~ Consensus.subtype, data=gc.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='OS (%)', palette = c('#609EA2','#C27664'),
    legend.title = '',legend.lab=paste0('Cluster',1:2),
    ggtheme = theme_classic2(),
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12
  )

survfit(Surv(DFI.time, DFI) ~ Consensus.subtype, data=gc.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='DFS (%)',
    legend.title = '',legend.lab=paste0('Cluster',1:2),
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

survfit(Surv(OS.m, Death) ~ Consensus.subtype, data=GSE62254.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='OS (%)',
    legend.title = '',legend.lab=paste0('Cluster',1:2),
    ggtheme = theme_classic2(),
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12,palette = c('#609EA2','#C27664')
  )

survfit(Surv(DFS.m, Recur) ~ Consensus.subtype, data=GSE62254.subtype) %>%
  ggsurvplot(
    pval = TRUE,conf.int = F,
    risk.table = T,risk.table.col = "strata",size=1,pval.size=8,
    xlab='Follow up (months)', ylab='DFS (%)',
    legend.title = '',legend.lab=paste0('Cluster',1:2),
    ggtheme = theme_classic2(),
    font.x = 15,font.y=15,font.main=18,font.legend=15,font.tickslab=12,palette = c('#609EA2','#C27664')
  )
