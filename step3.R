## ======== 1. Perform the differential analysis ===========
identical(colnames(gc.expr), gc.subtype$sample)
consen.subtype.limma <- limma::lmFit(object = gc.expr, design = model.matrix(~ifelse(gc.subtype$Consensus.subtype == 'Cluster2',1,0))) %>% 
  limma::eBayes() %>% limma::topTable(coef = 2, number = Inf, adjust.method = 'BH')

## ======== 2. Prioritize the potential regulons and target genes ===========
load('Datasets/TF_EMT_signature.rda')
mRNA.limma <- consen.subtype.limma[setdiff(rownames(gc.expr), tf.emt.sgt$TF$Symbol),]
tf.limma <- consen.subtype.limma[tf.emt.sgt$TF$Symbol,]
data4enrich <- mRNA.limma$logFC
names(data4enrich) <- rownames(mRNA.limma)

TF.regulon <- rownames(tf.limma %>% filter(adj.P.Val < 0.001 & abs(logFC) > 0.5))
mRNA.target <- rownames(mRNA.limma %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1))

## ======== 3. Establish the integrative network inference ===========
rtni <- RTN::tni.constructor(expData = gc.expr[TF.regulon,] %>% t() %>% scale() %>% t() %>% as.data.frame() %>% 
                               bind_rows(
                                 gc.expr[mRNA.target,] %>% t() %>% scale() %>% t() %>% as.data.frame()
                               ) %>% as.matrix(), regulatoryElements = TF.regulon) %>% 
  RTN::tni.permutation(nPermutations = 1000) %>% RTN::tni.bootstrap() %>% RTN::tni.dpi.filter()

## ======== 4. Perform the master regulatory analysis ===========
rtna <- RTN::tni2tna.preprocess(object = rtni, phenotype = data4enrich, hits = tf.emt.sgt$EMT.gene$GeneSymbol) %>% 
  RTN::tna.mra() %>% RTN::tna.get(what="mra", ntop = -1)

## ======== 5. Prioritize the master regulatory TFs ===========
regulon <- rtna$Regulon[rtna$Pvalue < 0.05]
tmp <- gc.expr[regulon,] %>% t() %>% as.data.frame() %>% bind_cols(gc.subtype)
regulon <- lapply(regulon, function(x){
  res.cox<-base::eval(parse(text = sprintf('coxph(Surv(OS.time, OS) ~ %s, data=tmp)',x)))
  return(res.cox %>% broom::tidy())
}) %>% do.call(rbind, .) %>% filter(p.value < 0.05) %>% pull(term)
