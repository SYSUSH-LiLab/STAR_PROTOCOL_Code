## ======== 1. Perform the differential analysis using limma ===========
### Confirm that the sample names in the gene expression matrix match the sample metadata
identical(colnames(gc.expr), gc.subtype$sample)

### Use the limma pipeline to perform differential expression analysis between Cluster2 vs. others
### eBayes moderates the variance, and topTable extracts the DE results for the contrast (coef=2)
### All genes are returned (number = Inf), with multiple testing correction (Benjamini-Hochberg)
consen.subtype.limma <- limma::lmFit(object = gc.expr, design = model.matrix(~ifelse(gc.subtype$Consensus.subtype == 'Cluster2',1,0))) %>% 
  limma::eBayes() %>% limma::topTable(coef = 2, number = Inf, adjust.method = 'BH')

## ======== 2. Prioritize the potential regulons and target genes ===========
### Load a pre-defined list of transcription factors (TFs) and EMT-related signature genes
load('Datasets/TF_EMT_signature.rda')
### Extract differentially expressed genes excluding transcription factors — assumed to be targets
mRNA.limma <- consen.subtype.limma[setdiff(rownames(gc.expr), tf.emt.sgt$TF$Symbol),]
### Extract differentially expressed transcription factors from the limma result
tf.limma <- consen.subtype.limma[tf.emt.sgt$TF$Symbol,]

### Create a named vector of log fold changes for target genes, 
### used later in enrichment or regulatory analysis
data4enrich <- mRNA.limma$logFC
names(data4enrich) <- rownames(mRNA.limma)

### === Prioritize the regulons and target genes based on the differential expression analysis ====
TF.regulon <- rownames(tf.limma %>% filter(adj.P.Val < 0.001 & abs(logFC) > 0.5))
mRNA.target <- rownames(mRNA.limma %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1))

## ======== 3. Establish the integrative network inference ===========
### Build the transcriptional network using RTN:
### - Input matrix includes scaled expression of candidate TFs and target genes
### - TFs are defined as regulators
### - RTN performs 1000 permutations to assess statistical significance,
###   then bootstrapping and DPI (data processing inequality) filtering to infer direct TF-target interactions
rtni <- RTN::tni.constructor(expData = gc.expr[TF.regulon,] %>% t() %>% scale() %>% t() %>% as.data.frame() %>% 
                               bind_rows(
                                 gc.expr[mRNA.target,] %>% t() %>% scale() %>% t() %>% as.data.frame()
                               ) %>% as.matrix(), regulatoryElements = TF.regulon) %>% 
  RTN::tni.permutation(nPermutations = 1000) %>% RTN::tni.bootstrap() %>% RTN::tni.dpi.filter()

## ======== 4. Perform the master regulatory analysis ===========
### Perform Master Regulator Analysis (MRA) on the inferred network:
### - Use log fold changes of target genes as phenotype signature
### - Use EMT-related genes as the hit set
### - Identify TFs whose regulons are significantly enriched in phenotype signature (MRA)
### - Extract all results (ntop = -1 means no filtering on top N)
rtna <- RTN::tni2tna.preprocess(object = rtni, phenotype = data4enrich, hits = tf.emt.sgt$EMT.gene$GeneSymbol) %>% 
  RTN::tna.mra() %>% RTN::tna.get(what="mra", ntop = -1)

## ======== 5. Prioritize the master regulatory TFs ===========
### Select TFs with statistically significant enrichment (P < 0.05) from MRA output
regulon <- rtna$Regulon[rtna$Pvalue < 0.05]
tmp <- gc.expr[regulon,] %>% t() %>% as.data.frame() %>% bind_cols(gc.subtype)

### For each candidate master regulator TF:
### - Perform univariate Cox regression on its expression and OS.
### - Collect TFs with P < 0.05 (significant prognostic value).
regulon <- lapply(regulon, function(x){
  res.cox<-base::eval(parse(text = sprintf('coxph(Surv(OS.time, OS) ~ %s, data=tmp)',x)))
  return(res.cox %>% broom::tidy())
}) %>% do.call(rbind, .) %>% filter(p.value < 0.05) %>% pull(term)
