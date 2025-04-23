# Protocol for deciphering tumor prognosis and heterogeneity by transcriptome using machine learning

## Summary
Identifying high-risk patients and dissecting tumor heterogeneity are critical for the advancement of individualized cancer therapeutics. Thus, leveraging specific biological or clinical parameters, such as age, is essential for developing clinical indicators that guide therapy and stratify patients into distinct subgroups. This approach further enables the prioritization of key druggable targets for those subgroups exhibiting the most aggressive biological behavior.

## Major work
In this study, we utilized aging as a biological criterion to identify high-risk gastric tumors and classify them into clinically relevant subtypes. We gathered 445 aging-associated genes from the Aging Atlas database and constructed an independent prognostic indicator model, termed the Aging-Associated Index (AAI), using robustly prioritized signature genes through a machine learning approach. Subsequently, molecular stratification was performed using an unsupervised clustering method, and a machine learning-based classifier was developed to predict patient subtypes in an independent cohort. Additionally, we inferred a subtype-specific regulatory network of transcription factors (TFs) and prioritized key TFs using master regulatory analysis. Finally, a drug sensitivity analysis was conducted to repurpose drugs targeting the master regulatory TFs. The pipeline described here is applicable to any cancer type and based on any biological process.

## Codes used in the protocol
* step1.R: Develop a prognostic model using robustly prioritized signature gene.
* step2.R: Identify molecular subtypes via unsupervised clustering and train a subtype classifier.
* step3.R: Establish a subtype-specific regulatory network and perform master regulatory analysis.
* step4.R: Prioritize therapeutic candidates with drug sensitivity analysis.

## Data used in the protocol
* Datasets/TCGA-STAD.rda: Processed TCGA-STAD gene expression profiles, clinical information, and aging-associated genes.
* Datasets/GSE62254.rda: Processed gene expression profiles and clinical information of validation dataset GSE62254.
* Datasets/TF_EMT_signature.rda: An R list contains human transcription factors and epithelial-to-mesenchymal transition (EMT) siganture genes.
* Datasets/CTRP_GC.rda: Processing gene expression profiles, drug activity score profiles of human cell lines, and drug information.


