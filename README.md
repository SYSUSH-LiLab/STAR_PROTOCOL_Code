# Protocol for deciphering tumor prognosis and heterogeneity by transcriptome using machine learning

## Summary
Identifying high-risk patients and dissecting tumor heterogeneity are critical for the advancement of individualized cancer therapeutics. Thus, leveraging specific biological or clinical parameters, such as age, is essential for developing clinical indicators that guide therapy and stratify patients into distinct subgroups. This approach further enables the prioritization of key druggable targets for those subgroups exhibiting the most aggressive biological behavior.

## Major work
In this study, we utilized aging as a biological criterion to identify high-risk gastric tumors and classify them into clinically relevant subtypes. We gathered 445 aging-associated genes from the Aging Atlas database and constructed an independent prognostic indicator model, termed the Aging-Associated Index (AAI), using robustly prioritized signature genes through a machine learning approach. Subsequently, molecular stratification was performed using an unsupervised clustering method, and a machine learning-based classifier was developed to predict patient subtypes in an independent cohort. Additionally, we inferred a subtype-specific regulatory network of transcription factors (TFs) and prioritized key TFs using master regulatory analysis. Finally, a drug sensitivity analysis was conducted to repurpose drugs targeting the master regulatory TFs. The pipeline described here is applicable to any cancer type and based on any biological process.


