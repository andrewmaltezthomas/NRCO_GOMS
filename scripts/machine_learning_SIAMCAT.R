library(SIAMCAT)
library(phyloseq)
metaphlan_species <- readRDS("data/R_vs_NR_meta_analysis.rds")
metaphlan_species <- subset_samples(metaphlan_species, !(metaphlan_species@sam_data$Cohort %in% c("Leeds_LeeK_2022", "Barcelona_LeeK_2022", "WindTT_2020")))
metaphlan_species@sam_data$Cohort <- gsub("RoutyB_2022_Lung", "RoutyB_2018_Lung", metaphlan_species@sam_data$Cohort)
metaphlan_species <- prune_taxa(taxa_sums(metaphlan_species)>0, metaphlan_species)
print("Running Leave-one-centre-out...")
datasets <- as.character(unique(sample_data(metaphlan_species)$Cohort))
labels_to_predict <- "ORR"
siamcat_ml_results <- expand.grid(datasets, labels_to_predict)
siamcat_ml_results <- as.data.frame(siamcat_ml_results)
siamcat_ml_results$AUC_ROC <- NA
siamcat_ml_results$AUC_PRC <- NA
colnames(siamcat_ml_results) <- c("Cohort", "Label", "AUC_ROC")
for (j in 1:length(datasets)) {
  all_datasets <- subset_samples(metaphlan_species, sample_data(metaphlan_species)$Cohort != datasets[j])
  leftout_dataset <- subset_samples(metaphlan_species, sample_data(metaphlan_species)$Cohort == datasets[j])
  if ( any(table(sample_data(all_datasets)[,labels_to_predict]) <= 5)  | any(table(sample_data(leftout_dataset)[,labels_to_predict]) <= 5) | dim(table(sample_data(leftout_dataset)[,labels_to_predict])) 
       < 2 | dim(table(sample_data(all_datasets)[,labels_to_predict])) < 2 ) {
    next()
  }
  else {
    labels.response <- create.label(meta=sample_data(all_datasets), label=labels_to_predict, case="R")
    siamcat <- siamcat(feat=otu_table(all_datasets) / 100, label=labels.response, meta=sample_data(all_datasets))
    siamcat <- normalize.features(siamcat, norm.method = "log.clr", norm.param = list(log.n0 = 1e-06, sd.min.q = 0.1), verbose = 2, feature.type = "original")
    siamcat <-  create.data.split(siamcat, num.folds = 2, num.resample = 1)
    siamcat <- train.model(siamcat, method = "lasso")
    siamcat <- make.predictions(siamcat)
    siamcat <-  evaluate.predictions(siamcat)
    
    #LeftOut Dataset
    labels.response_leftout <- create.label(meta=sample_data(leftout_dataset), label=labels_to_predict, case="R")
    siamcat_leftout <- siamcat(feat=otu_table(leftout_dataset) / 100, label=labels.response_leftout, meta=sample_data(leftout_dataset))
    siamcat_leftout <- normalize.features(siamcat_leftout, verbose = 2, norm.param = norm_params(siamcat), feature.type = "original")
    siamcat_leftout <- make.predictions(siamcat, siamcat.holdout = siamcat_leftout, normalize.holdout = FALSE)
    siamcat_leftout <-  evaluate.predictions(siamcat_leftout)
    position <- which(siamcat_ml_results$Cohort == datasets[j] & siamcat_ml_results$Label == labels_to_predict)
    siamcat_ml_results[position,3] <- eval_data(siamcat_leftout)$auroc
  }
  #saveRDS(siamcat_leftout, paste("LODO_",  datasets[j], "_", labels_to_predict, ".rds", sep = ""))
}
write_table(siamcat_ml_results, "results/SIAMCAT_lodo.txt", row.names = F, col.names = F, sep = "\t", quote = F)
