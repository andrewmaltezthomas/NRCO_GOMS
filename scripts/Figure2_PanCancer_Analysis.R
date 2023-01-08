######### Pairwise Pibble Analysis #######################################
set.seed(1986)
write_pibble <- function(pibble_fit=pibble_fit, focus_var=focus_var, siglevel=siglevel) {
  p <- as.data.frame.table(pibble_fit$Lambda) %>% 
    filter(Var2==focus_var) %>% 
    mutate(idx=as.numeric(as.factor(Var1))) %>% 
    select(idx, Freq, x=Var2, y=Var1) %>% 
    group_by(idx) %>%
    ggdist::median_qi(Freq, .width=c(0.5, 0.75, 0.90, 0.95, 0.97)) %>%
    mutate(species=rep(rownames(pibble_fit$Y),5)) %>% 
    filter(.width == 0.95) 
  
  return(p)    
}
library(reshape2)
library(dplyr)
library(ggplot2)
library(tidyr)
library(phyloseq)
library(fido)
library(stringr)
library(metagMisc)
library(zCompositions)
library(driver)
library(ggdist)
library(purrr)
cancers <- c("Breast", "Melanoma", "Ovary", "Pancreatic", "Prostate", "Colorectal", "Kidney", "Lung")
for (cancer in cancers) {
  dir.create(paste("results/", cancer, "_PIBBLE", sep = ""))
  metaphlan_species <- readRDS("data/Cancer_vs_Controls_meta_analysis.rds")
  datasets_control <- as.character(unique(sample_data(metaphlan_species)$Cohort[metaphlan_species@sam_data$Cancer == "Control"]))
  datasets_cancer <- as.character(unique(sample_data(metaphlan_species)$Cohort[metaphlan_species@sam_data$Cancer == cancer]))
  results <- expand.grid(datasets_cancer, datasets_control)
  colnames(results) <- c("Cancer_Dataset", "Control_Dataset")
  
  if(length(which(as.character(results$Cancer_Dataset) == as.character(results$Control_Dataset))) != 0) {
    results <- results[-which(as.character(results$Cancer_Dataset) == as.character(results$Control_Dataset)),]
  }
  
  for (i in 1:nrow(results)) {
    
    dataset <- as.character(results$Control_Dataset[i])
    j <- as.character(results$Cancer_Dataset[i])
    print(paste(dataset, "_v_", j, sep = ""))
    
    subset_dataset1 <- subset_samples(metaphlan_species, rownames(metaphlan_species@sam_data) %in% rownames(metaphlan_species@sam_data)[which(sample_data(metaphlan_species)$Cohort == dataset)])
    subset_dataset2 <- subset_samples(metaphlan_species, rownames(metaphlan_species@sam_data) %in% rownames(metaphlan_species@sam_data)[which(sample_data(metaphlan_species)$Cohort == j)])
    
    subset_dataset1 <- metagMisc::phyloseq_filter_prevalence(subset_dataset1, prev.trh=0.05)
    subset_dataset1 <- subset_taxa(subset_dataset1, rowSums(otu_table(subset_dataset1)) != 0)
    
    subset_dataset2 <- metagMisc::phyloseq_filter_prevalence(subset_dataset2, prev.trh=0.05)
    subset_dataset2 <- subset_taxa(subset_dataset2, rowSums(otu_table(subset_dataset2)) != 0)
    
    subset_dataset <- merge_phyloseq(subset_dataset1, subset_dataset2)
    subset_dataset@sam_data$age <- as.numeric(subset_dataset@sam_data$age)
    subset_dataset@sam_data$age <- scale(subset_dataset@sam_data$age)[,1]
    subset_dataset@sam_data$SampleID <- rownames(subset_dataset@sam_data)
    
    if (dim(table(subset_dataset@sam_data$gender)) == 2 ) {
      mdat_fido <- as(sample_data(subset_dataset), "data.frame") %>% mutate(SampleID = as.character(SampleID), gender = relevel(factor(gender, ordered = FALSE), ref="male"), age = age, Cohort = relevel(factor(Cohort, ordered = FALSE), ref = dataset)) %>% select(SampleID, Cohort, gender, age)  
      ps_fido <- phyloseq(otu_table(subset_dataset, taxa_are_rows=T), sample_data(mdat_fido))
      ps_fido <- prune_taxa(taxa_sums(ps_fido) > 0, ps_fido)
      ps_fido <- prune_samples(sample_sums(ps_fido) > 0, ps_fido)
      if (j == "NagataN_2022") { # Nagata has age as categorical
        f <- reformulate(termlabels=c("Cohort", "gender"))
        X <- t(model.matrix(f, data=as(sample_data(ps_fido), "data.frame")))  
      }
      else {
        f <- reformulate(termlabels=c("Cohort", "age", "gender"))
        X <- t(model.matrix(f, data=as(sample_data(ps_fido), "data.frame")))
      } 
    }
    
    else if (dim(table(subset_dataset@sam_data$gender)) == 1 ) {
      mdat_fido <- as(sample_data(subset_dataset), "data.frame") %>% mutate(SampleID = as.character(SampleID), gender = gender, age = age, Cohort = relevel(factor(Cohort, ordered = FALSE), ref = dataset)) %>% select(SampleID, Cohort, gender, age)  
      ps_fido <- phyloseq(otu_table(subset_dataset, taxa_are_rows=T), sample_data(mdat_fido))
      ps_fido <- prune_taxa(taxa_sums(ps_fido) > 0, ps_fido)
      ps_fido <- prune_samples(sample_sums(ps_fido) > 0, ps_fido)
      f <- reformulate(termlabels=c("Cohort", "age"))
      X <- t(model.matrix(f, data=as(sample_data(ps_fido), "data.frame")))
    }
    
    ## Run the pibble model
    Y <- driver::miniclo(t(otu_table(ps_fido))) 
    Y <- as.matrix(zCompositions::cmultRepl(Y)) 
    Y <- t(Y) 
    N <- ncol(Y)
    D <- nrow(Y)
    upsilon <- D+3
    Omega <- diag(D)
    G <- cbind(diag(D-1), -1)
    Xi <- (upsilon-ntaxa(ps_fido))*G%*%Omega%*%t(G)
    Theta <- matrix(0, D-1, nrow(X))
    Gamma <- diag(nrow(X))
    priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)
    eta_init <- t(driver::alr(t(Y)))
    eta_array <- array(eta_init, dim=c(nrow(eta_init), ncol(eta_init), 2000))
    posterior <- uncollapsePibble(eta_array, priors$X, priors$Theta, priors$Gamma, priors$Xi, priors$upsilon, seed=1)
    dimnames(posterior$Lambda)[[2]] <- rownames(X)
    dimnames(posterior$Lambda)[[1]] <- rownames(Y)[-length(rownames(Y))]
    dimnames(posterior$Sigma)[[1]] <- dimnames(posterior$Sigma)[[2]] <- rownames(Y)[-length(rownames(Y))]
    posterior <- pibblefit(D=D,N=N,Q=nrow(X),coord_system="alr", iter=2000L, alr_base=D,Eta=eta_array,Lambda=posterior$Lambda,Sigma=posterior$Sigma,Y=Y,X=X,names_categories=rownames(Y),names_samples=colnames(Y),names_covariates=rownames(X))
    posterior <- to_clr(posterior)
    fit1 <- posterior
    
    ## Save results of pibble model
    statInfo <- write_pibble(pibble_fit=fit1, focus_var=paste("Cohort", j, sep = ""), siglevel="p95")
    colnames(statInfo) <- c("Index", "LogFC", "CI.L", "CI.UP", "width", "point", "interval", "Taxa")
    statInfo <- statInfo[,c("LogFC", "CI.L", "CI.UP", "Taxa")]
    statInfo$rawP <- NA
    assign(paste(dataset, "_v_", j, "_PIBBLE", sep = ""), statInfo)
    saveRDS(statInfo, paste("results/", cancer, "_PIBBLE/", dataset, "_v_", j, "_PIBBLE.rds", sep = ""))
    rm(posterior)
    rm(priors)
    rm(statInfo)
    all_pibbles <- ls(pat = "*_PIBBLE")
    rm(list = all_pibbles)
}
}
######### Meta analysis of pairwise comparisons #######################################
library(MetaVolcanoR)
library(xlsx)
cancers <- c("Breast", "Melanoma", "Ovary", "Pancreatic", "Prostate", "Colorectal", "Kidney", "Lung")
for (cancer in cancers) {
  cancer_files <- list.files(paste("results/", cancer, "_PIBBLE/", sep = ""))
  for (cancer_rds in cancer_files) {
    pairwise_rds <- readRDS(paste("results/", cancer, "_PIBBLE/", cancer_rds, sep = ""))
    assign(gsub(".rds", "", cancer_rds), pairwise_rds)
  }
  pibble_meta <- mget(ls(pat = "*_PIBBLE"))
  names(pibble_meta) <- gsub("_PIBBLE", "", names(pibble_meta))
  meta_degs_rem <- rem_mv(diffexp=pibble_meta,
                          pcriteria="rawP",
                          foldchangecol='logFC', 
                          genenamecol='Taxa',
                          geneidcol=NULL,
                          collaps=FALSE,
                          llcol='CI.L',
                          rlcol='CI.UP',
                          vcol=NULL, 
                          cvar=TRUE,
                          metathr=0.05,
                          jobname="MetaVolcano",
                          outputfolder=".", 
                          draw='HTML',
                          ncores=1)
  meta_degs_rem@metaresult$AdjP <- p.adjust(meta_degs_rem@metaresult[, "randomP"], method = "fdr")
  write.xlsx(meta_degs_rem@metaresult, "results/PanCancer_AllvsAll.xlsx", sheetName = cancer, col.names = T, row.names = F, append = T)
  all_pibbles <- ls(pat = "*_PIBBLE")
  rm(list = all_pibbles)
}
######### Filter and rank the results #######################################
library(readxl)
library(ComplexHeatmap)
breast <- read_xlsx("results/PanCancer_AllvsAll.xlsx", sheet = 1)
# Keeping 45% of pairwise comparisons in the same direction (1 cancer cohort, 17 control cohorts)
breast <- breast[-which(abs(breast$signcon) < round((17 * 45 / 100))),]
breast <- breast[order(breast$randomP, decreasing = F),]
breast$cancer_rank <- NA 
breast$cancer_rank[which(breast$randomSummary > 0)] <- 1:length(which(breast$randomSummary > 0))
breast$control_rank <- NA 
breast$control_rank[which(breast$randomSummary < 0)] <- -c(1:length(which(breast$randomSummary < 0)))

melanoma <- read_xlsx("results/PanCancer_AllvsAll.xlsx", sheet = 2)
# Keeping 45% of pairwise comparisons in the same direction (11 cancer cohort, 17 control cohorts)
melanoma <- melanoma[-which(abs(melanoma$signcon) < round((187 * 45 / 100))),]
melanoma <- melanoma[order(melanoma$randomP, decreasing = F),]
melanoma$cancer_rank <- NA 
melanoma$cancer_rank[which(melanoma$randomSummary > 0)] <- 1:length(which(melanoma$randomSummary > 0))
melanoma$control_rank <- NA 
melanoma$control_rank[which(melanoma$randomSummary < 0)] <- -c(1:length(which(melanoma$randomSummary < 0)))

ovary <- read_xlsx("results/PanCancer_AllvsAll.xlsx", sheet = 3)
# Keeping 45% of pairwise comparisons in the same direction (1 cancer cohort, 17 control cohorts)
ovary <- ovary[-which(abs(ovary$signcon) < round((17 * 45 / 100))),]
ovary <- ovary[order(ovary$randomP, decreasing = F),]
ovary$cancer_rank <- NA 
ovary$cancer_rank[which(ovary$randomSummary > 0)] <- 1:length(which(ovary$randomSummary > 0))
ovary$control_rank <- NA 
ovary$control_rank[which(ovary$randomSummary < 0)] <- -c(1:length(which(ovary$randomSummary < 0)))

pancreatic <- read_xlsx("results/PanCancer_AllvsAll.xlsx", sheet = 4)
# Keeping 45% of pairwise comparisons in the same direction (3 cancer cohorts, 17 control cohorts)
pancreatic <- pancreatic[-which(abs(pancreatic$signcon) < round((51 * 45 / 100))),]
pancreatic <- pancreatic[order(pancreatic$randomP, decreasing = F),]
pancreatic$cancer_rank <- NA 
pancreatic$cancer_rank[which(pancreatic$randomSummary > 0)] <- 1:length(which(pancreatic$randomSummary > 0))
pancreatic$control_rank <- NA 
pancreatic$control_rank[which(pancreatic$randomSummary < 0)] <- -c(1:length(which(pancreatic$randomSummary < 0)))

prostate <- read_xlsx("results/PanCancer_AllvsAll.xlsx", sheet = 5)
# Keeping 45% of pairwise comparisons in the same direction (3 cancer cohorts, 17 control cohorts)
prostate <- prostate[-which(abs(prostate$signcon) < round((51 * 45 / 100))),]
prostate <- prostate[order(prostate$randomP, decreasing = F),]
prostate$cancer_rank <- NA 
prostate$cancer_rank[which(prostate$randomSummary > 0)] <- 1:length(which(prostate$randomSummary > 0))
prostate$control_rank <- NA 
prostate$control_rank[which(prostate$randomSummary < 0)] <- -c(1:length(which(prostate$randomSummary < 0)))

crc <- read_xlsx("results/PanCancer_AllvsAll.xlsx", sheet = 6)
# Keeping 45% of pairwise comparisons in the same direction (8 cancer cohorts, 17 control cohorts)
crc <- crc[-which(abs(crc$signcon) < round((136 * 45 / 100))),]
crc <- crc[order(crc$randomP, decreasing = F),]
crc$cancer_rank <- NA 
crc$cancer_rank[which(crc$randomSummary > 0)] <- 1:length(which(crc$randomSummary > 0))
crc$control_rank <- NA 
crc$control_rank[which(crc$randomSummary < 0)] <- -c(1:length(which(crc$randomSummary < 0)))

kidney <- read_xlsx("results/PanCancer_AllvsAll.xlsx", sheet = 7)
# Keeping 45% of pairwise comparisons in the same direction (1 cancer cohort, 17 control cohorts)
kidney <- kidney[-which(abs(kidney$signcon) < round((17 * 45 / 100))),]
kidney <- kidney[order(kidney$randomP, decreasing = F),]
kidney$cancer_rank <- NA 
kidney$cancer_rank[which(kidney$randomSummary > 0)] <- 1:length(which(kidney$randomSummary > 0))
kidney$control_rank <- NA 
kidney$control_rank[which(kidney$randomSummary < 0)] <- -c(1:length(which(kidney$randomSummary < 0)))

lung <- read_xlsx("results/PanCancer_AllvsAll.xlsx", sheet = 8)
# Keeping 45% of pairwise comparisons in the same direction (2 cancer cohorts, 17 control cohorts)
lung <- lung[-which(abs(lung$signcon) < round((34 * 45 / 100))),]
lung <- lung[order(lung$randomP, decreasing = F),]
lung$cancer_rank <- NA 
lung$cancer_rank[which(lung$randomSummary > 0)] <- 1:length(which(lung$randomSummary > 0))
lung$control_rank <- NA 
lung$control_rank[which(lung$randomSummary < 0)] <- -c(1:length(which(lung$randomSummary < 0)))

metaphlan_species <- readRDS("data/Cancer_vs_Controls_meta_analysis.rds")
## Join the control ranks from the different cancers
pan_cancer_ranks_controls <- full_join(full_join(full_join(full_join(full_join(full_join(full_join(breast[,c("Taxa", "control_rank")], melanoma[,c("Taxa", "control_rank")], by = "Taxa"), ovary[,c("Taxa", "control_rank")], by = "Taxa"), pancreatic[,c("Taxa", "control_rank")], by = "Taxa"), prostate[,c("Taxa", "control_rank")], by = "Taxa"), crc[,c("Taxa", "control_rank")],  by = "Taxa"), lung[,c("Taxa", "control_rank")],  by = "Taxa"), kidney[,c("Taxa", "control_rank")], by = "Taxa")
pan_cancer_ranks_controls <- as.data.frame(pan_cancer_ranks_controls)
## Find and remove SGBs that were not ranked in at most 2 cancers
remove <- c()
for (i in 1:nrow(pan_cancer_ranks_controls)) {
  row_to_get <- pan_cancer_ranks_controls[i,-1]
  if(length(which(is.na(row_to_get))) >= 3) {
    remove <- append(i, remove)
  }
}
pan_cancer_ranks_controls <- pan_cancer_ranks_controls[-remove,]
colnames(pan_cancer_ranks_controls) <- c("Taxa", "Breast", "Melanoma", "Ovary", "Pancreatic", "Prostate", "CRC", "Lung", "Kidney")
# Get the average control rank from the different cancers not considering unranked SGBs in the max 2 cancers
pan_cancer_ranks_controls$Rank_Mean <- NA
for (i in 1:nrow(pan_cancer_ranks_controls)) {
  row_to_get <- pan_cancer_ranks_controls[i,c(-1,-10)]
  pan_cancer_ranks_controls$Rank_Mean[i] <- sum(row_to_get[,!(is.na(row_to_get))] / length(row_to_get[!(is.na(row_to_get))]))
}
pan_cancer_ranks_controls <- pan_cancer_ranks_controls[order(abs(pan_cancer_ranks_controls$Rank_Mean)),]
pan_cancer_ranks_controls <- pan_cancer_ranks_controls[1:30,]
pan_cancer_ranks_controls$Taxonomy <- apply(tax_table(metaphlan_species)[pan_cancer_ranks_controls$Taxa,c("Family", "Genus", "Species", "Strain")], 1, function(x) {paste(x, collapse = ";")})

## Join the cancer ranks from the different cancers
pan_cancer_ranks_cancer <- full_join(full_join(full_join(full_join(full_join(full_join(full_join(breast[,c("Taxa", "cancer_rank")], melanoma[,c("Taxa", "cancer_rank")], by = "Taxa"), ovary[,c("Taxa", "cancer_rank")], by = "Taxa"), pancreatic[,c("Taxa", "cancer_rank")], by = "Taxa"), prostate[,c("Taxa", "cancer_rank")], by = "Taxa"), crc[,c("Taxa", "cancer_rank")],  by = "Taxa"), lung[,c("Taxa", "cancer_rank")],  by = "Taxa"), kidney[,c("Taxa", "cancer_rank")], by = "Taxa")
pan_cancer_ranks_cancer <- as.data.frame(pan_cancer_ranks_cancer)
## Find and remove SGBs that were not ranked in at most 2 cancers
remove <- c()
for (i in 1:nrow(pan_cancer_ranks_cancer)) {
  row_to_get <- pan_cancer_ranks_cancer[i,-1]
  if(length(which(is.na(row_to_get))) >= 3) {
    remove <- append(i, remove)
  }
}
pan_cancer_ranks_cancer <- pan_cancer_ranks_cancer[-remove,]
colnames(pan_cancer_ranks_cancer) <- c("Taxa", "Breast", "Melanoma", "Ovary", "Pancreatic", "Prostate", "CRC", "Lung", "Kidney")
# Get the average cancer rank from the different cancers not considering unranked SGBs in the max 2 cancers
pan_cancer_ranks_cancer$Rank_Mean <- NA
for (i in 1:nrow(pan_cancer_ranks_cancer)) {
  row_to_get <- pan_cancer_ranks_cancer[i,c(-1,-10)]
  pan_cancer_ranks_cancer$Rank_Mean[i] <- sum(row_to_get[,!(is.na(row_to_get))] / length(row_to_get[!(is.na(row_to_get))]))
}
pan_cancer_ranks_cancer <- pan_cancer_ranks_cancer[order(pan_cancer_ranks_cancer$Rank_Mean),]
pan_cancer_ranks_cancer <- pan_cancer_ranks_cancer[1:30,]
pan_cancer_ranks_cancer$Taxonomy <- apply(tax_table(metaphlan_species)[pan_cancer_ranks_cancer$Taxa,c("Family", "Genus", "Species", "Strain")], 1, function(x) {paste(x, collapse = ";")})

######### Plot the results #######################################
pan_cancer <- rbind(pan_cancer_ranks_cancer, pan_cancer_ranks_controls)
pan_cancer <- as.data.frame(pan_cancer)
row_annot <- rowAnnotation(Average_Rank = anno_text(abs(round(pan_cancer[,10], 2)), location = 0.5, just = "center", show_name = TRUE, gp = gpar(fill = c(rep("#BC544B", 30), rep("#A5D7E9", 30)), border = "black", fontsize = 8, fontface = "bold"), width = max_text_width(abs(round(pan_cancer[,10], 2)))*1.2), show_annotation_name = TRUE)
rownames(pan_cancer) <- pan_cancer$Taxonomy
pan_cancer_plot <- as.matrix(pan_cancer[,c(-1,-10,-11)])
pan_cancer_plot[is.na(pan_cancer_plot)] <- 0
rownames(pan_cancer_plot) <- gsub("f__|g__|s__|t__", "", rownames(pan_cancer_plot))
rownames(pan_cancer_plot) <- gsub("_unclassified", "", rownames(pan_cancer_plot))
Heatmap(matrix = pan_cancer_plot,
        cluster_columns = T,
        cluster_rows = F,
        show_column_names = T,
        rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill, gp = gpar(col = "grey", fill = NA)) {
          if(pan_cancer_plot[i, j] > 0) {
            grid.rect(x, y, width, height, gp = gpar(fill = "#BC544B"))
            grid.text(sprintf("%s", pan_cancer_plot[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "bold"))}
          if(pan_cancer_plot[i, j] < 0) {
            grid.rect(x, y, width, height, gp = gpar(fill = "#A5D7E9"))
            grid.text(sprintf("%s", pan_cancer_plot[i, j] * -1), x, y, gp = gpar(fontsize = 8, fontface = "bold"))}
          if(pan_cancer_plot[i, j] == 0) {
            grid.rect(x, y, width, height, gp = gpar(fill = NA))}
        },
        show_row_names = T, column_dend_height = unit(0.2, "cm"),
        row_split = c(rep("Cancer", 30), rep("Controls", 30)),
        left_annotation = row_annot,
        column_names_gp = gpar(fontsize = 9),
        column_names_rot = 60,
        row_names_gp = gpar(fontsize = 8, fontface = "italic"),
        show_row_dend = F, show_heatmap_legend = F, width = unit(8, "cm"), 
        height = unit(17, "cm"))

