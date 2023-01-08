#### Meta Analysis with methods i-vi ########
set.seed(1986)
normEdgeR <- function(physeq, method = c('TMM', 'RLE', 'upperquartile')) {
  # require(edgeR)
  otuTab <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq))
  {
    otuTab <- t(otuTab)
  } else {}
  
  if (method == "upperquartile")
  {
    scaledCounts <- t(otuTab) / colSums(otuTab)
    tmpNF <- apply(scaledCounts, MARGIN = 1L, FUN = function(x)
      quantile(x[x != 0], probs = .75))
    normFacts <- tmpNF/exp(mean(log(tmpNF)))
    method <- "UQ"
  } else {
    normFacts <- edgeR:::calcNormFactors(otuTab, method = "TMM")
  }
  if (all(is.na(normFacts)))
  {
    normFacts = sample_sums(physeq)
  }
  physeq@sam_data@.Data <- c(physeq@sam_data@.Data, list(normFacts))
  aux <- physeq@sam_data@names
  aux[length(aux)] <- paste("NF", "TMM", sep = ".")
  physeq@sam_data@names <- aux
  physeq
}
limma_voom <- function(physeq, design = as.formula("~ ORR + age + gender"), normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS")){
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  counts <- as(otu_table(physeq), "matrix")
  if( normFacts=="TSS"){
    NFs = 1
  } else {
    normFacts <- paste("NF", normFacts, sep = ".")
    NFs = get_variable(physeq, normFacts)
    NFs = NFs/exp(mean(log(NFs)))
  }
  
  dge <- DGEList(counts = counts) #, remove.zeros=TRUE)
  dge$samples$norm.factors <- NFs
  design <- model.matrix(design, data.frame(sample_data(physeq)))
  v <- voom(dge, design, plot=FALSE, lib.size = colSums(counts)*NFs)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, n = nrow(dge), sort.by="none", confint = T)
  pval <- tt$P.Value
  padj <- p.adjust(pval,method="BH")
  pValMat <- cbind("rawP" = pval, "adjP" = padj)
  rownames(pValMat) = rownames(tt)
  statInfo <- cbind("logFC" = tt$logFC, "t" = tt$t, "CI.L" = tt$CI.L, "CI.R" = tt$CI.R)
  list("pValMat" = pValMat, "statInfo" = statInfo)
}
limma_voom_zinbweights <- function(physeq, design = as.formula("~ ORR + age + gender"), normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS"), weights) {
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  counts <- as(otu_table(physeq), "matrix")
  if( normFacts=="TSS"){
    NFs = 1
  } else {
    normFacts <- paste("NF", normFacts, sep = ".")
    NFs = get_variable(physeq, normFacts)
    NFs = NFs/exp(mean(log(NFs)))
  }
  
  dge <- DGEList(counts = counts) #, remove.zeros=TRUE)
  dge$samples$norm.factors <- NFs
  design <- model.matrix(design, data.frame(sample_data(physeq)))
  v <- voom(dge, design, plot = FALSE, weights = weights, lib.size = colSums(counts)*NFs)
  v$weights <- v$weights * weights
  fit <- lmFit(v, design, weights = v$weights)
  fit$df.residual <- rowSums(weights) - ncol(design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, n = nrow(dge), sort.by="none", confint = T)
  pval <- tt$P.Value
  padj <- p.adjust(pval,method="BH")
  pValMat <- cbind("rawP" = pval, "adjP" = padj)
  rownames(pValMat) = rownames(tt)
  statInfo <- cbind("logFC" = tt$logFC, "t" = tt$t, "CI.L" = tt$CI.L, "CI.R" = tt$CI.R)
  list("pValMat" = pValMat, "statInfo" = statInfo)
}
normDESeq2 <- function(physeq, whichOTUs = NULL, method = c("poscounts","ratio")){
  # require(DESeq2)
  method <- match.arg(method)
  
  ### Coerce count data to vanilla matrix of integers and check if there are zeroes
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ## select which OTUs to analyse
  if (!missing(whichOTUs) || !is.null(whichOTUs))
  {
    physeq <- prune_taxa(taxa_names(physeq)[whichOTUs], physeq)
  } else
    ## Calculate size factors
    if (method == "poscounts")
    {
      obj <- phyloseq_to_deseq2(physeq, design = ~ ORR + age + gender)
      normFacts <- sizeFactors(DESeq2::estimateSizeFactors(obj,type = "poscounts"))
    } else {
      otuTab <- as(otu_table(physeq), "matrix")
      if (any(otuTab == 0))
      {
        otuTab <- otuTab + 1L
      } else {}
      normFacts <- DESeq2::estimateSizeFactorsForMatrix(otuTab)
    }
  physeq@sam_data@.Data <- c(physeq@sam_data@.Data, list(normFacts))
  aux <- physeq@sam_data@names
  aux[length(aux)] <- paste("NF", method, sep = ".")
  physeq@sam_data@names <- aux
  physeq
}
negBinTestDESeq2 <- function(physeq, design = as.formula("~ ORR + age + gender"), IndepFilter = NULL,
                             normFacts = c("TMM", "RLE", "poscounts","ratio", "CSS", "UQ", "none", "TSS"), returnDispEsts = FALSE) {
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  dds <- phyloseq_to_deseq2(physeq, design = design)
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs = get_variable(physeq, normFacts)
  if(normFacts == "NF.TMM"){
    NFs = NFs * sample_sums(physeq)
  }

  sizeFactors(dds) <- NFs/exp(mean(log(NFs))) #This has no impact on the results but facilitates fitting
  
  ### Run DESeq
  ddsRes <- DESeq(object = dds, test = "LRT", reduced = ~ 1, parallel = FALSE)
  dispEsts <- dispersions(ddsRes)
  ddsRes <- results(ddsRes, alpha = 0.05)
  
  ### Independent Filtering, should be before everything
  if(!is.null(IndepFilter))
  {
    toKeep <- ddsRes$baseMean >= IndepFilter & !is.na(ddsRes$pvalue)
    ddsResFilt <- ddsRes[toKeep, ]
    ddsResFilt$padj <- p.adjust(ddsResFilt$pvalue, method = "BH")
    ddsRes <- as(ddsResFilt, "data.frame")
    ddsRes[order(ddsRes$padj), ]
  } else {}
  
  #  ddsRes$id <- rownames(ddsRes)
  pValMat <- as.matrix(ddsRes[, c("pvalue", "padj")])
  colnames(pValMat) <- c("rawP", "adjP")
  statInfo <- cbind("logFC" = ddsRes$log2FoldChange, "CI.U" = ddsRes$log2FoldChange + ddsRes$lfcSE, "CI.L" = ddsRes$log2FoldChange - ddsRes$lfcSE, "LRT" = ddsRes$stat)
  list("pValMat" = pValMat, "dispEsts" = dispEsts, "statInfo" = statInfo)
}
negBinTestDESeq2_zinbweights <- function(physeq, design = as.formula("~ ORR + age + gender"), IndepFilter = NULL,
                                         normFacts = c("TMM", "RLE", "poscounts", "ratio", "CSS", "UQ", "none", "TSS"), returnDispEsts = FALSE, weights){
  register(SerialParam())
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  dds <- phyloseq_to_deseq2(physeq, design = design)
  
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs = get_variable(physeq, normFacts)
  
  if(normFacts == "NF.TMM"){
    NFs = NFs *sample_sums(physeq)
  }
  
  sizeFactors(dds) <- NFs/exp(mean(log(NFs))) #This has no impact on the results but facilitates fitting
  
  ### ZINB-WaVE weights
  counts <- as(otu_table(physeq), "matrix")
  weights[which(weights<1e-6)] <- 1e-06
  assays(dds)[["weights"]] = weights
  ### Run DESeq
  
  ddsRes <- DESeq(object = dds, test = "LRT", reduced = ~1, parallel = FALSE)
  dispEsts <- dispersions(ddsRes)
  ddsRes <- results(ddsRes, alpha = 0.05)#, cooksCutoff = FALSE)
  
  ### Independent Filtering, should be before everything
  if(!is.null(IndepFilter))
  {
    toKeep <- ddsRes$baseMean >= IndepFilter & !is.na(ddsRes$pvalue)
    ddsResFilt <- ddsRes[toKeep, ]
    ddsResFilt$padj <- p.adjust(ddsResFilt$pvalue, method = "BH")
    ddsRes <- as(ddsResFilt, "data.frame")
    ddsRes[order(ddsRes$padj), ]
  } else {}
  
  #  ddsRes$id <- rownames(ddsRes)
  pValMat <- as.matrix(ddsRes[, c("pvalue", "padj")])
  colnames(pValMat) <- c("rawP", "adjP")
  statInfo <- cbind("logFC" = ddsRes$log2FoldChange, "CI.U" = ddsRes$log2FoldChange + ddsRes$lfcSE, "CI.L" = ddsRes$log2FoldChange - ddsRes$lfcSE, "LRT" = ddsRes$stat)
  list("pValMat" = pValMat,"dispEsts" = dispEsts, "statInfo" = statInfo)
}
computeExactWeights <- function (model, x) {
  mu <- getMu(model)
  pi <- getPi(model)
  theta <- getTheta(model)
  theta <- matrix(rep(theta, each = ncol(x)), ncol = nrow(x))
  nb_part <- dnbinom(t(x), size = theta, mu = mu)
  zinb_part <- pi * ( t(x) == 0 ) + (1 - pi) *  nb_part
  zinbwg <- ( (1 - pi) * nb_part ) / zinb_part
  zinbwg <- t(zinbwg)
  zinbwg[x > 0] <- 1
  zinbwg[zinbwg < 1e-15] <- 1e-15
  zinbwg
}
library(Maaslin2)
library(ANCOMBC)
library(MetaVolcanoR)
library(metafor)
library(reshape2)
library(plyr)
library(limma)
library(curatedMetagenomicData)
library(ComplexHeatmap)
library(edgeR)
library(BiocParallel)
library(zinbwave)
library(DESeq2)
library(phyloseq)
library(xlsx)
metaphlan_species <- readRDS("data/R_vs_NR_meta_analysis.rds")
metaphlan_species2 <- metaphlan_species
taxa <- sweep(otu_table(metaphlan_species), 2, metaphlan_species@sam_data$Number_Reads, FUN = "*") / 100
taxa <- round(taxa, digits = 0)
otu_table(metaphlan_species) <- otu_table(taxa, taxa_are_rows = taxa_are_rows(metaphlan_species))
datasets <- as.character(unique(sample_data(metaphlan_species)$Cohort))
metaphlan_species <- subset_taxa(metaphlan_species, rowSums(otu_table(metaphlan_species)) != 0)
for (j in 1:length(datasets)) {
  print(datasets[j])
  subset_dataset <- subset_samples(metaphlan_species, sample_data(metaphlan_species)$Cohort == datasets[j])
  subset_dataset <- metagMisc::phyloseq_filter_prevalence(subset_dataset, prev.trh=0.05)
  subset_dataset <- subset_taxa(subset_dataset, rowSums(otu_table(subset_dataset)) != 0)
  physeq <- subset_dataset
  physeq@sam_data$ORR <- factor(physeq@sam_data$ORR, levels = c("NR", "R"))
  physeq@sam_data$gender <- factor(physeq@sam_data$gender, levels = c("female", "male"))
  physeq@sam_data$age <- as.numeric(physeq@sam_data$age)
  physeq@sam_data$age <- scale(physeq@sam_data$age)[,1]
  
  physeq@otu_table@.Data[which(physeq@otu_table@.Data>.Machine$integer.max)] <- .Machine$integer.max
  
  physeq <- normDESeq2(physeq = physeq)  # poscounts, similar to RLE
  physeq <- normEdgeR(physeq = physeq, method = "TMM")
  
  zinbmodel <- zinbFit(Y = physeq@otu_table@.Data, 
                       X = model.matrix(~ physeq@sam_data$ORR + physeq@sam_data$age + physeq@sam_data$gender), K = 0,
                       epsilon = 1e10, commondispersion = TRUE, verbose = FALSE, BPPARAM = SerialParam())
  weights <- computeExactWeights(model = zinbmodel, x = physeq@otu_table@.Data)
  colnames(weights) <- colnames(physeq@otu_table)
  rownames(weights) <- rownames(physeq@otu_table)
  DESeq2_poscounts <- negBinTestDESeq2(physeq, normFacts = "poscounts")
  DESeq2_poscounts_zinbwave <- negBinTestDESeq2_zinbweights(physeq, normFacts = "poscounts", weights = weights)
  DESeq2_TMM <- negBinTestDESeq2(physeq, normFacts = "TMM")
  
  DESeq2_poscounts_results <- cbind.data.frame(as.character(rownames(DESeq2_poscounts$pValMat)), DESeq2_poscounts$pValMat, DESeq2_poscounts$statInfo)
  DESeq2_poscounts_zinbwave_results <- cbind.data.frame(as.character(rownames(DESeq2_poscounts_zinbwave$pValMat)), DESeq2_poscounts_zinbwave$pValMat, DESeq2_poscounts_zinbwave$statInfo)
  DESeq2_TMM_results <- cbind.data.frame(as.character(rownames(DESeq2_TMM$pValMat)), DESeq2_TMM$pValMat, DESeq2_TMM$statInfo)
  
  colnames(DESeq2_poscounts_results)[1] <- "Taxa"
  colnames(DESeq2_poscounts_zinbwave_results)[1] <- "Taxa"
  colnames(DESeq2_TMM_results)[1] <- "Taxa"
  
  assign(paste(datasets[j], "_DESeq2_poscounts", sep = ""), DESeq2_poscounts_results)
  assign(paste(datasets[j], "_DESeq2_poscounts_zinbwave", sep = ""), DESeq2_poscounts_zinbwave_results)
  assign(paste(datasets[j], "_DESeq2_TMM", sep = ""), DESeq2_TMM_results)
  
  if (!(datasets[j] %in% c("Barcelona_LeeK_2022", "GopalakrishnanV_2018")) ) {
    limma_voom_TMM_zinbwave <- limma_voom_zinbweights(physeq, normFacts = "TMM", weights = weights)
    limma_voom_TMM_zinbwave_results <- cbind.data.frame(rownames(limma_voom_TMM_zinbwave$pValMat), limma_voom_TMM_zinbwave$pValMat, limma_voom_TMM_zinbwave$statInfo)
    colnames(limma_voom_TMM_zinbwave_results)[1] <- "Taxa"
    assign(paste(datasets[j], "_limma_voom_TMM_zinbwave", sep = ""), limma_voom_TMM_zinbwave_results)
  }
  
  limma_voom_TMM <- limma_voom(physeq, normFacts = "TMM")
  limma_voom_TMM_results <- cbind.data.frame(rownames(limma_voom_TMM$pValMat), limma_voom_TMM$pValMat, limma_voom_TMM$statInfo)
  colnames(limma_voom_TMM_results)[1] <- "Taxa"
  assign(paste(datasets[j], "_limma_voom_TMM", sep = ""), limma_voom_TMM_results)
  
  #ANCOMBC
  pseq <- subset_dataset
  pseq@sam_data$ORR <- factor(pseq@sam_data$ORR, levels = c("NR", "R"))
  pseq@sam_data$gender <- factor(pseq@sam_data$gender, levels = c("female", "male"))
  pseq@sam_data$age <- as.numeric(pseq@sam_data$age)
  pseq@sam_data$age <- scale(pseq@sam_data$age)[,1]
  res <- ancombc(phyloseq =  pseq,
                 formula = "ORR + age + gender", p_adj_method = "fdr",
                 zero_cut = 0.95, lib_cut = 1000, group = "ORR", struc_zero = F)
  statInfo <- cbind("rawP" = res$res$p_val$ORRR, "adjP" = res$res$q_val$ORRR, "LogFC" = res$res$beta$ORRR, "CI.UP" = res$res$beta$ORRR + res$res$se$ORRR, "CI.L" = res$res$beta$ORRR - res$res$se$ORRR)
  statInfo <- as.data.frame(statInfo)
  rownames(statInfo) <- rownames(res$res$p_val)
  statInfo$Taxa <- rownames(res$res$p_val)
  assign(paste(datasets[j], "_ANCOMBC", sep = ""), statInfo)
  
  # Maaslin
  subset_dataset2 <- subset_samples(metaphlan_species2, sample_data(metaphlan_species2)$Cohort == datasets[j])
  subset_dataset2 <- subset_taxa(subset_dataset, rowSums(otu_table(subset_dataset)) != 0)
  subset_dataset2@sam_data$ORR <- factor(subset_dataset2@sam_data$ORR, levels = c("R", "NR"))
  subset_dataset2@sam_data$gender <- factor(subset_dataset2@sam_data$gender, levels = c("female", "male"))
  subset_dataset2@sam_data$age <- as.numeric(subset_dataset2@sam_data$age)
  subset_dataset2@sam_data$age <- scale(subset_dataset2@sam_data$age)[,1]
  
  fit_data <- Maaslin2(
    as(t(subset_dataset2@otu_table), "matrix"), as(subset_dataset2@sam_data, "data.frame"), 'demo_output', transform = "LOGIT",
    fixed_effects = c('ORR', "age", "gender"),
    normalization = 'TSS', correction = "BH", min_prevalence = 0.05,
    standardize = FALSE)
  maaslin2_statInfo <- cbind("rawP" = fit_data$results$pval[fit_data$results$metadata == "ORR"], "adjP" = fit_data$results$qval[fit_data$results$metadata == "ORR"], "LogFC" = fit_data$results$coef[fit_data$results$metadata == "ORR"], "CI.UP" = fit_data$results$coef[fit_data$results$metadata == "ORR"] + fit_data$results$stderr[fit_data$results$metadata == "ORR"], "CI.L" = fit_data$results$coef[fit_data$results$metadata == "ORR"] - fit_data$results$stderr[fit_data$results$metadata == "ORR"])
  maaslin2_statInfo <- as.data.frame(maaslin2_statInfo)
  rownames(maaslin2_statInfo) <- fit_data$results$feature[fit_data$results$metadata == "ORR"]
  maaslin2_statInfo$Taxa <- fit_data$results$feature[fit_data$results$metadata == "ORR"]
  assign(paste(datasets[j], "_maaslin2", sep = ""), maaslin2_statInfo)
}

limma_zb_TMM <- mget(ls(pat = "*_limma_voom_TMM_zinbwave$"))
names(limma_zb_TMM) <- gsub("_limma_voom_TMM_zinbwave", "", names(limma_zb_TMM))

limma_TMM <- mget(ls(pat = "*_limma_voom_TMM$"))
names(limma_TMM) <- gsub("_limma_voom_TMM", "", names(limma_TMM))

deseq2_poscounts <-  mget(ls(pat = "*_DESeq2_poscounts$"))
names(deseq2_poscounts) <- gsub("_DESeq2_poscounts", "", names(deseq2_poscounts))

deseq2_poscounts_zb <- mget(ls(pat = "*_DESeq2_poscounts_zinbwave"))
names(deseq2_poscounts_zb) <- gsub("_DESeq2_poscounts_zinbwave", "", names(deseq2_poscounts_zb))

deseq2_tmm <- mget(ls(pat = "*_DESeq2_TMM"))
names(deseq2_tmm) <- gsub("_DESeq2_TMM", "", names(deseq2_tmm))

ancombc_meta <- mget(ls(pat = "*_ANCOMBC"))
names(ancombc_meta) <- gsub("_ANCOMBC", "", names(ancombc_meta))

maaslin_meta <- mget(ls(pat = "*_maaslin2"))
names(maaslin_meta) <- gsub("_maaslin2", "", names(maaslin_meta))

meta_degs_rem <- rem_mv(diffexp=limma_zb_TMM,
                        pcriteria="rawP",
                        foldchangecol='logFC', 
                        genenamecol='Taxa',
                        geneidcol=NULL,
                        collaps=FALSE,
                        llcol='CI.L',
                        rlcol='CI.R',
                        vcol=NULL, 
                        cvar=TRUE,
                        metathr=0.05,
                        jobname="MetaVolcano",
                        outputfolder=".", 
                        draw='HTML',
                        ncores=1)
meta_degs_rem@metaresult$AdjP <- p.adjust(meta_degs_rem@metaresult[, "randomP"], method = "fdr")
limma_zb_TMM_sig <- meta_degs_rem@metaresult[which(abs(meta_degs_rem@metaresult$signcon) >= 3 ),]
limma_zb_TMM_sig <- limma_zb_TMM_sig[which(limma_zb_TMM_sig$randomP < 0.05),]
write.xlsx(limma_zb_TMM_sig, "results/R_vs_NR_MetaAnalysis.xlsx", sheetName = "limma_voom_TMM_zinbwave", col.names = T, row.names = F, append = T)

meta_degs_rem <- rem_mv(diffexp=limma_TMM,
                        pcriteria="rawP",
                        foldchangecol='logFC', 
                        genenamecol='Taxa',
                        geneidcol=NULL,
                        collaps=FALSE,
                        llcol='CI.L',
                        rlcol='CI.R',
                        vcol=NULL, 
                        cvar=TRUE,
                        metathr=0.05,
                        jobname="MetaVolcano",
                        outputfolder=".", 
                        draw='HTML',
                        ncores=1)
meta_degs_rem@metaresult$AdjP <- p.adjust(meta_degs_rem@metaresult[, "randomP"], method = "fdr")
limma_zb_TMM_sig <- meta_degs_rem@metaresult[which(abs(meta_degs_rem@metaresult$signcon) >= 3 ),]
limma_zb_TMM_sig <- limma_zb_TMM_sig[which(limma_zb_TMM_sig$randomP < 0.05),]
write.xlsx(limma_zb_TMM_sig, "results/R_vs_NR_MetaAnalysis.xlsx", sheetName = "limma_voom_TMM", col.names = T, row.names = F, append = T)

deseq2_poscounts_rem <- rem_mv(diffexp=deseq2_poscounts,
                               pcriteria="rawP",
                               foldchangecol='logFC', 
                               genenamecol='Taxa',
                               geneidcol=NULL,
                               collaps=FALSE,
                               llcol='CI.L',
                               rlcol='CI.U',
                               vcol=NULL, 
                               cvar=TRUE,
                               metathr=0.05,
                               jobname="MetaVolcano",
                               outputfolder=".", 
                               draw='HTML',
                               ncores=1)
deseq2_poscounts_rem@metaresult$AdjP <- p.adjust(deseq2_poscounts_rem@metaresult[, "randomP"], method = "fdr")
deseq2_poscounts_rem_sig <- deseq2_poscounts_rem@metaresult[which(abs(deseq2_poscounts_rem@metaresult$signcon) >= 3 ),]
deseq2_poscounts_rem_sig <- deseq2_poscounts_rem_sig[which(deseq2_poscounts_rem_sig$randomP < 0.05 ),]
write.xlsx(deseq2_poscounts_rem_sig, "results/R_vs_NR_MetaAnalysis.xlsx", sheetName = "deseq2_poscounts", col.names = T, row.names = F, append = T)

deseq2_poscounts_zb_rem <- rem_mv(diffexp=deseq2_poscounts_zb,
                                  pcriteria="rawP",
                                  foldchangecol='logFC', 
                                  genenamecol='Taxa',
                                  geneidcol=NULL,
                                  collaps=FALSE,
                                  llcol='CI.L',
                                  rlcol='CI.U',
                                  vcol=NULL, 
                                  cvar=TRUE,
                                  metathr=0.05,
                                  jobname="MetaVolcano",
                                  outputfolder=".", 
                                  draw='HTML',
                                  ncores=1)
deseq2_poscounts_zb_rem@metaresult$AdjP <- p.adjust(deseq2_poscounts_zb_rem@metaresult[, "randomP"], method = "fdr")
deseq2_poscounts_zb_rem_sig <- deseq2_poscounts_zb_rem@metaresult[which(abs(deseq2_poscounts_zb_rem@metaresult$signcon) >= 3 ),]
deseq2_poscounts_zb_rem_sig <- deseq2_poscounts_zb_rem_sig[which(deseq2_poscounts_zb_rem_sig$randomP < 0.05 ),]
write.xlsx(deseq2_poscounts_zb_rem_sig, "results/R_vs_NR_MetaAnalysis.xlsx", sheetName = "deseq2_poscounts_zb", col.names = T, row.names = F, append = T)

deseq2_tmm_rem <- rem_mv(diffexp=deseq2_tmm,
                         pcriteria="rawP",
                         foldchangecol='logFC', 
                         genenamecol='Taxa',
                         geneidcol=NULL,
                         collaps=FALSE,
                         llcol='CI.L',
                         rlcol='CI.U',
                         vcol=NULL, 
                         cvar=TRUE,
                         metathr=0.05,
                         jobname="MetaVolcano",
                         outputfolder=".", 
                         draw='HTML',
                         ncores=1)
deseq2_tmm_rem@metaresult$AdjP <- p.adjust(deseq2_tmm_rem@metaresult[, "randomP"], method = "fdr")
deseq2_tmm_rem_sig <- deseq2_tmm_rem@metaresult[which(abs(deseq2_tmm_rem@metaresult$signcon) >= 3 ),]
deseq2_tmm_rem_sig <- deseq2_tmm_rem_sig[which(deseq2_tmm_rem_sig$randomP < 0.05),]
write.xlsx(deseq2_tmm_rem_sig, "results/R_vs_NR_MetaAnalysis.xlsx", sheetName = "deseq2_tmm", col.names = T, row.names = F, append = T)

ancombc_meta_rem <- rem_mv(diffexp=ancombc_meta,
                           pcriteria="rawP",
                           foldchangecol='LogFC', 
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
ancombc_meta_rem@metaresult$AdjP <- p.adjust(ancombc_meta_rem@metaresult[, "randomP"], method = "fdr")
ancombc_meta_rem_rem_sig <- ancombc_meta_rem@metaresult[which(abs(ancombc_meta_rem@metaresult$signcon) >= 3 ),]
ancombc_meta_rem_rem_sig <- ancombc_meta_rem_rem_sig[which(ancombc_meta_rem_rem_sig$randomP < 0.05 ),]
write.xlsx(ancombc_meta_rem_rem_sig, "results/R_vs_NR_MetaAnalysis.xlsx", sheetName = "ANCOMBC", col.names = T, row.names = F, append = T)

maaslin_meta_rem <- rem_mv(diffexp=maaslin_meta,
                           pcriteria="rawP",
                           foldchangecol='LogFC', 
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
maaslin_meta_rem@metaresult$AdjP <- p.adjust(maaslin_meta_rem@metaresult[, "randomP"], method = "fdr")
maaslin_meta_rem_sig <- maaslin_meta_rem@metaresult[which(abs(maaslin_meta_rem@metaresult$signcon) >= 3 ),]
maaslin_meta_rem_sig <- maaslin_meta_rem_sig[which(maaslin_meta_rem_sig$randomP < 0.05 ),]
write.xlsx(maaslin_meta_rem_sig, "results/R_vs_NR_MetaAnalysis.xlsx", sheetName = "Maaslin2_Logit", col.names = T, row.names = F, append = T)

#### Meta Analysis with SMD (method vii) ########
library(metafor)
library(reshape2)
library(plyr)
metaphlan_species <- readRDS("data/R_vs_NR_meta_analysis.rds")
metaphlan_species <- metagMisc::phyloseq_filter_prevalence(metaphlan_species, prev.trh=0.05)
metaphlan_species <- prune_taxa(taxa_sums(metaphlan_species)>0, metaphlan_species)
taxa <- apply(otu_table(metaphlan_species), 2, function(x) asin(sqrt(x/sum(x))))
otu_table(metaphlan_species) <- otu_table(taxa, taxa_are_rows = taxa_are_rows(metaphlan_species))
meta_analysis_results <- matrix(NA, ncol = 9, nrow=nrow(otu_table(metaphlan_species)))
colnames(meta_analysis_results) <- c("Species", "Number_Samples", "I^2", "P-value Fit", "Pvalue", "Standard Error", "Coefficient", "Confidence Interval Lower Limit", "Confidence Interval Upper Limit")
meta_analysis_results <- as.data.frame(meta_analysis_results)
for (k in 1:nrow(otu_table(metaphlan_species))) {
  species <- otu_table(metaphlan_species)[k,]
  species_table <- matrix(NA, ncol = 8, nrow = length(unique(metaphlan_species@sam_data$Cohort)))
  colnames(species_table) <- c("Study", "Dataset", "SampleSize_Cancer", "Mean_Cancer", "SD_Cancer", "SampleSize_Control", "Mean_Control", "SD_Control")
  species_table <- as.data.frame(species_table)
  datasets <- as.character(unique(sample_data(metaphlan_species)$Cohort))
  for (i in 1:length(datasets)) {
    dataset <- datasets[i]
    species_table$Study[i] <- i
    species_table$Dataset[i] <- dataset
    species_table$SampleSize_Cancer[i] <- length(which(sample_data(metaphlan_species)$ORR == "R" & sample_data(metaphlan_species)$Cohort == dataset))
    species_table$Mean_Cancer[i] <- mean(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$ORR == "R" & sample_data(metaphlan_species)$Cohort == dataset]])
    species_table$SD_Cancer[i] <- sd(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$ORR == "R" & sample_data(metaphlan_species)$Cohort == dataset]])
    species_table$SampleSize_Control[i] <- length(which(sample_data(metaphlan_species)$ORR == "NR" & sample_data(metaphlan_species)$Cohort == dataset))
    species_table$Mean_Control[i] <- mean(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$ORR == "NR" & sample_data(metaphlan_species)$Cohort == dataset]])
    species_table$SD_Control[i] <- sd(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$ORR == "NR" & sample_data(metaphlan_species)$Cohort == dataset]])
  }
  dat1 <- escalc(measure="SMD", m1i=Mean_Cancer, sd1i=SD_Cancer, n1i=SampleSize_Cancer,
                 m2i=Mean_Control, sd2i=SD_Control, n2i=SampleSize_Control, data=species_table, vtype="UB")
  if (length(which((dat1$Mean_Cancer == dat1$Mean_Control) == TRUE)) > 3) {
    meta_analysis_results$Species[k] <- rownames(otu_table(metaphlan_species))[k]
    meta_analysis_results$Number_Samples[k] <- 0  
    next()
  }
  else{
    res1 <- rma(yi, vi, data=dat1, control=list(stepadj=0.5, maxiter = 10000), verbose=TRUE, digits=5)
    meta_analysis_results$Species[k] <- rownames(otu_table(metaphlan_species))[k]
    meta_analysis_results$Number_Samples[k] <- sum(res1$ni)
    meta_analysis_results$`I^2`[k] <- res1$I2
    meta_analysis_results$`P-value Fit`[k] <- res1$QEp
    meta_analysis_results$Pvalue[k] <- res1$pval
    meta_analysis_results$`Standard Error`[k] <- res1$b
    meta_analysis_results$Coefficient[k] <- res1$b
    meta_analysis_results$`Confidence Interval Lower Limit`[k] <- res1$ci.lb
    meta_analysis_results$`Confidence Interval Upper Limit`[k] <- res1$ci.ub
  }
}
meta_analysis_results$PvalueAdjusted <- p.adjust(meta_analysis_results$Pvalue, method = "fdr")
meta_analysis_results <- meta_analysis_results[order(abs(meta_analysis_results$Coefficient), decreasing = T), ]
rownames(meta_analysis_results) <- meta_analysis_results$Species
meta_analysis_results <- meta_analysis_results[!(is.na(meta_analysis_results$`I^2`)),]
meta_analysis_results_filt <- meta_analysis_results[meta_analysis_results$Pvalue < 0.05,]
write.xlsx(meta_analysis_results_filt, "results/R_vs_NR_MetaAnalysis.xlsx", sheetName = "SMD", col.names = T, row.names = F, append = T)

#### Meta Analysis with Pibble (method viii) ########
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
metaphlan_species <- readRDS("data/R_vs_NR_meta_analysis.rds")
metaphlan_species@sam_data$SampleID <- rownames(metaphlan_species@sam_data)
datasets <- as.character(unique(sample_data(metaphlan_species)$Cohort))
for (j in 1:length(datasets)) {
  if(datasets[j] %in% c("WindTT_2020", "Barcelona_LeeK_2022", "Leeds_LeeK_2022")) {
    next()
  }
  else {
    subset_dataset <- subset_samples(metaphlan_species, sample_data(metaphlan_species)$Cohort == datasets[j])
    subset_dataset <- subset_taxa(subset_dataset, rowSums(otu_table(subset_dataset)) != 0)
    subset_dataset <- metagMisc::phyloseq_filter_prevalence(subset_dataset, prev.trh=0.05)
    subset_dataset@sam_data$age <- as.numeric(subset_dataset@sam_data$age)
    subset_dataset@sam_data$age <- scale(subset_dataset@sam_data$age)[,1]
    subset_dataset@sam_data$SampleID <- rownames(subset_dataset@sam_data)
    mdat_fido <-
      as(sample_data(subset_dataset), "data.frame") %>% 
      mutate(SampleID = as.character(SampleID),
             gender = relevel(factor(gender, ordered = FALSE), ref="male"), 
             age = age,
             ORR = relevel(factor(ORR, ordered = FALSE), ref="NR")) %>% 
      select(SampleID, ORR, gender, age)
    ps_fido <- phyloseq(otu_table(subset_dataset, taxa_are_rows=T),
                        sample_data(mdat_fido)) %>% 
      metagMisc::phyloseq_filter_prevalence(prev.trh=0.05)
    
    ps_fido <- prune_taxa(taxa_sums(ps_fido) > 0, ps_fido)
    ps_fido <- prune_samples(sample_sums(ps_fido) > 0, ps_fido)
    f <- reformulate(termlabels=c("ORR", "age", "gender"))
    X <- t(model.matrix(f, data=as(sample_data(ps_fido), "data.frame")))
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
    
    # Clean up output
    dimnames(posterior$Lambda)[[2]] <- rownames(X)
    dimnames(posterior$Lambda)[[1]] <- rownames(Y)[-length(rownames(Y))]
    dimnames(posterior$Sigma)[[1]] <- dimnames(posterior$Sigma)[[2]] <- rownames(Y)[-length(rownames(Y))]
    posterior <- pibblefit(D=D,
                           N=N,
                           Q=nrow(X),
                           coord_system="alr",
                           iter=2000L,
                           alr_base=D,
                           Eta=eta_array,
                           Lambda=posterior$Lambda,
                           Sigma=posterior$Sigma,
                           Y=Y,
                           X=X,
                           names_categories=rownames(Y),
                           names_samples=colnames(Y),
                           names_covariates=rownames(X))
    posterior <- to_clr(posterior)
    fit1 <- posterior
    statInfo <- write_pibble(pibble_fit=fit1, focus_var="ORRR", siglevel="p95")
    colnames(statInfo) <- c("Index", "LogFC", "CI.L", "CI.UP", "width", "point", "interval", "Taxa")
    statInfo <- statInfo[,c("LogFC", "CI.L", "CI.UP", "Taxa")]
    statInfo$rawP <- NA
    assign(paste(datasets[j], "_PIBBLE", sep = ""), statInfo)
    print(datasets[j])
  }
}
pibble_meta <- mget(ls(pat = "*_PIBBLE"))
names(pibble_meta) <- gsub("_PIBBLE", "", names(pibble_meta))
saveRDS(pibble_meta, "results/Meta_analysis_ORR_pibble.rds")

library(MetaVolcanoR)
library(xlsx)
pibble_meta <- readRDS("results/Meta_analysis_ORR_pibble.rds")

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
pibble_sig <- meta_degs_rem@metaresult[which(abs(meta_degs_rem@metaresult$signcon) >= 3 ),]
pibble_sig <- pibble_sig[which(pibble_sig$randomP < 0.05),]
write.xlsx(pibble_sig, "results/R_vs_NR_MetaAnalysis.xlsx", sheetName = "Pibble", col.names = T, row.names = F, append = T)


#### Read in and filter the results ########
library(metafor)
library(reshape2)
library(plyr)
library("readxl")
library(curatedMetagenomicData)

limma_voom_TMM_zinbwave <- read_excel("results/R_vs_NR_MetaAnalysis.xlsx", sheet = 1)
limma_voom_TMM <- read_excel("results/R_vs_NR_MetaAnalysis.xlsx", sheet = 2)
DESeq2_poscounts <- read_excel("results/R_vs_NR_MetaAnalysis.xlsx", sheet = 3)
DESeq2_poscounts_zinbwave <- read_excel("results/R_vs_NR_MetaAnalysis.xlsx", sheet = 4)
DESeq2_TMM <- read_excel("results/R_vs_NR_MetaAnalysis.xlsx", sheet = 5)
ANCOMBC <- read_excel("results/R_vs_NR_MetaAnalysis.xlsx", sheet = 6)
maaslin <- read_excel("results/R_vs_NR_MetaAnalysis.xlsx", sheet = 7)
SMD <- read_excel("results/R_vs_NR_MetaAnalysis.xlsx", sheet = 8)
pibble <- read_excel("results/R_vs_NR_MetaAnalysis.xlsx", sheet = 9)
limma_voom_TMM_zinbwave <- limma_voom_TMM_zinbwave[abs(limma_voom_TMM_zinbwave$signcon) >= 4,]
limma_voom_TMM <- limma_voom_TMM[abs(limma_voom_TMM$signcon) >= 4,]
DESeq2_poscounts <- DESeq2_poscounts[abs(DESeq2_poscounts$signcon) >= 4,]
DESeq2_poscounts_zinbwave <- DESeq2_poscounts_zinbwave[abs(DESeq2_poscounts_zinbwave$signcon) >= 4,]
DESeq2_TMM <- DESeq2_TMM[abs(DESeq2_TMM$signcon) >= 4,]
ANCOMBC <- ANCOMBC[abs(ANCOMBC$signcon) >= 4,]
maaslin <- maaslin[abs(maaslin$signcon) >= 4,]
pibble <- pibble[abs(pibble$signcon) >= 4,]
listInput_crc <- list(DESeq2_poscounts = DESeq2_poscounts$Taxa, 
                      DESeq2_poscounts_zinbwave = DESeq2_poscounts_zinbwave$Taxa, 
                      DESeq2_TMM = DESeq2_TMM$Taxa, 
                      limma_voom_TMM_zinbwave = limma_voom_TMM_zinbwave$Taxa,
                      limma_voom_TMM = limma_voom_TMM$Taxa, 
                      SMD = SMD$Species, 
                      ANCOMBC = ANCOMBC$Taxa, 
                      Maaslin2 = maaslin$Taxa,
                      Pibble = pibble$Taxa)
allsig_taxa_up <- table(c(ANCOMBC$Taxa[ANCOMBC$signcon > 0],
                          DESeq2_poscounts_zinbwave$Taxa[DESeq2_poscounts_zinbwave$signcon > 0],
                          DESeq2_poscounts$Taxa[DESeq2_poscounts$signcon > 0],
                          DESeq2_TMM$Taxa[DESeq2_TMM$signcon > 0],
                          limma_voom_TMM_zinbwave$Taxa[limma_voom_TMM_zinbwave$signcon > 0],
                          limma_voom_TMM$Taxa[limma_voom_TMM$signcon > 0],
                          SMD$Species[SMD$Coefficient > 0],
                          maaslin$Taxa[maaslin$signcon > 0],
                          pibble$Taxa[pibble$signcon > 0]))
allsig_taxa_up <- allsig_taxa_up[allsig_taxa_up >= 3]
allsig_taxa_up <- names(sort(allsig_taxa_up, decreasing = T))
allsig_taxa_down <- table(c(ANCOMBC$Taxa[ANCOMBC$signcon < 0],
                            DESeq2_poscounts_zinbwave$Taxa[DESeq2_poscounts_zinbwave$signcon < 0],
                            DESeq2_poscounts$Taxa[DESeq2_poscounts$signcon < 0],
                            DESeq2_TMM$Taxa[DESeq2_TMM$signcon < 0],
                            limma_voom_TMM_zinbwave$Taxa[limma_voom_TMM_zinbwave$signcon < 0],
                            limma_voom_TMM$Taxa[limma_voom_TMM$signcon < 0],
                            SMD$Species[SMD$Coefficient < 0], 
                            maaslin$Taxa[maaslin$signcon < 0],
                            pibble$Taxa[pibble$signcon < 0]))
allsig_taxa_down <- allsig_taxa_down[allsig_taxa_down >= 3]
allsig_taxa_down <- names(sort(allsig_taxa_down, decreasing = T))
allsig_taxa <- c(allsig_taxa_up, allsig_taxa_down)

#### Plot the results ########
library(xlsx)
metaphlan_species <- readRDS("data/R_vs_NR_meta_analysis.rds")
metaphlan_species <- subset_taxa(metaphlan_species, rownames(metaphlan_species@otu_table) %in% allsig_taxa)
meta_analysis_results <- matrix(NA, ncol = length(unique(metaphlan_species@sam_data$Cohort)), nrow=nrow(otu_table(metaphlan_species)))
colnames(meta_analysis_results) <- unique(metaphlan_species@sam_data$Cohort)
meta_analysis_results <- as.data.frame(meta_analysis_results)
datasets <- unique(metaphlan_species@sam_data$Cohort)
for (k in 1:nrow(otu_table(metaphlan_species))) {
  species <- otu_table(metaphlan_species)[k,]
  for (j in 1:length(datasets)) {
    meta_analysis_results[k,datasets[j]] <- wilcox.test(as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$ORR == "R" & sample_data(metaphlan_species)$Cohort == datasets[j]]]),
                                                        as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$ORR == "NR" & sample_data(metaphlan_species)$Cohort == datasets[j]]]))$p.value
    
    if (mean(as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$ORR == "R" & sample_data(metaphlan_species)$Cohort == datasets[j]]])) < mean(as.numeric(species[,rownames(sample_data(metaphlan_species))[sample_data(metaphlan_species)$ORR == "NR" & sample_data(metaphlan_species)$Cohort == datasets[j]]]))) {
      meta_analysis_results[k,datasets[j]] <- meta_analysis_results[k,datasets[j]]  * -1
      
    }
  }
}
rownames(meta_analysis_results) <- rownames(otu_table(metaphlan_species))
meta_analysis_results <- meta_analysis_results[rownames(meta_analysis_results) %in% allsig_taxa,]
meta_analysis_results <- meta_analysis_results[allsig_taxa,]
library(circlize)
library(ComplexHeatmap)
prev <- apply(otu_table(metaphlan_species)[allsig_taxa,], 1, function(x) (length(which(x > 0)) / length(x)) * 100 )
meta_analysis_results[meta_analysis_results == 'NaN'] <- NA
meta_analysis_results[is.na(meta_analysis_results)] <- 99
meta_analysis_results <- as.matrix(meta_analysis_results)
col_fun = circlize::colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
anno_df = data.frame(DESeq2_poscounts = allsig_taxa %in% listInput_crc$DESeq2_poscounts,
                     DESeq2_poscounts_zinbwave = allsig_taxa %in% listInput_crc$DESeq2_poscounts_zinbwave,
                     DESeq2_TMM = allsig_taxa %in% listInput_crc$DESeq2_TMM, 
                     limma_voom_TMM_zinbwave = allsig_taxa %in% listInput_crc$limma_voom_TMM_zinbwave,
                     limma_voom_TMM = allsig_taxa %in% listInput_crc$limma_voom_TMM,
                     SMD = allsig_taxa %in% listInput_crc$SMD, 
                     ANCOMBC = allsig_taxa %in% listInput_crc$ANCOMBC,
                     Maaslin2 = allsig_taxa %in% listInput_crc$Maaslin2,
                     Pibble = allsig_taxa %in% listInput_crc$Pibble)
anno_df[anno_df == TRUE] <- "1"
anno_df[anno_df == FALSE] <- "0"
rownames(anno_df) <- allsig_taxa
anno_df <- anno_df[c(allsig_taxa_up, allsig_taxa_down),]
for (i in 1:nrow(meta_analysis_results)) {
  rownames(meta_analysis_results)[i] <- paste(as.character(tax_table(metaphlan_species)[rownames(meta_analysis_results)[i],c("Family", "Genus", "Species", "Strain")]), collapse = ";")
}
rownames(meta_analysis_results) <- gsub("f__|g__|s__|t__", "", rownames(meta_analysis_results))
rownames(meta_analysis_results) <- gsub("_unclassified", "", rownames(meta_analysis_results))
row_annot <- rowAnnotation(Prevalence = anno_barplot(prev, axis_param = list(direction = "reverse", side = "top"), 
                                                     width = unit(1.2, "cm"), height = unit(0.8, "cm")),
                           df = anno_df,
                           col = list(DESeq2_poscounts = c("1"="black", "0"="#FFFFFF"), 
                                      DESeq2_poscounts_zinbwave =c("1"="black", "0"="#FFFFFF"), 
                                      DESeq2_TMM =c("1"="black", "0"="#FFFFFF"), 
                                      limma_voom_TMM_zinbwave =c("1"="black", "0"="#FFFFFF"), 
                                      limma_voom_TMM =c("1"="black", "0"="#FFFFFF"), 
                                      SMD =c("1"="black", "0"="#FFFFFF"), 
                                      ANCOMBC =c("1"="black", "0"="#FFFFFF"), 
                                      Maaslin2 = c("1"="black", "0"="#FFFFFF"),
                                      Pibble = c("1"="black", "0"="#FFFFFF")),
                           simple_anno_size = unit(0.15, "cm"), 
                           gp = gpar(fontsize = 2), 
                           annotation_name_gp = gpar(fontsize = 6), 
                           show_legend = FALSE, gap = unit(c(1, rep(0.3, 8)), "mm"))

colnames(meta_analysis_results) <- gsub("RoutyB_2022_Lung", "RoutyB_2018_Lung", colnames(meta_analysis_results))
Heatmap(matrix = meta_analysis_results,
        cluster_columns = T,
        cluster_rows = T,
        show_column_names = T,
        rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill, gp = gpar(col = "grey", fill = NA)) {
          if(meta_analysis_results[i, j] > 0 & meta_analysis_results[i, j] > 0.05 & meta_analysis_results[i, j] < 99) {
            grid.rect(x, y, width, height, gp = gpar(fill = "#BC544B"))}
          if(meta_analysis_results[i, j] < 0 & meta_analysis_results[i, j] < -0.05 ) {
            grid.rect(x, y, width, height, gp = gpar(fill = "#A5D7E9"))}
          if(meta_analysis_results[i, j] == 99) {
            grid.rect(x, y, width, height, gp = gpar(fill = NA))}
          if(meta_analysis_results[i, j] < 0.05 & meta_analysis_results[i, j] > 0 & meta_analysis_results[i, j] < 99) {
            grid.rect(x, y, width, height, gp = gpar(fill = "#BC544B"))
            grid.text(sprintf("%.3f", meta_analysis_results[i, j]), x, y, gp = gpar(fontsize = 7, fontface = "bold"))}
          if(meta_analysis_results[i, j] > -0.05 & meta_analysis_results[i, j] < 0) {
            grid.rect(x, y, width, height, gp = gpar(fill = "#A5D7E9"))
            grid.text(sprintf("%.3f", meta_analysis_results[i, j] * -1), x, y, gp = gpar(fontsize = 7, fontface = "bold"))}
        },
        show_row_names = T, name = "P-value", column_dend_height = unit(0.2, "cm"),
        row_split = c(rep("R", 46), rep("NR", 15)),
        left_annotation = row_annot,
        column_names_gp = gpar(fontsize = 8),
        column_names_rot = 90,
        row_names_gp = gpar(fontsize = 7.7, fontface = "italic"),
        show_row_dend = F, show_heatmap_legend = FALSE, width = unit(7.8, "cm"), 
        height = unit(16, "cm"))