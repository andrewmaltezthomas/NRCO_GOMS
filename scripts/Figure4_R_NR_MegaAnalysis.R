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
metaphlan_species <- readRDS("data/R_vs_NR_meta_analysis.rds")
metaphlan_species@sam_data$age <- as.numeric(metaphlan_species@sam_data$age)
metaphlan_species@sam_data$age <- scale(metaphlan_species@sam_data$age)[,1]
metaphlan_species@sam_data$gender <- factor(metaphlan_species@sam_data$gender, levels = c("female", "male"))
metaphlan_species@sam_data$SampleID <- rownames(metaphlan_species@sam_data)
table(metaphlan_species@sam_data$Cohort)
mdat_fido <-
  as(sample_data(metaphlan_species), "data.frame") %>% 
  mutate(SampleID = as.character(SampleID),
         ORR = relevel(factor(ORR, ordered = FALSE), ref="NR"),
         Cohort = Cohort,
         gender = gender, 
         age = age) %>% 
  select(SampleID, ORR, Cohort, gender, age)

# Filter SGBs with less than 5% abundance
ps_fido <- phyloseq(otu_table(metaphlan_species, taxa_are_rows=T),
                    sample_data(mdat_fido), 
                    tax_table(metaphlan_species)) %>% 
  metagMisc::phyloseq_filter_prevalence(prev.trh=0.05)

ps_fido <- prune_taxa(taxa_sums(ps_fido) > 0, ps_fido)
ps_fido <- prune_samples(sample_sums(ps_fido) > 0, ps_fido)

f <- reformulate(termlabels=c("ORR", "Cohort", "age", "gender"))
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

# Fitting the model
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
for (i in 1:nrow(fit1$Y)) {
  rownames(fit1$Y)[i] <- paste(as.character(tax_table(metaphlan_species)[rownames(fit1$Y)[i],c("Family", "Genus", "Species", "Strain")]), collapse = ";")
}
p3 <- as.data.frame.table(fit1$Lambda) %>% 
  filter(Var2=="ORRR") %>% 
  mutate(idx=as.numeric(as.factor(Var1))) %>% 
  select(idx, Freq, x=Var2, y=Var1) %>% 
  group_by(idx) %>%
  ggdist::median_qi(Freq, .width=c(0.5, 0.75, 0.90, 0.95, 0.97)) %>%
  mutate(species=rep(rownames(fit1$Y),5))
saveRDS(p3, "results/Mega_analysis_ORR.rds")

fit1 <- readRDS("results/Mega_analysis_ORR_2.rds")
fit1$species <- gsub("f__|g__|s__|t__", "", fit1$species)
fit1$species <- gsub("_unclassified", "", fit1$species)
sig_decreasing <- fit1 %>%
  mutate(species=reorder(factor(species),Freq),) %>% 
  pivot_wider(species, names_from=.width, values_from=.upper) %>%
  dplyr::select(species, p50=`0.5`, p75=`0.75`, p90=`0.9`, p95=`0.95`, p97=`0.97`) %>%
  filter(!! rlang::sym('p95') < -0.2) %>%
  mutate(species=factor(species)) %>% 
  dplyr::select(species) %>% 
  pull() %>% 
  levels()

sig_increasing <- fit1 %>%
  mutate(species=reorder(factor(species),Freq),) %>% 
  pivot_wider(species, names_from=.width, values_from=.lower) %>%
  dplyr::select(species, p50=`0.5`, p75=`0.75`, p90=`0.9`, p95=`0.95`, p97=`0.97`) %>%
  filter(!! rlang::sym('p95') > 0.2) %>%
  mutate(species=factor(species)) %>% 
  dplyr::select(species) %>% 
  pull() %>% 
  levels()

p <- fit1 %>% 
  filter(species %in% c(sig_decreasing, sig_increasing)) %>%
  ggplot(aes(y=reorder(factor(species),Freq), x=Freq, xmin=.lower, xmax=.upper)) +
  ggdist::geom_interval(aes(alpha=.width), color="orange3") +
  scale_alpha_continuous("Credible interval", range=c(.7, .15), breaks=c(0.5, 0.75, 0.90, 0.95, 0.97)) +
  geom_point() +
  theme(
    legend.key=element_rect(fill='white'),
    legend.text=element_text(size=10, color="black"),
    legend.position=c(0.5, 0,8),
    strip.background=element_blank(),
    strip.text=element_blank(),
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(), 
    panel.background=element_rect(fill="white"),
    panel.border=element_rect(colour="black", fill=NA, size=1),
    axis.ticks.length=unit(0.25,"cm"), 
    axis.text.x=element_text(size=8.5, color="black"),
    axis.text.y=element_text(size=7.5, color="black")) +
  labs(x="Log-Ratio Value", y=NULL) +
  geom_vline(xintercept=0, linetype="dashed", color="darkgray")
ggsave(plot = p, filename = "results/Mega_analysis_ORR.pdf", width = 10, height = 5, device = "pdf")
  