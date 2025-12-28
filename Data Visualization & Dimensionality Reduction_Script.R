

# Analysis of invertebrate diversity and community composition
# in relation to algal cover and algal type across stations and seasons
#
# This script explores:
#  - how invertebrate alpha diversity (Shannon, Simpson) varies with
#    total algal cover and algal richness,
#  - whether algal structure is associated with changes in community
#    composition using NMDS,
#  - and how specific algal types and their coverage constrain
#    invertebrate assemblages using Canonical Correspondence Analysis (CCA).
#
# The analysis combines univariate statistics, unconstrained ordination,
# and constrained ordination to obtain a comprehensive ecological picture.



############################################
## 0. Load packages
############################################

# Assuming the packages are already installed.
# If not, install them once from the console with install.packages().
library(vegan)
library(ggplot2)
library(ggrepel)
library(grid)   # for arrow() in the CCA plot

############################################
## 1. Load data and set up basic variables
############################################

dat <- read.csv(
  "C:/Users/Ariol/Downloads/inv_algae_2stations.csv",
  stringsAsFactors = TRUE
)

# Make sure Station and Season are treated as factors
dat$Station <- factor(dat$Station)
dat$Season  <- factor(dat$Season,
                      levels = c("Winter", "Spring", "Summer", "Autumn"))

# Pick out the invertebrate and algae columns by name pattern
inv_cols <- grep("^Inv", names(dat), value = TRUE)   # invertebrate species
alg_cols <- grep("^Alg", names(dat), value = TRUE)   # algae (% cover)

inv <- dat[, inv_cols, drop = FALSE]
alg <- dat[, alg_cols, drop = FALSE]

# Total invertebrate abundance per sample (sum of all species)
dat$SpeciesAbundance <- rowSums(inv)

# Total algal cover per sample (sum of all algae)
dat$Alg_Total <- rowSums(alg)

# Number of algal taxa present (how many have cover > 0)
dat$TotalAlgaeSpecies <- rowSums(alg > 0)

# Create an ID for each replicate: Station_Season_Rep
dat$Label <- paste(dat$Station, dat$Season, dat$Rep, sep = "_")

############################################
## 2. Alpha diversity (Shannon & Simpson)
############################################

# Calculate Shannon and Simpson diversity for each replicate
dat$Shannon <- diversity(inv, index = "shannon")
dat$Simpson <- diversity(inv, index = "simpson")

############################################
## 3. Correlations and multiple linear model
############################################
## - Correlate Shannon/Simpson with total algal cover
## - Fit LM: Shannon ~ total algal cover + algal richness
############################################

cat("\n===== Shannon vs Total Algal Cover =====\n")
shan_pearson  <- cor.test(dat$Shannon, dat$Alg_Total, method = "pearson")
shan_spearman <- cor.test(dat$Shannon, dat$Alg_Total, method = "spearman")
print(shan_pearson)
print(shan_spearman)

cat("\n===== Simpson vs Total Algal Cover =====\n")
simp_pearson  <- cor.test(dat$Simpson, dat$Alg_Total, method = "pearson")
simp_spearman <- cor.test(dat$Simpson, dat$Alg_Total, method = "spearman")
print(simp_pearson)
print(simp_spearman)

# Linear model: does Shannon change with algal cover and algal richness?
lm_shannon_multi <- lm(Shannon ~ Alg_Total + TotalAlgaeSpecies, data = dat)

cat("\n===== Multiple LM: Shannon ~ Alg_Total + TotalAlgaeSpecies =====\n")
print(summary(lm_shannon_multi))

############################################
## 4. Plots: Shannon & Simpson vs total algal cover
############################################

## 4.1 Shannon vs Alg_Total

shan_r <- unname(shan_pearson$estimate)
shan_p <- shan_pearson$p.value

p_shannon <- ggplot(dat, aes(x = Alg_Total, y = Shannon, label = Label)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  geom_text_repel(size = 3) +
  labs(
    title   = "Shannon vs Total Algal Cover",
    x       = "Total Algal Cover (%)",
    y       = "Shannon Diversity",
    caption = paste0("Pearson r = ", round(shan_r, 3),
                     ", p = ", signif(shan_p, 3),
                     "\nLM: Shannon ~ Alg_Total + TotalAlgaeSpecies")
  ) +
  theme_minimal()

print(p_shannon)

## 4.2 Simpson vs Alg_Total

simp_r <- unname(simp_pearson$estimate)
simp_p <- simp_pearson$p.value

p_simpson <- ggplot(dat, aes(x = Alg_Total, y = Simpson, label = Label)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  geom_text_repel(size = 3) +
  labs(
    title   = "Simpson vs Total Algal Cover",
    x       = "Total Algal Cover (%)",
    y       = "Simpson Diversity",
    caption = paste0("Pearson r = ", round(simp_r, 3),
                     ", p = ", signif(simp_p, 3))
  ) +
  theme_minimal()

print(p_simpson)

############################################
## 5. NMDS (Bray–Curtis) on invertebrate communities
############################################

# Hellinger transform before ordination (standard for community data)
inv_hel <- decostand(inv, method = "hellinger")

set.seed(123)
nmds <- metaMDS(inv_hel,
                distance       = "bray",
                k              = 2,
                trymax         = 200,
                autotransform  = FALSE)

# Check stress in console as a quick quality check
cat("\n===== NMDS Stress =====\n")
print(nmds$stress)
# stressplot(nmds)  # can be used to inspect the fit in more detail

# Extract site scores for plotting
nmds_points <- as.data.frame(scores(nmds, display = "sites"))
names(nmds_points)[1:2] <- c("NMDS1", "NMDS2")

nmds_points$Station           <- dat$Station
nmds_points$Season            <- dat$Season
nmds_points$Alg_Total         <- dat$Alg_Total
nmds_points$TotalAlgaeSpecies <- dat$TotalAlgaeSpecies
nmds_points$Label             <- paste0(nmds_points$Station,
                                        " (", nmds_points$Season, ")")

############################################
## 6. NMDS plots
############################################

## 6.1 NMDS coloured by Season

season_cols <- c(
  "Winter" = "#984EA3",
  "Spring" = "#377EB8",
  "Summer" = "#4DAF4A",
  "Autumn" = "#E41A1C"
)

p_nmds_season <- ggplot(nmds_points,
                        aes(x = NMDS1, y = NMDS2, color = Season)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Label),
                  size = 3,
                  max.overlaps = Inf,
                  show.legend = FALSE) +
  scale_color_manual(values = season_cols, drop = FALSE) +
  labs(
    title = "NMDS Ordination of Invertebrate Communities by Season",
    x     = "NMDS Axis 1",
    y     = "NMDS Axis 2",
    color = "Season"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position   = "right",
    panel.grid.major  = element_line(color = "grey85"),
    panel.grid.minor  = element_blank()
  )

print(p_nmds_season)

## 6.2 NMDS with algal cover and richness overlay

p_nmds_algae <- ggplot(nmds_points,
                       aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Alg_Total,
                 size  = TotalAlgaeSpecies),
             alpha = 0.9) +
  geom_text_repel(aes(label = Label),
                  size = 3,
                  show.legend = FALSE) +
  scale_color_gradient(name = "Algae Coverage (%)",
                       low  = "blue",
                       high = "red") +
  scale_size_continuous(name = "Total Algae Species",
                        range = c(3, 10)) +
  labs(
    title = "NMDS of Invertebrates with Algal Cover and Richness",
    x     = "NMDS Axis 1",
    y     = "NMDS Axis 2"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position   = "right",
    panel.grid.major  = element_line(color = "grey85"),
    panel.grid.minor  = element_blank()
  )

print(p_nmds_algae)

############################################
## 7. CCA: invertebrate communities constrained by algae
############################################

# CCA model: invertebrate community ~ all algae variables
cca_mod <- cca(inv ~ ., data = alg)

cat("\n===== CCA summary =====\n")
print(summary(cca_mod))

############################################
## 8. CCA permutation tests
############################################

cat("\n===== Global CCA test (all algae) =====\n")
print(anova(cca_mod))             # overall test

cat("\n===== CCA test by axis =====\n")
print(anova(cca_mod, by = "axis"))

cat("\n===== CCA test by algae variable (terms) =====\n")
print(anova(cca_mod, by = "terms"))

############################################
## 9. CCA biplot (species, algae, sites)
############################################

# Scores for sites (samples)
sites_df <- as.data.frame(scores(cca_mod, display = "sites"))
sites_df$Station <- dat$Station
sites_df$Season  <- dat$Season
sites_df$SiteLab <- paste0(dat$Station, " (", dat$Season, ")")

# Scores for invertebrate species
spec_df <- as.data.frame(scores(cca_mod, display = "species"))
spec_df$Species <- rownames(spec_df)

# Scores for algae (environmental variables)
env_df <- as.data.frame(scores(cca_mod, display = "bp"))
env_df$Variable <- rownames(env_df)

# Percentage of constrained variance explained by CCA1 and CCA2
eig      <- cca_mod$CCA$eig
var_exp  <- eig / sum(eig)
cca1_pct <- round(var_exp[1] * 100, 1)
cca2_pct <- round(var_exp[2] * 100, 1)

p_cca <- ggplot() +
  # Invertebrate species as red arrows
  geom_segment(data = spec_df,
               aes(x = 0, y = 0,
                   xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.15, "cm")),
               colour = "red",
               linewidth = 0.4,
               alpha = 0.7) +
  geom_text_repel(data = spec_df,
                  aes(x = CCA1, y = CCA2, label = Species),
                  colour = "red",
                  size = 2,
                  max.overlaps = Inf,
                  show.legend = FALSE) +
  
  # Algal variables as green arrows
  geom_segment(data = env_df,
               aes(x = 0, y = 0,
                   xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "darkgreen",
               linewidth = 1.0) +
  geom_text_repel(data = env_df,
                  aes(x = CCA1, y = CCA2, label = Variable),
                  colour = "darkgreen",
                  size = 3.5,
                  max.overlaps = Inf,
                  show.legend = FALSE) +
  
  # Sites coloured by station
  geom_point(data = sites_df,
             aes(x = CCA1, y = CCA2, colour = Station),
             size = 3) +
  geom_text_repel(data = sites_df,
                  aes(x = CCA1, y = CCA2,
                      label = SiteLab,
                      colour = Station),
                  size = 3,
                  max.overlaps = 50,
                  show.legend = FALSE) +
  
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "CCA: Invertebrate Communities by Station\nconstrained by Algal Coverage",
    x     = paste0("CCA1 (", cca1_pct, "% of constrained variance)"),
    y     = paste0("CCA2 (", cca2_pct, "% of constrained variance)"),
    colour = "Station"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position   = "right",
    panel.grid.major  = element_line(color = "grey85"),
    panel.grid.minor  = element_blank()
  )

print(p_cca)

