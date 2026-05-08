# Bee communities and restoration
# Charles Thrift
# 4 February 2026

setwd("~/bee_communities/scripts")

# 0. Libraries
library(tidyverse)
library(rareNMtests)
library(vegan)

# 0. Data
## 0A. Bee Occurrences (All)
occ <- read.csv("../data/bee_occ.csv")
## 02. Bee Species Lists By Site (Including traits)
bee_list <- read.csv("../data/bee_lists.csv")
## 03. Species Matrix (Bees only)
sp_mat_bee <- read.csv("../data/sp_matrix_bees.csv")
## 04. Species Matrix (Plants only)
sp_mat_plant <- read.csv("../data/sp_matrix_plants.csv")
## 05. Species Matrix (Bees and plants together)
sp_mat_all <- read.csv("../data/sp_matrix_all.csv")
## 06. Plant Lists By Site (Including traits)
plant_list <- read.csv("../data/plant_lists.csv")


# 1. Overview 
## Spatial Autocorrelation of Pan Traps
### Input: Species matrix (bees) or the bee occurrences?
library(ape)
spatial <- read.csv("../data/spatial.csv")
spatial_mat <- as.matrix(dist(cbind(spatial$lat, spatial$long)))

spatial_mat.inv <- 1/spatial_mat
diag(spatial_mat.inv) <- 0
spatial_mat.inv[1:5, 1:5]

moran_rich <- Moran.I(spatial$beeRichness, spatial_mat.inv)
moran_rich$observed # Moran's I observed value
moran_rich$expected # Moran's I expected value
(moran_rich$observed - moran_rich$expected) / moran_rich$sd # Z-score
moran_rich$p.value # p-value

moran_divers <- Moran.I(spatial$beeDiversity, spatial_mat.inv)
moran_divers$observed # Moran's I value
moran_divers$expected # Moran's I expected value
(moran_divers$observed - moran_divers$expected) / moran_divers$sd # Z-score
moran_divers$p.value # p-value

moran_abund <- Moran.I(spatial$beeAbundance, spatial_mat.inv)
moran_abund$observed # Moran's I value
moran_abund$expected # Moran's I expected value
(moran_abund$observed - moran_abund$expected) / moran_abund$sd # Z-score
moran_abund$p.value # p-value


# 2. Question 1: Do bee communities vary depending on plant communities?

## 2a. Test similarity by site

### 2a.i. Visualize similarity by site: Euler diagrams
#### Input: Bee species lists by site; plant species lists by site
#### Output: figures
library(eulerr)
bee_euler <- bee_list[,c(2:5)]
bee_euler <- bee_euler %>% rename("UCSB Lagoon" = "UCSB_Lagoon",
                                  "Coal Oil Point Reserve" = "Coal_Oil_Point_Reserve",
                                  "North Campus Open Space" = "North_Campus_Open_Space")
bee_plot <- plot(euler(bee_euler[, 2:4], shape = "ellipse"), 
                 quantities = list(type = c("counts"),col="black",fontsize=12),
                 labels=list(col="black", fontsize = 12),
                 edges = T,
                 main = c("(A) Bee Species",col="black"))

bee_plot

ggsave(bee_plot, 
       filename = "../figures/euler_bee.jpeg",
       device = "jpeg",
       width=10,units="in")

plant_euler <- plant_list %>% 
  distinct(scientificName, restoration_site) %>% 
  mutate(present = 1) %>% 
  pivot_wider(
    names_from = restoration_site,
    values_from = present,
    values_fill = 0
  ) %>% 
  relocate("UCSB Lagoon", .after = "North Campus Open Space")
plant_plot <- plot(euler(plant_euler[, 2:4], shape = "ellipse"), 
                 quantities = list(type = c("counts"),col="black",fontsize=12),
                 labels=list(col="black", fontsize = 12),
                 edges = T,
                 main = c("(B) Plant Species",col="black"))

plant_plot

plots <- cowplot::plot_grid(bee_plot, plant_plot, ncol = 2, nrow = 1,
                   scale = 0.85)
plots

ggsave(plots, 
       filename = "../figures/euler_both.jpeg",
       width=28,
       units="cm",
       dpi = 300)



### 2a.ii. Measure Sorenson similarity by site
#### Input: Bee species lists by site; plant species lists by site
#### Output: values, plus matrix (table)
ncos_subset <- subset(bee_list, North_Campus_Open_Space > 0)
lago_subset <- subset(bee_list, UCSB_Lagoon > 0)
copr_subset <- subset(bee_list, Coal_Oil_Point_Reserve > 0)

# Sorenson's Beta is (2a)/(2a + b + c)
# where a is common taxa
# and b is taxa exclusive to site 1
# and c is taxa exclusive to site 2

# FIRST, NCOS and COPR
# a_NC is common to both NCOS and COPR
# b_NC is taxa exclusive to NCOS
# c_NC is taxa exclusive to COPR

a_NC <- length(intersect(ncos_subset$species, copr_subset$species))
b_NC <- length(setdiff(ncos_subset$species, copr_subset$species))
c_NC <- length(setdiff(copr_subset$species, ncos_subset$species))

Sorenson_NC <- ((2*a_NC)/((2*a_NC)+b_NC+c_NC))
Sorenson_NC

Jaccard_NC <- (a_NC)/(a_NC+b_NC+c_NC)

# SECOND, Lagoon and COPR
# a_LC is common to both Lagoon and COPR
# b_LC is taxa exclusive to Lagoon
# c_LC is taxa exclusive to COPR

a_LC <- length(intersect(lago_subset$species, copr_subset$species))
b_LC <- length(setdiff(lago_subset$species, copr_subset$species))
c_LC <- length(setdiff(copr_subset$species, lago_subset$species))

Sorenson_LC <- ((2*a_LC)/((2*a_LC)+b_LC+c_LC))
Sorenson_LC

# THIRD, Lagoon and NCOS
# a_LN is common to both Lagoon and NCOS
# b_LN is taxa exclusive to Lagoon
# c_LN is taxa exclusive to NCOS

a_LN <- length(intersect(lago_subset$species, ncos_subset$species))
b_LN <- length(setdiff(lago_subset$species, ncos_subset$species))
c_LN <- length(setdiff(ncos_subset$species, lago_subset$species))

Sorenson_LN <- ((2*a_LN)/((2*a_LN)+b_LN+c_LN))
Sorenson_LN

Sorenson_NC #NCOS and COPR
Sorenson_LC #Lagoon and COPR
Sorenson_LN #Lagoon and NCOS

Sorenson_mean <- (Sorenson_NC+Sorenson_LC+Sorenson_LN)/3

Sorenson_matrix <- data.frame("Sites" = c("North Campus Open Space vs. UCSB Lagoon",
                                                "North Campus Open Space vs. Coal Oil Point Reserve",
                                                "Coal Oil Point Reserve vs. UCSB Lagoon"),
                             "Sorenson_value" = c(Sorenson_LN,
                                                   Sorenson_NC,
                                                   Sorenson_LC))
Jaccard_NC <- (a_NC)/(a_NC+b_NC+c_NC)
Jaccard_LC <- (a_LC)/(a_LC+b_LC+c_LC)
Jaccard_LN <- (a_LN)/(a_LN+b_LN+c_LN)


###  Measure Sorenson similarity by site for PLANTS
plant_euler$species <- plant_euler$scientificName
ncos_subset <- subset(plant_euler, `North Campus Open Space` > 0)
lago_subset <- subset(plant_euler, `UCSB Lagoon` > 0)
copr_subset <- subset(plant_euler, `Coal Oil Point Reserve` > 0)

a_NC <- length(intersect(ncos_subset$species, copr_subset$species))
b_NC <- length(setdiff(ncos_subset$species, copr_subset$species))
c_NC <- length(setdiff(copr_subset$species, ncos_subset$species))

Sorenson_NC_plant <- ((2*a_NC)/((2*a_NC)+b_NC+c_NC))
Sorenson_NC_plant

a_LC <- length(intersect(lago_subset$species, copr_subset$species))
b_LC <- length(setdiff(lago_subset$species, copr_subset$species))
c_LC <- length(setdiff(copr_subset$species, lago_subset$species))

Sorenson_LC_plant <- ((2*a_LC)/((2*a_LC)+b_LC+c_LC))
Sorenson_LC_plant

a_LN <- length(intersect(lago_subset$species, ncos_subset$species))
b_LN <- length(setdiff(lago_subset$species, ncos_subset$species))
c_LN <- length(setdiff(ncos_subset$species, lago_subset$species))

Sorenson_LN_plant <- ((2*a_LN)/((2*a_LN)+b_LN+c_LN))
Sorenson_LN_plant

Sorenson_NC_plant #NCOS and COPR
Sorenson_LC_plant #Lagoon and COPR
Sorenson_LN_plant #Lagoon and NCOS

Sorenson_mean_plant <- (Sorenson_NC_plant+Sorenson_LC_plant+Sorenson_LN_plant)/3

Sorenson_matrix_plant <- data.frame("Sites" = c("North Campus Open Space vs. UCSB Lagoon",
                                          "North Campus Open Space vs. Coal Oil Point Reserve",
                                          "Coal Oil Point Reserve vs. UCSB Lagoon"),
                              "Sorenson_value" = c(Sorenson_LN_plant,
                                                   Sorenson_NC_plant,
                                                   Sorenson_LC_plant))
Sorenson_combined <- data.frame("Sites" = c("North Campus Open Space vs. UCSB Lagoon",
                                            "North Campus Open Space vs. Coal Oil Point Reserve",
                                            "Coal Oil Point Reserve vs. UCSB Lagoon"),
                                "Bees_Sorenson_value" = c(Sorenson_LN,
                                                     Sorenson_NC,
                                                     Sorenson_LC),
                                "Plants_Sorenson_value" = c(Sorenson_LN_plant,
                                                     Sorenson_NC_plant,
                                                     Sorenson_LC_plant))
Jaccard_NC_plant <- (a_NC)/(a_NC+b_NC+c_NC)
Jaccard_LC_plant <- (a_LC)/(a_LC+b_LC+c_LC)
Jaccard_LN_plant <- (a_LN)/(a_LN+b_LN+c_LN)

Jaccard_combined <- data.frame("Sites" = c("North Campus Open Space vs. UCSB Lagoon",
                                            "North Campus Open Space vs. Coal Oil Point Reserve",
                                            "Coal Oil Point Reserve vs. UCSB Lagoon"),
                                "Bees_Jaccard_value" = c(Jaccard_LN,
                                                          Jaccard_NC,
                                                          Jaccard_LC),
                                "Plants_Jaccard_value" = c(Jaccard_LN_plant,
                                                            Jaccard_NC_plant,
                                                            Jaccard_LC_plant))

## 2b. Bee community differences among sites

### 2b.i. Rarefaction-based ecological null hypothesis (Cayuela et al 2015)
#### Input: bee species matrix
sp_mat_bee_df <- sp_mat_bee[,c(4:57)]
#### Output: values
# # Ecological
set.seed(900)
ecoq3sites_m<-EcoTest.sample(sp_mat_bee_df, by=sp_mat_bee$site, MARGIN = 1, 
                             niter=200, q=0, method = "sample")
plot(ecoq3sites_m) #P(Obs <= null) = 0.31, Z = 0.29
# # Biogeographical -- won't run the biogeographical unless we reject the ecological
#biogq3sites_m<-BiogTest.sample(sp_mat_bee_df, by=sp_mat_bee$site, MARGIN = 1, 
#                         niter=200, q=0, method = "coverage")
#plot(biogq3sites_m) #P(Obs <= null) = 0.255



### 2b.ii. Chao1 diversity differences with ANOVA
#### Input: sp matrix bee
#### Output: values

col_counts <- tibble(
  fieldNumber = c("NCOS_01","NCOS_02","NCOS_03","NCOS_04","NCOS_05",
                  "NCOS_06","NCOS_07a","NCOS_07b",
                  "NCOS_A","NCOS_B","NCOS_C","NCOS_D","NCOS_E"),
  ncols = c(14,12,13,11,9,11,10,13,2,2,2,2,2)
)

chao_ncos_summary <- sp_mat_bee %>%
  inner_join(col_counts, by = "fieldNumber") %>%
  group_split(fieldNumber, ncols) %>%
  purrr::map_dfr(function(df) {
    
    subsite <- unique(df$fieldNumber)
    ncols   <- unique(df$ncols)
    
    df_t <- df %>%
      select(-ncols) %>%         # remove helper column
      t() %>%
      as.data.frame() %>%
      slice(-(1:3))              # remove first 3 rows
    
    df_t[, 1:ncols] <- lapply(df_t[, 1:ncols], function(x)
      as.numeric(as.character(x)))
    
    df_t$count <- rowSums(df_t[, 1:ncols], na.rm = TRUE)
    
    out <- chao1(df_t$count)
    out$subsite <- subsite
    out$site <- "NCOS"
    
    out
  })
# now, COPR
col_counts <- tibble(
  fieldNumber = c("COPR_01","COPR_02","COPR_03","COPR_04","COPR_05",
                  "COPR_06","COPR_A"),
  ncols = c(14,7,14,14,15,12,2)
)

chao_copr_summary <- sp_mat_bee %>%
  inner_join(col_counts, by = "fieldNumber") %>%
  group_split(fieldNumber, ncols) %>%
  purrr::map_dfr(function(df) {
    
    subsite <- unique(df$fieldNumber)
    ncols   <- unique(df$ncols)
    
    df_t <- df %>%
      select(-ncols) %>%         # remove helper column
      t() %>%
      as.data.frame() %>%
      slice(-(1:3))              # remove first 3 rows
    
    df_t[, 1:ncols] <- lapply(df_t[, 1:ncols], function(x)
      as.numeric(as.character(x)))
    
    df_t$count <- rowSums(df_t[, 1:ncols], na.rm = TRUE)
    
    out <- chao1(df_t$count)
    out$subsite <- subsite
    out$site <- "COPR"
    
    out
  })

# now, Lagoon
col_counts <- tibble(
  fieldNumber = c("Lagoon_01","Lagoon_02","Lagoon_03","Lagoon_04","Lagoon_05",
                  "Lagoon_A","Lagoon_C","Lagoon_D"),
  ncols = c(10,14,15,16,15,2,2,2)
)

chao_lago_summary <- sp_mat_bee %>%
  inner_join(col_counts, by = "fieldNumber") %>%
  group_split(fieldNumber, ncols) %>%
  purrr::map_dfr(function(df) {
    
    subsite <- unique(df$fieldNumber)
    ncols   <- unique(df$ncols)
    
    df_t <- df %>%
      select(-ncols) %>%         # remove helper column
      t() %>%
      as.data.frame() %>%
      slice(-(1:3))              # remove first 3 rows
    
    df_t[, 1:ncols] <- lapply(df_t[, 1:ncols], function(x)
      as.numeric(as.character(x)))
    
    df_t$count <- rowSums(df_t[, 1:ncols], na.rm = TRUE)
    
    out <- chao1(df_t$count)
    out$subsite <- subsite
    out$site <- "LAGO"
    
    out
  })
# combine
chao_overall_summary <- rbind(chao_ncos_summary,
                              chao_lago_summary,
                              chao_copr_summary)

# ANOVA on chao diversity metrics for bees
anova <- aov(S.chao1 ~ site, data = chao_overall_summary)
summary(anova) # F(2,25) = 1.22, p = 0.312


### 2b.iii. PERMANOVA of bees by site
#### Input: bee species matrix
#### Output: values
set.seed(731)
bee.dist<- vegdist(sp_mat_bee_df, method = "bray") 
set.seed(731)
PERM_bees <-adonis2(bee.dist~site,data=sp_mat_bee,permutations=999,method = "bray")
summary(PERM_bees) # p = 0.001
PERM_bees # F(2,259) = 3.89, R2 = 0.03, p = 0.001

set.seed(731)
bee.dist<- vegdist(sp_mat_bee_df, method = "jaccard") 
set.seed(731)
PERM_bees <-adonis2(bee.dist~site,data=sp_mat_bee,permutations=999,method = "jaccard")
summary(PERM_bees) # p = 0.001
PERM_bees # F(2,259) = 2.95, R2 = 0.02, p = 0.001

# permdisp bees
groups <- as.factor(sp_mat_bee$site)
permdisp_bees <- betadisper(bee.dist, groups)
permdisp_bees_test <- permutest(permdisp_bees, permutations = 999)
permdisp_bees_test

# proportions and bray curtis
comm_prop <- sp_mat_bee_df / rowSums(sp_mat_bee_df)
bee.dist<- vegdist(comm_prop, method = "bray") 
set.seed(731)
PERM_bees <-adonis2(bee.dist~site,data=sp_mat_bee,permutations=999,method = "bray")
summary(PERM_bees) # p = 0.001
PERM_bees # F(2,259) = 4.71, R2 = 0.04, p = 0.001


## 2c. Plant community differences among sites
sp_mat_plant_site <- sp_mat_plant %>% 
  mutate(site = ifelse(row_number() <= 33, 
                       "North Campus Open Space",
                       "UCSB Lagoon")) %>% 
  relocate(site,.after = year)
### 2c.i. Rarefaction-based ecological null hypothesis (Cayuela et al 2015)
ecoqPlant<-EcoTest.sample(sp_mat_plant, by=sp_mat_plant_site$site, 
                          MARGIN = 1, niter=200, q=0, method = "coverage")
plot(ecoqPlant) #P(Obs <= null) = 0.005; Zeco = 32.97
ecoqPlant
### 2c.ii. Rarefaction-based biogeographical null hypothesis (Cayuela et al 2015)

biogqPlant<-BiogTest.sample(sp_mat_plant, by=sp_mat_plant_site$site, 
                            MARGIN = 1, niter=150, q=0, method = "coverage")
plot(biogqPlant) #P(Obs <= null) = 0.41; Zbio = 1.49
biogqPlant

### 2c.iii. PERMANOVA of plants by site
sp_plant_mat_df <- sp_mat_plant[,c(3:174)]
set.seed(731)
plant.dist<- vegdist(sp_plant_mat_df, method = "bray")
set.seed(731)
PERM_plants<-adonis2(plant.dist~site,data=sp_mat_plant_site,permutations=999)
summary(PERM_plants) # p = 0.001
PERM_plants # F(1,51) = 18.77, R2 = 0.27, p = 0.001
set.seed(731)
plant.nmds <- metaMDS(sp_plant_mat_df, 
                      k=2, model = "global", distance = "bray")
plant.nmds$stress #0.182
colvec_p<-c("#984EA3", "orange")
plot(plant.nmds)
orditorp(plant.nmds, display = 'species')
with(sp_mat_plant_site,
     points(plant.nmds, pch=15,cex=1.0,
            col = colvec_p[site]))
mtext(text = "Stress = 0.182", side = 3, line = -1.5, adj = 0.95, 
      cex = 1.0)
legend("bottomright", c("North Campus Open Space", "UCSB Lagoon"), title = "Site", 
       col= colvec_p, pch = 15, cex=1)
ordiellipse(plant.nmds, groups = sp_mat_plant_site$site, draw = "polygon", 
            label = T, lty = 1, col = colvec_p)

# 3. Question 2: How important are native plants in determining bee community composition?

## 3a. NMDS
### Inputs: Combined species matrix (bees and plants)
### Outputs: supplemental figs
sp_mat_all_b <- sp_mat_all[,c(5:46)]
sp_mat_all_p <- sp_mat_all[,c(47:200)]

set.seed(413)
bee.nmds <- metaMDS(sp_mat_all_b, 
                    k=2, model = "global", distance = "bray")
ef_bee<-envfit(bee.nmds, sp_mat_all_p, permu=999, na.rm= TRUE )


# Site scores
sites <- as.data.frame(scores(bee.nmds, display = "sites"))
sites$site_group <- sp_mat_all$site

# Species scores
species <- as.data.frame(scores(bee.nmds, display = "species"))
species$species <- rownames(species)
vec <- as.data.frame(scores(ef_bee, display = "vectors"))
vec$var <- rownames(vec)

# Keep only significant vectors
vec$pval <- ef_bee$vectors$pvals
vec <- vec %>% filter(pval <= 0.05)

# Scale arrows properly
arrow_mult <- ordiArrowMul(ef_bee)
vec$NMDS1 <- vec$NMDS1 * arrow_mult
vec$NMDS2 <- vec$NMDS2 * arrow_mult
ggplot() +
  
  # Site points
  geom_point(data = sites,
             aes(x = NMDS1, y = NMDS2, color = site_group),
             size = 3) +
  
  # Ellipses (95% SE by default)
  stat_ellipse(data = sites,
               aes(x = NMDS1, y = NMDS2, fill = site_group),
               geom = "polygon",
               alpha = 0.2,
               color = NA) +
  
  # Species labels
  geom_text(data = species,
            aes(x = NMDS1, y = NMDS2, label = species),
            color = "grey40",
            size = 3) +
  
  # Envfit arrows
  geom_segment(data = vec,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "red") +
  
  # Envfit labels
  geom_text(data = vec,
            aes(x = NMDS1, y = NMDS2, label = var),
            color = "red",
            vjust = -0.5) +
  
  theme_minimal() +
  coord_equal() +
  labs(color = "Site", fill = "Site")

#######

list(ef_bee, p.max = 0.05)
#### 3a.i. NMDS Year 2018
sp_mat_all_y_2018 <- sp_mat_all %>% filter(year == 2018)
sp_mat_all_b_2018 <- sp_mat_all %>% filter(year == 2018) %>% select(5:46)
sp_mat_all_p_2018 <- sp_mat_all %>% filter(year == 2018) %>% select(47:200)
bee.nmds_2018 <- metaMDS(sp_mat_all_b_2018, k=2, model = "global", distance = "bray")
plot(bee.nmds_2018)
orditorp(bee.nmds_2018, display = 'species')
ordiellipse(bee.nmds_2018, groups = sp_mat_all_y_2018$site, draw = "polygon", 
            label = TRUE, lty = 1)
ef_bee_2018 <- envfit(bee.nmds_2018, sp_mat_all_p_2018, permu=999, na.rm= TRUE )
list(ef_bee_2018, p.max = 0.05)

#### 3a.ii. NMDS Year 2019
sp_mat_all_y_2019 <- sp_mat_all %>% filter(year == 2019)
sp_mat_all_b_2019 <- sp_mat_all %>% filter(year == 2019) %>% select(5:46)
sp_mat_all_p_2019 <- sp_mat_all %>% filter(year == 2019) %>% select(47:200)
bee.nmds_2019 <- metaMDS(sp_mat_all_b_2019, k=2, model = "global", distance = "bray")
plot(bee.nmds_2019)
orditorp(bee.nmds_2019, display = 'species')
ordiellipse(bee.nmds_2019, groups = sp_mat_all_y_2019$site, draw = "polygon", 
            label = TRUE, lty = 1)
ef_bee_2019 <- envfit(bee.nmds_2019, sp_mat_all_p_2019, permu=999, na.rm= TRUE )
list(ef_bee_2019, p.max = 0.05)

#### 3a.iii. NMDS Year 2020
sp_mat_all_y_2020 <- sp_mat_all %>% filter(year == 2020)
sp_mat_all_b_2020 <- sp_mat_all %>% filter(year == 2020) %>% select(5:46)
sp_mat_all_p_2020 <- sp_mat_all %>% filter(year == 2020) %>% select(47:200)
bee.nmds_2020 <- metaMDS(sp_mat_all_b_2020, k=2, model = "global", distance = "bray")
plot(bee.nmds_2020)
orditorp(bee.nmds_2020, display = 'species')
ordiellipse(bee.nmds_2020, groups = sp_mat_all_y_2020$site, draw = "polygon", 
            label = TRUE, lty = 1)
ef_bee_2020 <- envfit(bee.nmds_2020, sp_mat_all_p_2020, permu=999, na.rm= TRUE )
list(ef_bee_2020, p.max = 0.05)

#### 3a.iv. NMDS Year 2023
sp_mat_all_y_2023 <- sp_mat_all %>% filter(year == 2023)
sp_mat_all_b_2023 <- sp_mat_all %>% filter(year == 2023) %>% select(5:46)
sp_mat_all_p_2023 <- sp_mat_all %>% filter(year == 2023) %>% select(47:200)
bee.nmds_2023 <- metaMDS(sp_mat_all_b_2023, k=2, model = "global", distance = "bray")
plot(bee.nmds_2023)
orditorp(bee.nmds_2023, display = 'species')
ordiellipse(bee.nmds_2023, groups = sp_mat_all_y_2023$site, draw = "polygon", 
            label = TRUE, lty = 1)
ef_bee_2023 <- envfit(bee.nmds_2023, sp_mat_all_p_2023, permu=999, na.rm= TRUE )
list(ef_bee_2023, p.max = 0.05)

## 3b. NMDS envfit analysis
### Inputs: NMDS from each year
### Outputs: list of plants significant in each year
ef_bee_2018
ef_bee_2019
ef_bee_2020
ef_bee_2023


### Plot all 4 NMDS with each year separate

library(ggrepel)
###############
# 2018 NMDS ###
###############
# Site scores
sites <- as.data.frame(scores(bee.nmds_2018, display = "sites"))
sites$site_group <- sp_mat_all_y_2018$site

# Species scores
species <- as.data.frame(scores(bee.nmds_2018, display = "species"))
species$species <- rownames(species)
vec <- as.data.frame(scores(ef_bee_2018, display = "vectors"))
vec$var <- rownames(vec)

# Keep only significant vectors
vec$pval <- ef_bee_2018$vectors$pvals
vec <- vec %>% filter(pval <= 0.05)

# Scale arrows properly
arrow_mult <- ordiArrowMul(ef_bee_2018)
vec$NMDS1 <- vec$NMDS1 * arrow_mult
vec$NMDS2 <- vec$NMDS2 * arrow_mult
plot2018 <- ggplot() +
  
  # Site points
  geom_point(data = sites,
             aes(x = NMDS1, y = NMDS2, color = site_group),
             size = 3) +
  
  # Ellipses (95% SE by default)
  stat_ellipse(data = sites,
               aes(x = NMDS1, y = NMDS2, fill = site_group),
               geom = "polygon",
               alpha = 0.2,
               color = NA) +
  
  # Species labels
  geom_text_repel(data = species,
            aes(x = NMDS1, y = NMDS2, label = species),
            color = "grey40",
            size = 3) +
  
  # Envfit arrows
  geom_segment(data = vec,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "red") +
  
  # Envfit labels
  geom_text_repel(data = vec,
            aes(x = NMDS1, y = NMDS2, label = var),
            color = "red",
            vjust = -0.5) +
  
  theme_minimal() +
  coord_equal() +
  labs(color = "Site", fill = "Site",
       title = "(A) 2018")


###############
# 2019 NMDS ###
###############
# Site scores
sites <- as.data.frame(scores(bee.nmds_2019, display = "sites"))
sites$site_group <- sp_mat_all_y_2019$site

# Species scores
species <- as.data.frame(scores(bee.nmds_2019, display = "species"))
species$species <- rownames(species)
vec <- as.data.frame(scores(ef_bee_2019, display = "vectors"))
vec$var <- rownames(vec)

# Keep only significant vectors
vec$pval <- ef_bee_2019$vectors$pvals
vec <- vec %>% filter(pval <= 0.05)

# Scale arrows properly
arrow_mult <- ordiArrowMul(ef_bee_2019)
vec$NMDS1 <- vec$NMDS1 * arrow_mult
vec$NMDS2 <- vec$NMDS2 * arrow_mult
plot2019 <- ggplot() +
  
  # Site points
  geom_point(data = sites,
             aes(x = NMDS1, y = NMDS2, color = site_group),
             size = 3) +
  
  # Ellipses (95% SE by default)
  stat_ellipse(data = sites,
               aes(x = NMDS1, y = NMDS2, fill = site_group),
               geom = "polygon",
               alpha = 0.2,
               color = NA) +
  
  # Species labels
  geom_text_repel(data = species,
            aes(x = NMDS1, y = NMDS2, label = species),
            color = "grey40",
            size = 3) +
  
  # Envfit arrows
  geom_segment(data = vec,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "red") +
  
  # Envfit labels
  geom_text_repel(data = vec,
            aes(x = NMDS1, y = NMDS2, label = var),
            color = "red",
            vjust = -0.5) +
  
  theme_minimal() +
  coord_equal() +
  labs(color = "Site", fill = "Site",
       title = "(B) 2019")

###############
# 2020 NMDS ###
###############
# Site scores
sites <- as.data.frame(scores(bee.nmds_2020, display = "sites"))
sites$site_group <- sp_mat_all_y_2020$site

# Species scores
species <- as.data.frame(scores(bee.nmds_2020, display = "species"))
species$species <- rownames(species)
vec <- as.data.frame(scores(ef_bee_2020, display = "vectors"))
vec$var <- rownames(vec)

# Keep only significant vectors
vec$pval <- ef_bee_2020$vectors$pvals
vec <- vec %>% filter(pval <= 0.05)

# Scale arrows properly
arrow_mult <- ordiArrowMul(ef_bee_2020)
vec$NMDS1 <- vec$NMDS1 * arrow_mult
vec$NMDS2 <- vec$NMDS2 * arrow_mult
plot2020 <- ggplot() +
  
  # Site points
  geom_point(data = sites,
             aes(x = NMDS1, y = NMDS2, color = site_group),
             size = 3) +
  
  # Ellipses (95% SE by default)
  stat_ellipse(data = sites,
               aes(x = NMDS1, y = NMDS2, fill = site_group),
               geom = "polygon",
               alpha = 0.2,
               color = NA) +
  
  # Species labels
  geom_text_repel(data = species,
            aes(x = NMDS1, y = NMDS2, label = species),
            color = "grey40",
            size = 3) +
  
  # Envfit arrows
  geom_segment(data = vec,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "red") +
  
  # Envfit labels
  geom_text_repel(data = vec,
            aes(x = NMDS1, y = NMDS2, label = var),
            color = "red",
            vjust = -0.5) +
  
  theme_minimal() +
  coord_equal() +
  labs(color = "Site", fill = "Site",
       title = "(C) 2020")

###############
# 2023 NMDS ###
###############
# Site scores
sites <- as.data.frame(scores(bee.nmds_2023, display = "sites"))
sites$site_group <- sp_mat_all_y_2023$site

# Species scores
species <- as.data.frame(scores(bee.nmds_2023, display = "species"))
species$species <- rownames(species)
vec <- as.data.frame(scores(ef_bee_2023, display = "vectors"))
vec$var <- rownames(vec)

# Keep only significant vectors
vec$pval <- ef_bee_2023$vectors$pvals
vec <- vec %>% filter(pval <= 0.05)

# Scale arrows properly
arrow_mult <- ordiArrowMul(ef_bee_2023)
vec$NMDS1 <- vec$NMDS1 * arrow_mult
vec$NMDS2 <- vec$NMDS2 * arrow_mult
plot2023 <- ggplot() +
  
  # Site points
  geom_point(data = sites,
             aes(x = NMDS1, y = NMDS2, color = site_group),
             size = 3) +
  
  # Ellipses (95% SE by default)
  stat_ellipse(data = sites,
               aes(x = NMDS1, y = NMDS2, fill = site_group),
               geom = "polygon",
               alpha = 0.2,
               color = NA) +
  
  # Species labels
  geom_text_repel(data = species,
            aes(x = NMDS1, y = NMDS2, label = species),
            color = "grey40",
            size = 3) +
  
  # Envfit arrows
  geom_segment(data = vec,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "red") +
  
  # Envfit labels
  geom_text_repel(data = vec,
            aes(x = NMDS1, y = NMDS2, label = var),
            color = "red",
            vjust = -0.5) +
  
  theme_minimal() +
  coord_equal() +
  labs(color = "Site", fill = "Site",
       title = "(D) 2023")

# plots
plot2018
plot2019
plot2020
plot2023

library(patchwork)
finalplot <- (plot2018 | plot2019) /
(plot2020 | plot2023) +
  plot_layout(guides = "collect")+
theme(legend.position = "bottom")

finalplot

ggsave(finalplot, 
       filename = "../figures/nmds_plots.jpeg",
       width=20,
       units="cm",
       dpi = 300)

### REPEAT WITHOUT THE RED ARROWS IN EACH PANEL
# Site scores
sites <- as.data.frame(scores(bee.nmds_2018, display = "sites"))
sites$site_group <- sp_mat_all_y_2018$site

# Species scores
species <- as.data.frame(scores(bee.nmds_2018, display = "species"))
species$species <- rownames(species)
vec <- as.data.frame(scores(ef_bee_2018, display = "vectors"))
vec$var <- rownames(vec)

# Keep only significant vectors
vec$pval <- ef_bee_2018$vectors$pvals
vec <- vec %>% filter(pval <= 0.05)

# Scale arrows properly
arrow_mult <- ordiArrowMul(ef_bee_2018)
vec$NMDS1 <- vec$NMDS1 * arrow_mult
vec$NMDS2 <- vec$NMDS2 * arrow_mult
plot2018 <- ggplot() +
  geom_point(data = sites,
             aes(x = NMDS1, y = NMDS2, color = site_group),
             size = 3) +
  stat_ellipse(data = sites,
               aes(x = NMDS1, y = NMDS2, fill = site_group),
               geom = "polygon",
               alpha = 0.2,
               color = NA) +
  geom_text_repel(data = species,
            aes(x = NMDS1, y = NMDS2, label = species),
            color = "grey40",
            size = 3) +
  theme_minimal() +
  coord_equal() +
  labs(color = "Site", fill = "Site",
       title = "(A) 2018")


###############
# 2019 NMDS ###
###############
# Site scores
sites <- as.data.frame(scores(bee.nmds_2019, display = "sites"))
sites$site_group <- sp_mat_all_y_2019$site

# Species scores
species <- as.data.frame(scores(bee.nmds_2019, display = "species"))
species$species <- rownames(species)
vec <- as.data.frame(scores(ef_bee_2019, display = "vectors"))
vec$var <- rownames(vec)

# Keep only significant vectors
vec$pval <- ef_bee_2019$vectors$pvals
vec <- vec %>% filter(pval <= 0.05)

# Scale arrows properly
arrow_mult <- ordiArrowMul(ef_bee_2019)
vec$NMDS1 <- vec$NMDS1 * arrow_mult
vec$NMDS2 <- vec$NMDS2 * arrow_mult
plot2019 <- ggplot() +
  geom_point(data = sites,
             aes(x = NMDS1, y = NMDS2, color = site_group),
             size = 3) +
  stat_ellipse(data = sites,
               aes(x = NMDS1, y = NMDS2, fill = site_group),
               geom = "polygon",
               alpha = 0.2,
               color = NA) +
  geom_text_repel(data = species,
            aes(x = NMDS1, y = NMDS2, label = species),
            color = "grey40",
            size = 3) +
  theme_minimal() +
  coord_equal() +
  labs(color = "Site", fill = "Site",
       title = "(B) 2019")

###############
# 2020 NMDS ###
###############
# Site scores
sites <- as.data.frame(scores(bee.nmds_2020, display = "sites"))
sites$site_group <- sp_mat_all_y_2020$site

# Species scores
species <- as.data.frame(scores(bee.nmds_2020, display = "species"))
species$species <- rownames(species)
vec <- as.data.frame(scores(ef_bee_2020, display = "vectors"))
vec$var <- rownames(vec)

# Keep only significant vectors
vec$pval <- ef_bee_2020$vectors$pvals
vec <- vec %>% filter(pval <= 0.05)

# Scale arrows properly
arrow_mult <- ordiArrowMul(ef_bee_2020)
vec$NMDS1 <- vec$NMDS1 * arrow_mult
vec$NMDS2 <- vec$NMDS2 * arrow_mult
plot2020 <- ggplot() +
  geom_point(data = sites,
             aes(x = NMDS1, y = NMDS2, color = site_group),
             size = 3) +
  stat_ellipse(data = sites,
               aes(x = NMDS1, y = NMDS2, fill = site_group),
               geom = "polygon",
               alpha = 0.2,
               color = NA) +
  geom_text_repel(data = species,
            aes(x = NMDS1, y = NMDS2, label = species),
            color = "grey40",
            size = 3) +
  theme_minimal() +
  coord_equal() +
  labs(color = "Site", fill = "Site",
       title = "(C) 2020")

###############
# 2023 NMDS ###
###############
# Site scores
sites <- as.data.frame(scores(bee.nmds_2023, display = "sites"))
sites$site_group <- sp_mat_all_y_2023$site

# Species scores
species <- as.data.frame(scores(bee.nmds_2023, display = "species"))
species$species <- rownames(species)
vec <- as.data.frame(scores(ef_bee_2023, display = "vectors"))
vec$var <- rownames(vec)

# Keep only significant vectors
vec$pval <- ef_bee_2023$vectors$pvals
vec <- vec %>% filter(pval <= 0.05)

# Scale arrows properly
arrow_mult <- ordiArrowMul(ef_bee_2023)
vec$NMDS1 <- vec$NMDS1 * arrow_mult
vec$NMDS2 <- vec$NMDS2 * arrow_mult
plot2023 <- ggplot() +
  geom_point(data = sites,
             aes(x = NMDS1, y = NMDS2, color = site_group),
             size = 3) +
  stat_ellipse(data = sites,
               aes(x = NMDS1, y = NMDS2, fill = site_group),
               geom = "polygon",
               alpha = 0.2,
               color = NA) +
  geom_text_repel(data = species,
            aes(x = NMDS1, y = NMDS2, label = species),
            color = "grey40",
            size = 3) +
  theme_minimal() +
  coord_equal() +
  labs(color = "Site", fill = "Site",
       title = "(D) 2023")

# plots
plot2018
plot2019
plot2020
plot2023

library(patchwork)
finalplot <- (plot2018 | plot2019) /
(plot2020 | plot2023) +
  plot_layout(guides = "collect")+
theme(legend.position = "bottom")

finalplot

ggsave(finalplot, 
       filename = "../figures/nmds_plots_noVec.jpeg",
       width=20,
       units="cm",
       dpi = 300)


# 4. Supplemental Question: Estimate other bee species in this region
## Input: occurrence data from nearby coastal areas; 
### broader occurrence data to filter to pan traps
## Output: expected species list for this region with pan trap collections recorded
### Make maps, convex hulls, extract those that overlap with our sites. 

library(sf)
df <- read.csv("../data/occ_beelib_16mar2026.csv")
target <- df
# make spatial object
target_sp_all <- target[,c("scientificName",
                           "decimalLatitude",
                           "decimalLongitude")]
target.sf <- target_sp_all %>% 
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

target_hulls <- target.sf %>% 
  group_by(scientificName) %>% 
  summarise(geometry = st_combine(geometry)) %>% 
  st_convex_hull()

polygons <- target_hulls[st_geometry_type(target_hulls)=="POLYGON",]



region_sf <- st_as_sfc(st_bbox(c(xmin = -119.977989,
                                 xmax = -119.515938,
                                 ymin = 34.3857640,
                                 ymax = 34.460040),
                               crs = 4326))
st_crs(polygons)==st_crs(region_sf)
intersections <- st_intersects(polygons,region_sf)
hits <- lengths(intersections) > 0
species_in_region <- unique(polygons$scientificName[hits])
species_in_region
species_in_region_df <- as.data.frame(species_in_region)
species_in_region_df <- species_in_region_df %>% rename(scientificName = species_in_region)

# broader occurrence data to filter down to pan trappable species
df1 <- read.csv("../data/occ_beelib_ca_16mar2026.csv")
df1 <- df1 %>% filter(grepl("pan",samplingProtocol,ignore.case = TRUE))
df2 <- df1 %>% distinct(scientificName)
df2$pantrap <- "pan"

# combine the pantrappable list with hull list for region
df3 <- left_join(species_in_region_df,df2,by = "scientificName")
df4 <- subset(df3,pantrap == "pan")
df5 <- df4[grepl("\\s", df4$scientificName), ] # drop taxa IDed to genus only
df6 <- df5[!grepl("Bombus", df5$scientificName), ]
df6 <- subset(df6, scientificName != "Agapostemon subtilior") # species renamed from texanus in this region, but we are using texanus still in this project
df6 <- subset(df6, scientificName != "Colletes gaudialis")  # drop duplicate, retain subsp.
df7 <- subset(df6, !scientificName %in% c("Lasioglossum (Dialictus)",
                                          "Lasioglossum (Lasioglossum)",
                                          "Lasioglossum (Evylaeus)"))
