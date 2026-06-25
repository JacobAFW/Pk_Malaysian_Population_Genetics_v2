# packages
library(LEA)
library(tess3r)
library(tidyverse)
library(maps)

# Import data
genotype <- ped2lfmm('clusters_combined_fixed.ped')
tes_genotype <- read.table(genotype)

unique(as.vector(as.matrix(tes_genotype[-1]))) # need to get rid of the 9

## replace missing calls with NA 
tes_genotype_clean <- tes_genotype %>%
     mutate_all(~ifelse(. == 9, NA, .))

original <- read.table("clusters_combined.fam")
gis <- read.table("complete_gis_na_omit.fam") 
gis <- gis[order(match(gis$V1, original$V1)),]
gis <- gis %>%
     select(V7, V8) %>% 
     mutate(long = ifelse(V7 < V8, V8, V7)) %>%
     mutate(lat = ifelse(V8 > V7, V7, V8)) %>%
     select(long, lat) %>%
     as.matrix()

nrow(tes_genotype_clean) == nrow(gis)

# Visualise distribution of GIS coords - quick/dirty version
png("gis_coordinates.png")
plot(gis, pch = 19, cex = .5, 
     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
maps::map(add = T, interior = F)
dev.off()

# Estimating acestry coefficients
tess3.obj <- tess3(X = tes_genotype_clean, coord = gis, K = 1:10, 
                   method = "projected.ls", ploidy = 2, openMP.core.num = 10) 


# Visualise distribution of root mean-squared errors computed on subset of loci for cross validation - quick/dirty version
png("K_cross-validation_clean.png")
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")
dev.off()

# Ancestry matrix barplot

q.matrix <- qmatrix(tess3.obj, K = 3)

# STRUCTURE-like barplot for the Q-matrix 
png("barplot_clean2_k_4.png")
barplot(q.matrix, border = NA, space = 0, 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "Ancestry matrix") -> bp
dev.off()

# STRUCTURE-like map for the Q-matrix 
png("barplot_clean2_k_4.png")
plot(q.matrix, gis, method = "map.max", interpol = FieldsKrigModel(10),  
     main = "Ancestry coefficients",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4)
dev.off()









####################################### WIP - IMPUTATIONS
# Other matrix options


## imputation to replace missing calls


### Imputation 1

### estimates admixture coefficient
project.snmf = snmf(genotype, K = 3, 
        entropy = TRUE, 
        repetitions = 10,
        project = "new",
        CPU = 10)

### select the run with the lowest cross-entropy value for a given K
best = which.min(cross.entropy(project.snmf, K = 3))

### impute the missing genotypes
impute(project.snmf, genotype, method = 'mode', K = 3, run = best)

tes_genotype_imputed <- read.table("clusters_combined_fixed.lfmm_imputed.lfmm")


####################################

### Imputation 2

project.snmf_2 = snmf(genotype, K = 1:8, 
        entropy = TRUE, 
        repetitions = 10,
        project = "new",
        CPU = 10)

# plot cross-entropy criterion of all runs of the project
png("cross-entropy_for_imputation.png.png")
plot(project.snmf_2, lwd = 5, col = "red", pch=1)
dev.off()

### select the run with the lowest cross-entropy value for a given K - plateaus at 6
best = which.min(cross.entropy(project.snmf_2, K = 6))

### impute the missing genotypes
impute(project.snmf_2, genotype, method = 'mode', K = 6, run = best)

tes_genotype_imputed_2 <- read.table("clusters_combined_fixed.lfmm_imputed.lfmm")






#################################### WIP - NO MISS
## No missingess dataset

genotype_no_miss <- ped2lfmm('cleaned_fixed.ped')
tes_genotype_no_miss <- read.table(genotype_no_miss)

unique(as.vector(as.matrix(tes_genotype_no_miss[-1]))) # need to get rid of the 9

tes_genotype_no_miss <- tes_genotype_no_miss %>%
     mutate_all(~ifelse(. == 9, NA, .))

original <- read.table("cleaned.fam")
gis <- read.table("complete_gis_na_omit.fam") 
gis <- gis[order(match(gis$V1, original$V1)),]
gis <- gis %>%
     filter(V1 %in% original$V1) %>%
     select(V7, V8) %>% 
     mutate(long = ifelse(V7 < V8, V8, V7)) %>%
     mutate(lat = ifelse(V8 > V7, V7, V8)) %>%
     select(long, lat) %>%
     as.matrix()


nrow(tes_genotype_no_miss) == nrow(gis)

# Visualise distribution of GIS coords - quick/dirty version
png("gis_coordinates.png")
plot(gis, pch = 19, cex = .5, 
     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
maps::map(add = T, interior = F)
dev.off()

# Estimating acestry coefficients
tess3.obj <- tess3(X = tes_genotype_no_miss, coord = gis, K = 1:10, 
                   method = "projected.ls", ploidy = 2, openMP.core.num = 28) 


# Visualise distribution of root mean-squared errors computed on subset of loci for cross validation - quick/dirty version
png("K_cross-validation_clean.png")
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")
dev.off()



#################################### WIP











original <- read.table("clusters_combined.fam")
gis <- read.table("complete_gis_na_omit.fam") 
gis <- gis[order(match(gis$V1, original$V1)),]
gis <- gis %>%
     select(V7, V8) %>% 
     mutate(long = ifelse(V7 < V8, V8, V7)) %>%
     mutate(lat = ifelse(V8 > V7, V7, V8)) %>%
     select(long, lat) %>%
     as.matrix()

# Visualise distribution of GIS coords - quick/dirty version
png("gis_coordinates.png")
plot(gis, pch = 19, cex = .5, 
     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
maps::map(add = T, interior = F)
dev.off()

# Estimating acestry coefficients
tess3.obj <- tess3(X = tes_genotype_clean, coord = gis, K = 1:10, 
                   method = "projected.ls", ploidy = 2, openMP.core.num = 10) 


# Visualise distribution of root mean-squared errors computed on subset of loci for cross validation - quick/dirty version
png("K_cross-validation_imputed_2.png")
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")
dev.off()

png("K_cross-validation_clean_2.png")
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")
dev.off()




