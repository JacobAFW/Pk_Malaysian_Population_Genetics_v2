#Load packages
library(tidyverse)
library(ape)
library(ggtree)
library(viridis)

## Metadata 
metadata <- read_tsv("/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/data/metadata/full_metadata.tsv") %>% 
    mutate(sevemal = ifelse(grepl("UM", sevemal), "UM", .$sevemal)) %>% 
    mutate(sevemal = ifelse(sevemal == "1", "SM", ifelse(sevemal == "0", "UM", .$sevemal))) %>% 
    rename(Sample = studycode) %>% 
    mutate(district = str_to_title(district)) %>%
    add_column(Cluster = NA)

original_metadata <- readxl::read_xlsx("/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/data/metadata/PK_Sabah_Sample_naming_indexes.xlsx") %>%
    select(1:3) %>%
    mutate(subjectid = ifelse(str_length(subjectid) == 1, paste0("00", subjectid), # if subjectid length = 1 paste 00
        ifelse(str_length(subjectid) == 2, paste0("0", subjectid), # if not, then if subjectid length = 2 paste 0
        .$subjectid))) %>%
    unite(Sample, c("group", "subjectid"), sep ="") %>% 
    rename(Sample2 = sampleid)

metadata <- metadata %>%
    left_join(original_metadata) %>%
    mutate(Sample2 = ifelse(is.na(Sample2), Sample, Sample2)) %>%
    select(-Sample) %>%
    rename(Sample = Sample2) %>%
    relocate(Sample)
    
    
archived_data <- read_csv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/data/metadata/Pk_clusters_metadata.csv") %>%
        select(1, 3, 5) %>% 
        rename(district = area) %>%
        mutate(Cluster = str_remove(Group, "-Pk")) %>% 
        select(-Group) %>% 
        rbind(read_csv("/g/data/pq84/malaria/Pk_Malaysian_Population_Genetics/data/metadata/Pk_clusters_peninsular_metadata.csv") %>%
        select(2, 8) %>%
        rename(Sample = ENA_accession_no_ES) %>%
        rename(district = Location) %>%
        add_column(Cluster = "Peninsular")
        ) %>%
        add_column(enroldate=NA, age = NA, sex=NA, ethnicity=NA, village=NA, occupation=NA, Hb=NA, platelets=NA,
        parasitemia=NA, fever_days=NA, respiratory_rate=NA, oxygen_saturation=NA, systolic_BP=NA, diastolic_BP=NA, G6PD=NA,
        pregnant=NA, creatinine=NA, study=NA, bicarbonate=NA, PCR=NA, sevemal=NA, bilirubin=NA, glucose=NA, bleeding_severity=NA)


metadata <- metadata %>%
    rbind(archived_data) %>%
  mutate_at(c("district", "village"), str_to_title) %>%
  mutate(Sample2 = Sample) %>% 
  mutate(district = ifelse(district == "Kota Murudu", "Kota Marudu", .$district)) %>% 
  mutate(district = ifelse(district == "Sungai Siuit", "Sungai Siput", .$district))


# Bar plot of ADMIXTURE output

admixture_table <- read.table("admixture/cleaned.3.Q") %>%
    cbind(
        read.table("cleaned.fam") %>%
            select(1) %>%
            rename(Sample = V1)
    ) %>% 
    arrange(desc(V3), desc(V2), desc(V1)) %>% 
    add_column(Sample2 = as.factor(1:nrow(.))) %>%
    pivot_longer(1:3, 
        names_to = "Ancestry", 
        values_to = "Proportion") %>% 
  mutate(Cluster = ifelse(Ancestry == "V1", "Mn",
                          ifelse(Ancestry == "V2", "Mf", "Peninsular"))) 

admix_plot <- admixture_table %>% 
  mutate(Ancestry = Proportion) %>% 
  ggplot(aes(x = Sample2, y = Ancestry, fill = Cluster)) +
  geom_col() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("#440154FF", "#39568CFF", "#73D055FF")) +
  xlab("Sample")

ggsave("outputs/admix_bar_plot.png", dpi=300, admix_plot, width = 20)

# Clusters for population specific filters
admix_clusters <- admixture_table %>% 
    group_by(Sample) %>% 
    filter(Proportion == max(Proportion)) 

write_tsv(admix_clusters, "outputs/admix_clusters.tsv")

admix_clusters %>%
    filter(Cluster == "Mn" | Cluster == "Peninsular") %>%
    left_join(
        read.table("cleaned.fam") %>%
            mutate(Sample = V1) 
    ) %>%
    select(V1:V5) %>%
    write_tsv("Exclude_for_Mf.tsv")

admix_clusters %>%
    filter(Cluster == "Mf" | Cluster == "Peninsular") %>%
    left_join(
        read.table("cleaned.fam") %>%
            mutate(Sample = V1) 
    ) %>%
    select(V1:V5) %>%
    write_tsv("Exclude_for_Mn.tsv")

admix_clusters %>%
    filter(Cluster == "Mn" | Cluster == "Mf") %>%
    left_join(
        read.table("cleaned.fam") %>%
            mutate(Sample = V1) 
    ) %>%
    select(V1:V5) %>%
    write_tsv("Exclude_for_Pen.tsv")
    

# PCA

## Read in data
pca <- read_table("Pk.eigenvec", col_names=F) %>%
    select(-1) %>%
    rename("sampleid" = "X2") %>%
    rename_at(vars(starts_with("X")), ~str_replace(., "X.*", paste0("PC", seq_along(.))))

eigenval <- scan("Pk.eigenval")

pca <- pca %>%
    rename(Sample = sampleid) %>%
    mutate(Sample2 = str_remove(Sample, "_DK.*")) %>% # we need the original names downstream for the NJT
    left_join(
        metadata %>% mutate(Sample2 = Sample), 
        by = "Sample2"
        ) 


## Percentage variance explained by each PC
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

pve_plot <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100) %>%
    ggplot(aes(PC, pve)) + 
    geom_bar(stat = "identity") + 
    ylab("Percentage variance explained") + 
    theme_light()

ggsave("outputs/percentage_variance_explained.png", dpi=300, pve_plot)


## PCA - clusters
pca <- pca %>% 
    mutate(Cluster = ifelse(PC1 < 0 & PC2 < 0.04, "Mf", 
        ifelse(PC1 > 0.04 & PC2 < 0.04, "Mn", 
        ifelse(PC2 > 0.1, "Peninsular", 
        ifelse(PC1 > 0 & PC1 < 0.04 & PC2 < 0.04, "Mixed", .$Cluster))))) 

pca_plot <- ggplot(pca, aes(PC1, PC2, colour = Cluster)) + 
    geom_point(size = 3) +
    scale_color_manual(values = c("#440154FF", "#39568CFF", "#1F968BFF", "#73D055FF")) +
    coord_equal() + 
    theme_light() + 
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
    ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

ggsave("outputs/pca_cluster.png", dpi = 600, pca_plot)

## PCA - districts
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

pca_plot <- ggplot(pca, aes(PC1, PC2, colour = district)) + 
    geom_point(size = 3) +
    scale_color_viridis_d() +
    coord_equal() + 
    theme_light() + 
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
    ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

ggsave("outputs/pca_districts.png", dpi = 600, pca_plot)

## PCA - sevemal
pca_plot <- ggplot(pca, aes(PC1, PC2, colour = sevemal)) + 
    geom_jitter(height = 0.02, width = 0.02) +
    coord_equal() + 
    theme_light() + 
    scale_color_manual(values = c("#440154FF", "#1F968BFF")) +
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
    ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

ggsave("outputs/pca_sevemal.png", dpi = 600, pca_plot)

# Updated metadata with clusters
metadata <- pca %>%
    select(Sample.x, enroldate:ncol(.)) %>%
    rename(Sample = Sample.x)

# Summarise sample clusters
metadata  %>% 
    group_by(Cluster) %>%
    summarise(n = n()) %>% 
    write_tsv("outputs/sample_cluster_sumamry.tsv")



############################################################################################################################################################

# MDS - clusters

## Read in data
mds <- read_table("Pk.mds", col_names=T) %>%
    select(-c("FID", "X14")) %>%
    rename("sampleid" = "IID") %>%
    rename_at(vars(starts_with("C")), ~str_replace(., "C", "MDS"))

## Add metadata 
mds <- mds %>%
    rename(Sample = sampleid) %>% 
    left_join(metadata)

## Plot MDS
mds_plot <- ggplot(mds, aes(MDS1, MDS2, colour = Cluster)) + 
    geom_point(size = 3) +
    scale_color_viridis_d() +
    coord_equal() + 
    theme_light()

ggsave("outputs/mds_cluster.png", dpi=600, mds_plot)

############################################################################################################################################################


# Neighbour-joining tree

## Create distance matrix from PLINK data and build NJT
NJT_ID <- read_table("Pk.dist.id", col_names=F) %>%
    as.data.frame() %>%
    mutate_all(~str_remove(., "_DK.*")) 

NJT_matrix <- read_table("Pk.dist", col_names=NJT_ID$X1) %>%
        as.data.frame() %>%
        add_column(Row_Names = NJT_ID$X1) %>%
        column_to_rownames("Row_Names") %>%
        as.matrix()

## Plot tree

#Plot 
options(ignore.negative.edge=TRUE)

# Clusters
# Read in metadata and assign colours - based on clusters
NJT_metadata <- metadata %>% 
    relocate(Sample) %>% 
    mutate(Sample = str_remove(Sample, "_DK.*")) %>% # for compatbility of datasets
    mutate(Cluster = ifelse(is.na(Cluster), "Missing", .$Cluster)) %>% # if not reassigned the final argument of the next line does not work
    mutate(Colour = ifelse(Cluster == "Mf", "#440154FF", ifelse(Cluster == "Mn", "#39568CFF", ifelse(Cluster == "Peninsular", "#1F968BFF", "#808080"))))

# Generate NJT and reorder outer to inner edges
NJT_tree <- nj(NJT_matrix) %>% 
       reorder.phylo(order = "postorder") 

# Get tip colours by joining metadata and tip labels
tip_colours <- NJT_tree$tip.label %>% 
    as.data.frame() %>% 
    rename(Sample = ".") %>%
    left_join(
        NJT_metadata %>% 
            select(Sample, Colour)
        ) %>% 
    select(Colour) %>% 
    mutate(Colour = ifelse(is.na(Colour), "#808080", .$Colour)) %>%
    pull()

# Combine tip colours and edge IDs - gives list of colours for each edge, with #### representing missing values (ie internal edges)
edge_colours <- tip_colours[NJT_tree$edge[, 2]]
edge_colours[is.na(edge_colours)] <- "####"

# Loop through each edge and assign colours to internal edges dependent on the nodes that branch from it (ie if two nodes/edges of the same colour branch from the same node/edge, then assign that branch with the colour of the two nodes/edges)
for (edge in 1:nrow(NJT_tree$edge)) {
    if (edge_colours[edge] == "####") { 
        nodes <- which(NJT_tree$edge[, 1] == NJT_tree$edge[edge, 2]) # creates object of all the nodes that have multiple branches coming from them
        if (edge_colours[nodes[1]] == edge_colours[nodes[2]]) { # if the colour of the node in the first column = the second column 
            edge_colours[edge] <- edge_colours[nodes[1]] # make the edge that they branch from the same colour
            } else {
            edge_colours[edge] <- "#808080" # otherwise we make it grey - this object has ALL colours for ALL nodes
            }
        }
    }

NJT_tree_colours <- edge_colours %>% 
    as.data.frame() %>%
    cbind(as.data.frame(NJT_tree$edge[, 2])) %>% # joining the resulting colours (for ALL nodes) to the edges data in the original tree - so that when we plot it assigns the correct colour to the correct edge
    rename(Colour = ".", node = "NJT_tree$edge[, 2]") %>% 
    relocate(node) # might not be needed...

NJT_tree_plot <- ggtree(NJT_tree, layout = "daylight", size = 0.5) %<+% NJT_tree_colours + 
    aes(colour = I(Colour)) 

ggsave("outputs/njt_tree_cluster_unrooted.png", dpi = 300, NJT_tree_plot)


# Tree with sevemal

NJT_metadata <- metadata %>% 
    relocate(Sample) %>% 
    mutate(Sample = str_remove(Sample, "_DK.*")) %>% # for compatbility of datasets
    mutate(sevemal = ifelse(is.na(sevemal), "Missing", .$sevemal)) %>% # if not reassigned the final argument of the next line does not work
    mutate(Colour = ifelse(sevemal == "SM", "#440154FF", ifelse(sevemal == "UM", "#39568CFF", "#808080")))

# Generate NJT and reorder outer to inner edges
NJT_tree <- nj(NJT_matrix) %>% 
       reorder.phylo(order = "postorder") 

# Get tip colours by joining metadata and tip labels
tip_colours <- NJT_tree$tip.label %>% 
    as.data.frame() %>% 
    rename(Sample = ".") %>%
    left_join(
        NJT_metadata %>% 
            select(Sample, Colour)
        ) %>% 
    select(Colour) %>% 
    mutate(Colour = ifelse(is.na(Colour), "#808080", .$Colour)) %>%
    pull()

# Combine tip colours and edge IDs - gives list of colours for each edge, with #### representing missing values (ie internal edges)
edge_colours <- tip_colours[NJT_tree$edge[, 2]]
edge_colours[is.na(edge_colours)] <- "####"

# Loop through each edge and assign colours to internal edges dependent on the nodes that branch from it (ie if two nodes/edges of the same colour branch from the same node/edge, then assign that branch with the colour of the two nodes/edges)
for (edge in 1:nrow(NJT_tree$edge)) {
    if (edge_colours[edge] == "####") { 
        nodes <- which(NJT_tree$edge[, 1] == NJT_tree$edge[edge, 2]) # creates object of all the nodes that have multiple branches coming from them
        if (edge_colours[nodes[1]] == edge_colours[nodes[2]]) { # if the colour of the node in the first column = the second column 
            edge_colours[edge] <- edge_colours[nodes[1]] # make the edge that they branch from the same colour
            } else {
            edge_colours[edge] <- "#808080" # otherwise we make it grey - this object has ALL colours for ALL nodes
            }
        }
    }

NJT_tree_colours <- edge_colours %>% 
    as.data.frame() %>%
    cbind(as.data.frame(NJT_tree$edge[, 2])) %>% # joining the resulting colours (for ALL nodes) to the edges data in the original tree - so that when we plot it assigns the correct colour to the correct edge
    rename(Colour = ".", node = "NJT_tree$edge[, 2]") %>% 
    relocate(node) # might not be needed...

NJT_tree_plot <- ggtree(NJT_tree, layout = "daylight", size = 0.5) %<+% NJT_tree_colours + 
    aes(colour = I(Colour)) 

ggsave("outputs/njt_tree_sevemal_unrooted.png", dpi = 300, NJT_tree_plot)


# Tree with sevemal
NJT_metadata <- metadata %>%
    mutate(parasitemia = as.numeric(parasitemia)) %>%
    mutate(bin = cut_interval(parasitemia, 20)) %>%
    group_by(bin) %>%
    add_column(bins = 1:nrow(.)) %>%
    ungroup() 

NJT_tree_plot <- ggtree(NJT_tree, layout="daylight", size = 0.5, aes(colour = bins)) %<+% NJT_metadata +
    theme(legend.position = "right", 
    legend.title = element_blank(), 
    legend.key = element_blank()) +
    scale_color_viridis_c()
    
ggsave("outputs/njt_tree_para_unrooted.png", dpi = 300, NJT_tree_plot)
