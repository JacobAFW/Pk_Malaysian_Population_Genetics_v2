# Load packages
library(tidyverse)

# Read in metadata
metadata <- read_tsv("/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/outputs/05_Analyses/Population_genetics/outputs/admix_clusters.tsv") %>%
    dplyr::select(Sample, Cluster)

############################################################################## WORK IN PROGRESS - using VCF to identify dominant alleles
# Get allele call for each sample
# Add metadata/cluster information
# Create column that has dominant allele call for each cluster at each position in the chromosome 
# use ifelse statement - for each sample, does it have the Mn or Mf allele

exclude_samples <- read_tsv("/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/outputs/05_Analyses/Population_genetics/introgression/exclude.txt", col_names = FALSE) 
genotype_table <- read_tsv("/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/outputs/05_Analyses/Population_genetics/introgression/hmmIBD.tsv") %>%
    dplyr::select(!all_of(exclude_samples$X1))

# create long version of genotype data
genotype_table_long <- genotype_table %>%
    pivot_longer(3:ncol(.), names_to = "SAMPLE", values_to = "SNP") %>%
    left_join(
        metadata %>%
            rename(SAMPLE = Sample)
    ) #%>%
    #filter(CHROM == "08") # removed when finish


dominant_allele <- genotype_table_long %>% 
    filter(SNP >= 0) %>% # remove missing data
    mutate(SNP = as.factor(SNP)) %>% # convert to factor so we can create a summary of the observations
    group_by(CHROM, POS, Cluster, SNP) %>% 
    summarise(Allele_count = n()) %>% # get the counts for each unique combination you've grouped by
    #filter(POS==56510) %>%
    filter(Allele_count == max(Allele_count)) %>% # filter to get the allele with max count for each combination (including the SNP/allele)
    select(-5) %>% # remove the allele count column
    ungroup() %>% 
    group_by(CHROM, POS) %>% 
    mutate(row = cur_group_id()) %>% # create an id specific to each combo of CHROM & POS
    na.omit() %>% # remove rows with missing data
    group_by(row) %>% # group by the new id
    summarise(CHROM, POS, Cluster, SNP, n = n()) %>% # get the number of rows for this new id - any with 4 or more have duplicate clusters for the same position, ie they have a 50-50 split in allele frequency 
    filter(n < 4) %>% # filter out positions that have a cluster that has 50-50 allele split
    pivot_wider(names_from = Cluster, values_from = SNP) %>% 
    ungroup() %>%
    dplyr::select(-c("row", "n"))

# ASSUMPTION = if it is not equal to the dominant allele in any of the clusters, just make it the allele of the cluster it belongs to - this only occurs in ~0.4% of SNPs
introgression_table <- genotype_table_long %>%
    left_join(dominant_allele) %>%
    add_column(Intro_clust = .$SNP) %>%
    mutate(Intro_clust = ifelse(Cluster == "Mn" & Intro_clust == Mn, "Mn", # Mn
        ifelse(Cluster == "Mn" & Intro_clust == Mf & Intro_clust != Peninsular, "Mf", 
        ifelse(Cluster == "Mn" & Intro_clust == Peninsular & Intro_clust != Mf, "Peninsular", 
        ifelse(Cluster == "Mn" & Intro_clust == Mf & Intro_clust == Peninsular, "Mf_Peninsular",
        ifelse(Cluster == "Mn" & (Intro_clust != Mn | Intro_clust != Mf | Intro_clust != Peninsular), "Mn", .$Intro_clust)))))) %>% 
    mutate(Intro_clust = ifelse(Cluster == "Mf" & Intro_clust == Mf, "Mf", # Mf
        ifelse(Cluster == "Mf" & Intro_clust == Mn & Intro_clust != Peninsular, "Mn", 
        ifelse(Cluster == "Mf" & Intro_clust == Peninsular & Intro_clust != Mn, "Peninsular", 
        ifelse(Cluster == "Mf" & Intro_clust == Mn & Intro_clust == Peninsular, "Mn_Peninsular",
        ifelse(Cluster == "Mf" & (Intro_clust != Mn | Intro_clust != Mf | Intro_clust != Peninsular), "Mf", .$Intro_clust)))))) %>%
    mutate(Intro_clust = ifelse(Cluster == "Peninsular" & Intro_clust == Peninsular, "Peninsular", # Peninsular
        ifelse(Cluster == "Peninsular" & Intro_clust == Mf & Intro_clust != Mn, "Mf",
        ifelse(Cluster == "Peninsular" & Intro_clust == Mn & Intro_clust != Mf, "Mn", 
        ifelse(Cluster == "Peninsular" & Intro_clust == Mf & Intro_clust == Mn, "Mf_Mn",
        ifelse(Cluster == "Peninsular" & (Intro_clust != Mn | Intro_clust != Mf | Intro_clust != Peninsular), "Peninsular", .$Intro_clust)))))) 

write_tsv(introgression_table, "introgression_updated/introgression_table.tsv")
#introgression_table <- read_tsv("introgression_updated/introgression_table.tsv")

# genetic distance per sliding window of each sample to each cluster (dominant allele) - the number of times the SNP is not the dominant allele of that cluster

window_size <- 10000 # should give us ~180 windows per sample per chrom

#WINDOWS Table for downstream
PKA1H1_windows <- read.table("/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/data/ref_genomes/PKA1H1/strain_A1_H.1.Icor.fasta.bed") %>% 
    mutate(CHROM = str_remove(V1, "ordered_PKNH_")) %>% 
    mutate(CHROM = str_remove(CHROM, "_v2")) %>% 
    rename(Start = V2, End = V3) %>% 
    select(-V1) %>% 
    filter(CHROM != "PKNH_MIT", CHROM != "new_API_strain_A1_H.1") 

library(plyr)
PKA1H1_windows <- ddply(PKA1H1_windows, "CHROM", summarise, POS = seq(Start, End)) %>% 
    mutate(TMP_WINDOW = (floor(POS/window_size) * window_size) + (window_size/2)) %>% 
    mutate(WINDOW = group_indices(., CHROM, TMP_WINDOW)) %>% 
    select(-TMP_WINDOW)
detach("package:plyr", unload=TRUE)

write_tsv(PKA1H1_windows , "introgression_updated/PKA1H1_windows.tsv")
#PKA1H1_windows <- read_tsv("introgression_updated/PKA1H1_windows.tsv")


introgression_table_window <- introgression_table %>% 
    left_join(PKA1H1_windows) %>% 
    filter(SNP >= 0) %>% # filter out missing calls - messes with the distance calculation
    group_by(SAMPLE, WINDOW) %>%
    mutate(n = n()) %>% # filter windows by the number of SNPs in a window (5) - if < 5, how reliable is the genetic distance?
    filter(n > 5) %>%
    summarise(Mf_distance = sum(SNP != Mf, na.rm = TRUE)/n()*100, # the number of times the sample allele doesn't match the dominant allele = genetic distance
        Mn_distance = sum(SNP != Mn, na.rm = TRUE)/n()*100,
        Pen_distance = sum(SNP != Peninsular, na.rm = TRUE)/n()*100)  %>%
    left_join(
        metadata %>% rename(SAMPLE = Sample)
    ) %>%
    filter(SAMPLE %in% metadata$Sample) %>% 
    filter(Mn_distance != Mf_distance) # filter out windows where samples have the same distance for both of the major clusters

write_tsv(introgression_table_window, "introgression_updated/introgression_table_window.tsv")
#introgression_table_window <- read_tsv("introgression_updated/introgression_table_window.tsv")


######################################## CREATE POLYGON PLOTS TO REPRESENT INTROGRESSION PER SAMPLE ########################################



# Plot function
## for unfiltered data - entire genome (doesn't contain point labels - too messy)
genetic_distance_plot <- function(data, sample){
subset <- data %>%  # for title
    filter(SAMPLE == sample)

data %>%
    ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
    geom_point(data = introgression_table_window %>%
        filter(SAMPLE == sample),
        size = 0.75, alpha = 0.75) +
    geom_density_2d_filled(mapping = aes(x = Mf_distance, y = Mn_distance, alpha = (..level..), fill = Cluster),
        data = introgression_table_window,
        contour_var = "density") +
    scale_fill_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) +
    scale_alpha_discrete(guide = "none") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    ggtitle(paste0(sample, ", ", subset[1,"Cluster"], ", ", subset[1,"Location"])) 
}

## for data filtered to a particular CHROM (contain point labels with window number) - filtering must be done prior - line 97
genetic_distance_plot_chrom <- function(data, sample){

subset <- data %>%  # subset to deal with smaller dataset
    filter(SAMPLE == sample)

max_value <- subset %>% # get max values to use for x and y limits
    summarise(Mf = max(Mf_distance), Mn = max(Mn_distance)) 

subset %>% # create plot
    ggplot(aes(x = Mf_distance, y = Mn_distance, colour = Cluster)) +
        geom_point() +
        ggtitle(paste0(sample, ", ", subset[1,"Cluster"], ", ", subset[1,"Location"])) + # extract cluster and location information to include in title
        ylim(0, ifelse(max_value$Mf > max_value$Mn, max_value$Mf + 5, max_value$Mn + 5)) + # use the max values calculated above to inform the axis values
        xlim(0, ifelse(max_value$Mf > max_value$Mn, max_value$Mf + 5, max_value$Mn + 5)) +
        theme(legend.position = "none") +
        scale_colour_manual(values = ifelse(subset$Cluster == "Mf", "#440154FF", ifelse(subset$Cluster == "Mn", "#39568CFF", "#73D055FF") )) + # colour based on cluster
        geom_vline(xintercept = 100, linetype = "dashed", size = 1) +
        geom_hline(yintercept = 100, linetype = "dashed", size = 1) +
        geom_text(aes(label = ifelse(Mn_distance > 100, as.character(WINDOW), "")), hjust = 0, vjust = 2) + # add 'window' text to values in outlier quandrants
        geom_text(aes(label = ifelse(Mf_distance > 100, as.character(WINDOW), "")), hjust = 0, vjust = 2)
}


# Create a loop to plot all samples that are Mn or Mf
sample_names <- introgression_table_window %>% 
    filter(Cluster == "Mn" | Cluster == "Mf") %>%
    dplyr::select(SAMPLE) %>%
    unique() %>%
    mutate(SAMPLE = as.factor(SAMPLE)) 




######################################## USE CONTOURS TO IDENTIFY INTROGRESSED REGIONS ########################################






# Function for identifying regions of introgression from a contours and point plot 
find_introgressed_regions <- function(SAMPLENAME){
    # libaries
    library(sp) 
    library(MASS)
    library(tidyverse)

    # Create plot with points and polygons
    raster_plot <- introgression_table_window  %>%
        ggplot(aes(x = Mf_distance, y = Mn_distance, group = Cluster)) +
        geom_point(data = introgression_table_window %>%
            filter(SAMPLE == SAMPLENAME),
            size = 0.75, alpha = 0.75) +
        geom_density_2d(mapping = aes(x = Mf_distance, y = Mn_distance, colour = Cluster),
            data = introgression_table_window,
            contour_var = "density") +
        scale_colour_manual(values = c("#440154FF", "#73D055FF", "#39568CFF")) +
        coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) 

    # Extract info on points and polygons from ggplot
    polygon_data <- ggplot_build(raster_plot) # create matrix for contours based on plot above - includes all datapoints
    polygon_points <- polygon_data$data[[1]]
    polygon_contours <- polygon_data$data[[2]] # extract contour data

    # Filter to cluster-specific polygon/contours and identify whether or not points fall within that cluster
    Mf_contours <- polygon_contours %>% 
        filter(colour == "#440154FF") %>% 
        filter(level > 5e-4) # greater than the lowest contour

    Mf_contours <- point.in.polygon(pol.x = as.numeric(unlist(Mf_contours[3])), 
        pol.y = as.numeric(unlist(Mf_contours[4])),
        point.x = as.numeric(unlist(polygon_points[1])), 
        point.y = as.numeric(unlist(polygon_points[2]))) %>% 
        as.data.frame() %>% 
        rename("Mf" = ".")

    Mn_contours <- polygon_contours %>% 
        filter(colour == "#73D055FF") %>% 
        filter(level > 5e-4)

    Mn_contours <- point.in.polygon(pol.x = as.numeric(unlist(Mn_contours[3])), 
        pol.y = as.numeric(unlist(Mn_contours[4])),
        point.x = as.numeric(unlist(polygon_points[1])), 
        point.y = as.numeric(unlist(polygon_points[2]))) %>% 
        as.data.frame() %>% 
        rename("Mn" = ".")

    Pen_contours <- polygon_contours %>% 
        filter(colour == "#39568CFF")  %>% 
        filter(level > 5e-4)

    Pen_contours <- point.in.polygon(pol.x = as.numeric(unlist(Pen_contours[3])), 
        pol.y = as.numeric(unlist(Pen_contours[4])),
        point.x = as.numeric(unlist(polygon_points[1])), 
        point.y = as.numeric(unlist(polygon_points[2]))) %>% 
        as.data.frame() %>% 
        rename("Pen" = ".")


    # Combine newly generated point presence-abscence info together
    points_in_contours_data <- Mf_contours %>% 
        cbind(Mn_contours) %>% 
        cbind(Pen_contours) 

    # Combine with original introgression table to filter to "introgressed" points - those that fall within another clusters polygon but not its own
    subset <- introgression_table_window %>%
        filter(SAMPLE == SAMPLENAME) %>% 
        cbind(points_in_contours_data) %>% 
        filter(ifelse(Cluster == "Mn", Mf == 1 & Mn == 0 & Pen == 0, # Mf to Mn introgression
            ifelse(Cluster == "Mf", Mf == 0 & Mn == 1 & Pen == 0, # Mn to Mf introgression
            ifelse(Cluster == "Peninsular", Mf == 1 & Mn == 1 & Pen == 0)))) # Mn or Mf to Peninsular introgression

    return(subset)
}


# Create a loop to identify introgressed windows
## extract sample names
sample_names <- introgression_table_window %>% 
    dplyr::select(SAMPLE) %>%
    unique() %>%
    mutate(SAMPLE = as.factor(SAMPLE)) 
     
## create empty df
introgressed_windows <-  introgression_table_window %>% 
    add_column(Mf = 1, Mn = 1, Pen = 1) %>% 
    slice(-(1:nrow(.)))

## loop into df
for(i in levels(sample_names$SAMPLE)){
    introgressed_windows <- introgressed_windows %>%
        rbind(find_introgressed_regions(paste0(i)))
}      


write_tsv(introgressed_windows, "introgression_updated/introgressed_windows_updated.tsv")
#introgressed_windows <- read_tsv("introgression_updated/introgressed_windows.tsv")



######################################## APPLY FILTERS ########################################



# Filter out low n windows - full datatset

## shoulder plot to come up with n threshold 
shoulder_plot <- introgressed_windows %>%
    dplyr::select(WINDOW) %>% 
    group_by(WINDOW) %>% 
    summarise(n = n()) %>%
    arrange(n) %>% 
    add_column(WINDOWS = 1:nrow(.)) %>%
    ggplot(aes(x = WINDOWS, y = n)) +
    geom_col() + 
    scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,20,30,40))

ggsave("introgression_updated/shoulder_plot.png", dpi = 300,  shoulder_plot)

# filter out windows that only appear below a treshold of samples
introgressed_windows_filter <- introgressed_windows %>%
    dplyr::select(WINDOW) %>% 
    group_by(WINDOW) %>% 
    summarise(n = n()) %>% 
    #summarise(n = n()/152*100) %>% # HARD CODED 
    filter(n > 5) 

introgressed_windows <- introgressed_windows %>% 
    filter(WINDOW %in% introgressed_windows_filter$WINDOW) 

# Filter out low n windows - per cluster

Mf_filter <- introgressed_windows %>%
    filter(Cluster == "Mf") %>%
    dplyr::select(WINDOW) %>% 
    group_by(WINDOW) %>% 
    summarise(n = n()) %>% 
    #summarise(n = n()/81*100) %>%  # HARD CODED 
    filter(n > 5) 


Mn_filter <- introgressed_windows %>%
    filter(Cluster == "Mn") %>%
    dplyr::select(WINDOW) %>% 
    group_by(WINDOW) %>% 
    summarise(n = n()) %>% 
    #summarise(n = n()/81*100) %>%  # HARD CODED 
    filter(n > 5) 

Pen_filter <- introgressed_windows %>%
    filter(Cluster == "Peninsular") %>%
    dplyr::select(WINDOW) %>% 
    group_by(WINDOW) %>% 
    summarise(n = n()) %>%
    #summarise(n = n()/81*100) %>%  # HARD CODED 
    filter(n > 5) 

introgressed_windows <- introgressed_windows %>% 
    filter(Cluster == "Mf") %>% 
    filter(WINDOW %in% Mf_filter$WINDOW) %>% 
    rbind(
        introgressed_windows %>% 
            filter(Cluster == "Mn") %>%
            filter(WINDOW %in% Mn_filter$WINDOW)
    ) %>%
    rbind(
        introgressed_windows %>% 
            filter(Cluster == "Peninsular") %>%
            filter(WINDOW %in% Pen_filter$WINDOW)
    )


# Filter out regions with SICAvar or Kir genes - these should be missing already, deoending on hwo the data was subset before the inrogression analysis

## Read in GFF file
GFF_annotation <- read.table("/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/data/ref_genomes/PKA1H1/gff/strain_A1_H.1.Icor.gff3", skip = 18, sep = "\t") %>%  
    rename(CHROM = V1,
    source = V2,
    feature= V3,
    start = V4,
    end = V5,
    score = V6,
    strand = V7,
    frame = V8,
    attribute = V9) %>% 
    filter(grepl("ordered", CHROM)) %>% 
    mutate(CHROM = str_remove(CHROM, "ordered_PKNH_")) %>% 
    mutate(CHROM = str_remove(CHROM, "_v2"))

SICAvar_KIR <- GFF_annotation %>% 
    filter(grepl("SICA", attribute)) %>% 
    rbind(
        GFF_annotation %>% 
            filter(grepl("KIR", attribute))
    ) %>% 
    as_tibble() %>% 
    pivot_longer(c("start", "end"), names_to = "range", values_to = "POS") %>% 
    unite(CHROM_POS, c("CHROM", "POS"), sep = "_") %>% 
    dplyr::select(CHROM_POS) %>%
    unique()

SICAvar_KIR <- PKA1H1_windows %>% 
    unite(CHROM_POS, c("CHROM", "POS"), sep = "_") %>% 
    filter(CHROM_POS %in% SICAvar_KIR$CHROM_POS)

introgressed_windows <- introgressed_windows %>% 
    filter(!(WINDOW %in% SICAvar_KIR$WINDOW)) # filter out known hypervariable regions - removes another ~1K windows

# Filter windows that appear in multiple clusters as potential introgression events - suggests hypervariability

hypervariable_filter <- introgressed_windows %>%
    dplyr::select(WINDOW, Cluster) %>% 
    unique() %>% 
    group_by(WINDOW) %>% 
    summarise(n = n()) %>% 
    filter(n > 1)

introgressed_windows <- introgressed_windows %>% 
    filter(!(WINDOW %in% hypervariable_filter$WINDOW)) # filter out unknown hypervariable regions - those windows that appear in multiple clusters


write_tsv(introgressed_windows, "introgression_updated/introgressed_windows_filtered.tsv")


######################################## Explore windows of introgression ########################################


    
# Summary info for TOTAL introgression across clusters
introgressed_windows %>% 
    dplyr::select(SAMPLE, Cluster) %>%
    group_by(SAMPLE, Cluster) %>% 
    summarise(n = n()) %>% 
    ungroup() %>% 
    group_by(Cluster) %>% 
    summarise(mean = mean(n), median = median(n), sd = sd(n), n = n()) %>% # number of windows
    arrange(desc(median)) %>% write_tsv("introgression_updated/average_windows_for_clusters.tsv")

# Number of introgression events per sample
intro_per_sample_summary <- introgressed_windows %>% 
    group_by(SAMPLE) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n)) %>% 
    mutate(intro_level = ntile(n, 3)) %>% 
    mutate(intro_level = ifelse(intro_level == 3, "High", ifelse(intro_level == 2, "Medium", "Low"))) 

intro_per_sample_summary %>% 
    write_tsv("introgression_updated/intro_per_sample_summary.tsv")

######################## Introgressed windows for each cluster


introgressed_windows %>% 
    filter(Cluster == "Mf" & Mf == 0 & Mn == 1) %>%
    group_by(WINDOW) %>% 
    summarise(n = n()) %>% 
    left_join(
        PKA1H1_windows %>% 
            group_by(CHROM, WINDOW) %>% 
            summarise(start = min(POS), end = max(POS))
    ) %>%
    arrange(desc(n)) %>% write_tsv("introgression_updated/mf_windows.tsv")

introgressed_windows %>% 
    filter(Cluster == "Mf" & Mf == 0 & Mn == 1) %>%
    group_by(SAMPLE) %>% 
    summarise(n = n()) %>%
    arrange(desc(n)) %>% 
    write_tsv("introgression_updated/mf_samples_window_counts.tsv")

introgressed_windows %>% 
    filter(Cluster == "Mn" & Mf == 1 & Mn == 0) %>%
    group_by(WINDOW) %>% 
    summarise(n = n()) %>% 
    left_join(
        PKA1H1_windows %>% 
            group_by(CHROM, WINDOW) %>% 
            summarise(start = min(POS), end = max(POS))
    ) %>%
    arrange(desc(n)) %>% write_tsv("introgression_updated/mn_windows.tsv")

introgressed_windows %>% 
    filter(Cluster == "Mn" & Mf == 1 & Mn == 0) %>%
    group_by(SAMPLE) %>% 
    summarise(n = n()) %>%
    arrange(desc(n)) %>% 
    write_tsv("introgression_updated/mn_samples_window_counts.tsv")
    
# Mf    
introgression_of_Mn_into_Mf <- introgressed_windows %>% # position of windows
    filter(Cluster == "Mf" & Mf == 0 & Mn == 1) %>% 
    dplyr::select(WINDOW) %>% 
    unique() %>% 
    dplyr::left_join(PKA1H1_windows)

### introgressed windows per chrom
introgression_of_Mn_into_Mf %>% 
    dplyr::select(WINDOW, CHROM) %>% 
    unique() %>% 
    group_by(CHROM) %>% 
    summarise(n_windows = n()) %>% 
    arrange(desc(n_windows)) %>% write_tsv("introgression_updated/Mf_windows_across_chrom.tsv")

# Summary table of samples and windows - presence/absence of introgression
Mf_introgression_matrix <- introgressed_windows %>% 
    filter(Cluster == "Mf" & Mf == 0 & Mn == 1)  %>%
    dplyr::select(SAMPLE, WINDOW) %>% 
    add_column(PRESENCE = 1) %>% 
    mutate(WINDOW = as.factor(as.character(WINDOW))) %>%
    pivot_wider(names_from = SAMPLE, values_from = PRESENCE) %>% 
    replace(is.na(.), 0 )


# Mn    
introgression_of_Mf_into_Mn <- introgressed_windows %>% # position of windows
    filter(Cluster == "Mn" & Mf == 1 & Mn == 0) %>% 
    dplyr::select(WINDOW) %>% 
    unique() %>% 
    dplyr::left_join(PKA1H1_windows)

### introgressed windows per chrom
introgression_of_Mf_into_Mn %>% 
    dplyr::select(WINDOW, CHROM) %>% 
    unique() %>% 
    group_by(CHROM) %>% 
    summarise(n_windows = n()) %>% 
    arrange(desc(n_windows)) %>% write_tsv("introgression_updated/Mn_windows_across_chrom.tsv")

# Summary table of samples and windows - presence/absence of introgression
Mf_introgression_matrix <- introgressed_windows %>% 
    filter(Cluster == "Mn" & Mf == 1 & Mn == 0) %>% 
    dplyr::select(SAMPLE, WINDOW) %>% 
    add_column(PRESENCE = 1) %>% 
    mutate(WINDOW = as.factor(as.character(WINDOW))) %>%
    pivot_wider(names_from = SAMPLE, values_from = PRESENCE) %>% 
    replace(is.na(.), 0 )
