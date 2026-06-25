# Plotting
library(igraph)
library(tidyverse)
#library(MetamapsDB)
library(viridis)
library(janitor)

# Metadata

## Read in cluster data (based on ADMIXTURE)
metadata <- read_tsv("outputs/admix_clusters.tsv")

## Add full metadata
metadata_full <- read_tsv("/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/data/metadata/full_metadata.tsv") %>% 
    mutate(sevemal = ifelse(grepl("UM", sevemal), "UM", .$sevemal)) %>% 
    mutate(sevemal = ifelse(sevemal == "1", "SM", ifelse(sevemal == "0", "UM", .$sevemal))) %>% 
    rename(Sample = studycode) %>% 
    mutate(district = str_to_title(district)) %>%
    add_column(Cluster = NA)

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


metadata_full <- metadata_full %>%
    rbind(archived_data) %>%
  mutate_at(c("district", "village"), str_to_title) 


Pk.IBD.file <- "hmmIBD/no_lab_strains/Pk.hmm_fract.txt"
Pk.IBD <- read.delim(Pk.IBD.file) 
  #mutate(sample1 = str_remove(sample1, "_DK.*")) %>% 
  #mutate(sample2 = str_remove(sample2, "_DK.*")) 

# Plot Fraction of sites that are IBD
IBD_meta_combined <- Pk.IBD %>%
    as.data.frame() %>%
    select(sample1, sample2, fract_sites_IBD) %>%
    left_join(
      metadata %>% mutate(sample1 = Sample) %>% select(sample1, Cluster, Proportion),
      by = "sample1"
    ) %>% 
    left_join(
      metadata %>% mutate(sample2 = Sample) %>% select(sample2, Cluster, Proportion),
      by = "sample2" 
    )
  
############################### WIP

# Investigate samples of high relatedness within clusters 
clonal_samples <- function(CLUSTER, THRESHOLD){
  IBD_meta_combined %>% 
    filter(Cluster.x == CLUSTER & Cluster.y == CLUSTER & fract_sites_IBD >= THRESHOLD) %>% 
    select(sample1, sample2, fract_sites_IBD) %>% 
    mutate(sample1 = str_remove(sample1, "_DK.*"), sample2 = str_remove(sample2, "_DK.*")) %>% 
    mutate(sample1 = str_replace(sample1, "_", "-"), sample2 = str_replace(sample2, "_", "-")) %>% 
    left_join(
      metadata_full %>%
        select(Sample, enroldate, district, village, occupation, parasitemia, age, sex) %>% 
        rename_all(~tolower(.)) %>% 
        rename_all(~paste0(., "1"))
    ) %>% 
    left_join(
      metadata_full %>%
        select(Sample, enroldate, district, village, occupation, parasitemia, age, sex) %>% 
        rename_all(~tolower(.)) %>% 
        rename_all(~paste0(., "2"))
    ) %>%
    add_column(clonal_cluster = CLUSTER)
}

clonal_samples("Mn", 0.3189) 


clonal_samples("Mn", 0.95) %>% 
  rbind(
    clonal_samples("Mf", 0.95)
  ) %>%
  rbind(
    clonal_samples("Peninsular", 0.95)
  ) %>% 
  write_tsv("outputs/IBD_clonal_transmission.tsv")