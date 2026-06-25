# Load packages
library(tidyverse)
library(data.table)

# Sample-based missing data
## Read in data
sample_mis <- read_table("Pk.imiss", col_names=T) 

# Density of miss by sample
sample_miss_density <- sample_mis %>% 
    rename(Sample = IID,
        Geno_Miss = F_MISS) %>%
    select(2:ncol(.)) %>%
    ggplot(aes(x = Geno_Miss)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept = 0.05, colour = "#1F968BFF") 

ggsave("filtering/sample_miss_density_plot.png", width = 12, dpi = 600, sample_miss_density)

# Mean missingness
metadata <- read_tsv("/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/data/metadata/full_metadata.tsv") %>% 
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

metadata <- metadata %>%
    rbind(archived_data)


sample_mis %>%
    left_join(
        metadata %>%
            rename(FID = Sample)
    ) %>%
    summarise(mean = mean(F_MISS), sd = sd(F_MISS), min = min(F_MISS), max = max(F_MISS)) %>% 
    write_csv("filtering/missingness.csv")


# Variant-based missing data
## Read in data
var_miss <- read_table("Pk.lmiss", col_names=T) %>%
    mutate(CHR = str_remove(SNP, ":.*")) %>%
    rename(Geno_Miss = F_MISS,
            Chr = CHR)

# Density of miss by variant
var_miss_density <- var_miss %>% 
    ggplot(aes(x = Geno_Miss)) +
    geom_density(alpha=.3) +
    geom_vline(xintercept = 0.2, colour = "#1F968BFF") +
    geom_vline(xintercept = 0.1, colour = "#1F968BFF") +
    geom_vline(xintercept = 0.05, colour = "#1F968BFF") 

ggsave("filtering/variant_miss_density_plot.png", width = 12, dpi = 600, var_miss_density)


# Summary stats

miss_summary <- var_miss %>%
    summarise(Miss_mean = mean(Geno_Miss),
                Miss_SD = sd(Geno_Miss),
                Miss_SE = (sd(Geno_Miss))/sqrt(n()))

write_csv(miss_summary, "filtering/Pk_miss_summary.csv")