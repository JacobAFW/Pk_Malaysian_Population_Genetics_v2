# Load packages
library(tidyverse)

# Read in data
Pk_maf <- read_table("Pk.frq", col_names=T) %>%
    mutate(CHR = str_remove(SNP, ":.*"))

# Density of MAF by SNP
Pk_maf_density <- Pk_maf  %>%
    ggplot(aes(x = MAF)) +
    geom_density(alpha=.3) 

ggsave("filtering/Pk_maf_density.png", dpi = 600, Pk_maf_density)

# Plot MAF by chr
Pk_maf_chr <- Pk_maf %>% 
    mutate(SNP = str_remove(SNP, "ordered_PKNH_"),
        SNP = str_remove(SNP,"_v2.*")) %>% 
    group_by(SNP) %>%
    summarise(MAF = mean(MAF)) %>%
    ggplot(aes(x = SNP, y = MAF)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 45), legend.position = "none") +
    scale_fill_viridis_d() 

ggsave("filtering/Pk_maf_plot_chr.png", dpi = 600, Pk_maf_chr)

# Summary stats

Pk_maf %>% 
    mutate(SNP = str_remove(SNP, "ordered_PKNH_"),
        SNP = str_remove(SNP,"_v2.*")) %>% 
    group_by(SNP) %>%
    summarise(MAF_mean = mean(MAF),
            MAF_SD = sd(MAF),
            MAF_SE = (sd(MAF))/sqrt(n())) %>%
rbind(
    Pk_maf %>%
        summarise(MAF_mean = mean(MAF),
                MAF_SD = sd(MAF),
                MAF_SE = (sd(MAF))/sqrt(n())) %>%
        add_column(SNP = "Total")
    ) %>%
    write_csv(col_names = T, "filtering/Pk_maf_summary.csv")
