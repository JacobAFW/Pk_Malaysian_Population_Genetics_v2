# Packages
library(dplyr)
library(readr)
library(stringr)
library(data.table)

# Read in and wrangle VCF

read_tsv("hmmIBD/sample_miss_5/grep_patterns.tsv", col_names = c("CHROM", "POS")) %>%
        mutate(POS = as.numeric(POS)) %>%
        left_join(
                read_table("hmmIBD/Consensus_SNPs_subset_no_MOI.vcf", skip = 72) %>%
                rename("CHROM" = `#CHROM`) %>% # this first section selects only variants that passed the PLINK filters we applied
                mutate(POS = as.numeric(POS)) 
        ) %>%
        mutate(CHROM = str_remove(CHROM, "ordered_PKNH_"),
        CHROM = str_remove(CHROM, "_v2")) %>%
        mutate_at(c(10:ncol(.)), ~str_remove(., ":.*")) %>% 
        mutate_at(c(10:ncol(.)), ~ifelse(. %like% "1/1" | . %like% "1/0" | . %like% "0/1" , "1", .)) %>% 
        mutate_at(c(10:ncol(.)), ~ifelse(. %like% "2/2" | . %like% "2/0" | . %like% "0/2" , "2", .)) %>% 
        mutate_at(c(10:ncol(.)), ~ifelse(. %like% "3/3" | . %like% "3/0" | . %like% "0/3" , "3", .)) %>% 
        mutate_at(c(10:ncol(.)), ~ifelse(. %like% "4/4" | . %like% "4/0" | . %like% "0/4" , "4", .)) %>% 
        mutate_at(c(10:ncol(.)), ~ifelse(. %like% "0/0" , "0", .)) %>% # 0 for reference allele
        mutate_at(c(10:ncol(.)), ~ifelse(. %like% "./." , "-1", .)) %>% # for missing 
        select(-c(3:9)) %>%
        arrange(CHROM, POS) %>%
        write_tsv("hmmIBD/sample_miss_5/hmmIBD.tsv")   