# Packages
library(dplyr)
library(readr)
library(stringr)
library(data.table)
library(tibble)
library(janitor)

# Define Edwins Function
.add.label.to.IBD <-
  function(IBD, metadata, label.col, sample.col = "Sample") {
    metadata <- metadata[, c(sample.col, label.col)]
    
    IBD <- merge(
      IBD,
      metadata,
      by.x = "sample1",
      by.y = sample.col,
      all.x = TRUE,
      sort = FALSE
    )
    colnames(IBD)[length(IBD)] <- paste0(label.col, 1)
    
    IBD <- merge(
      IBD,
      metadata,
      by.x = "sample2",
      by.y = sample.col,
      all.x = TRUE,
      sort = FALSE
    )
    colnames(IBD)[length(IBD)] <- paste0(label.col, 2)
    
    IBD
  }


.create.pairwise.matrix <- function(metadata, label.col) {
  unique.labels <- sort(unique(metadata[, label.col]))
  labels.length <- length(unique.labels)
  matrix(
    nrow = labels.length,
    ncol = labels.length,
    dimnames = list(unique.labels, unique.labels)
  )
}


.populate.pairwise.matrix.with.statistic <-
  function(pairwise.matrix,
           IBD,
           label.col,
           statistic) {
    label.cols1 <- paste0(label.col, 1)
    label.cols2 <- paste0(label.col, 2)
    labels <- rownames(pairwise.matrix)
    
    for (i in 1:nrow(pairwise.matrix)) {
      label1 <- labels[i]
      
      for (j in 1:ncol(pairwise.matrix)) {
        if (i < j)
          break
        label2 <- labels[j]
        
        e11 <- IBD[, label.cols1] == label1
        e12 <- IBD[, label.cols1] == label2
        e21 <- IBD[, label.cols2] == label1
        e22 <- IBD[, label.cols2] == label2
        
        pairwise.matrix[i, j] <-
          statistic(IBD[(e11 & 
                            e22) | (e12 & e21), "fract_sites_IBD"]) 
      }
    }
    
    pairwise.matrix
  }


calculate.pairwise.matrix <-
  function(IBD,
           metadata,
           label.col,
           statistic,
           sample.col = "Sample") {
    IBD <-
      .add.label.to.IBD(IBD, metadata, label.col, sample.col = sample.col)
    
    pairwise.matrix <-
      .create.pairwise.matrix(IBD, paste0(label.col, 1))
    
    .populate.pairwise.matrix.with.statistic(pairwise.matrix, IBD, label.col, statistic)
  }


save.pairwise.matrix <- function(pairwise.matrix, file) {
  write.table(pairwise.matrix,
              file,
              quote = FALSE,
              sep = "\t",
              na = "")
}


.add.sentinel.to.upper.triangle <- function(pairwise.matrix) {
  random <- rnorm(1)
  while (any(pairwise.matrix == random, na.rm = TRUE))
    random <- rnorm(1)
  
  pairwise.matrix[upper.tri(pairwise.matrix)] <- random
  
  pairwise.matrix
}


.transform.pairwise.matrix.to.table <- function(pairwise.matrix) {
  pairwise.matrix <- .add.sentinel.to.upper.triangle(pairwise.matrix)
  
  randoms <- pairwise.matrix[upper.tri(pairwise.matrix)]
  random <- randoms[1]
  if (any(is.na(randoms)) ||
      !all(randoms == random))
    warning("Matrix upper triangle does not contain sentinels.")
  
  pairwise.table <- as.data.frame(as.table(pairwise.matrix))
  pairwise.table[!(pairwise.table[, "Freq"] %in% random), ]
}


.concatenate.pairwise.labels <-
  function(pairwise.table, sep = "-") {
    labels <-
      do.call(paste, c(pairwise.table[, c("Var1", "Var2")], sep = sep))
    as.data.frame(cbind(labels, pairwise.table[, "Freq"]))
  }


calculate.pairwise.table <- function(IBD,
                                     metadata,
                                     label.col,
                                     statistics,
                                     sample.col = "Sample") {
  pairwise.tables <- NULL
  # safeguard against non-vector
  statistics <- c(statistics)
  stats <- names(statistics)
  
  
  for (i in seq_along(statistics)) {
    pairwise.matrix <- calculate.pairwise.matrix(IBD,
                                                 metadata,
                                                 label.col,
                                                 statistics[[i]],
                                                 sample.col = sample.col)
    pairwise.table <-
      .transform.pairwise.matrix.to.table(pairwise.matrix)
    
    pairwise.table <- .concatenate.pairwise.labels(pairwise.table)
    if (is.null(stats)) {
      names(pairwise.table) <- c(label.col, LETTERS[i])
    } else {
      names(pairwise.table) <- c(label.col, stats[i])
    }
    
    if (is.null(pairwise.tables)) {
      pairwise.tables <- pairwise.table
      next
    }
    
    pairwise.tables <-
      merge(pairwise.tables,
            pairwise.table,
            by = label.col,
            all = TRUE)
  }
  
  pairwise.tables
}


save.pairwise.table <- function(pairwise.table, file) {
  write.table(
    pairwise.table,
    file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
}


statistics <- c(min, function(x)
  quantile(x, probs = 0.25),
  median, mean, function(x)
    quantile(x, probs = 0.75), max)

names(statistics) <- c("Min.", "1st Qu.", "Median", "Mean",
                       "3rd Qu.", "Max.")

# Metadata
#Load packages
library(tidyverse)
library(ape)
library(ggtree)
library(viridis)
# Metadata

## Read in data

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

# Read in cluster data
metadata <- read_tsv("/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/outputs/05_Analyses/Population_genetics/outputs/admix_clusters.tsv")


# stat medians

Pk.IBD <- read_table("Pk.hmm_fract.txt") 

Pk.IBD.file <- "Pk.hmm_fract.txt"
Pk.IBD.tmp <- read.delim(Pk.IBD.file) 

Pk.pairwise.median <- calculate.pairwise.matrix(Pk.IBD, metadata, "Cluster", median)
save.pairwise.matrix(Pk.pairwise.median, "Pk_pairwise_IBD_median.txt")

## stat summary
Pk.pairwise.summary <- calculate.pairwise.table(Pk.IBD, metadata, "Cluster", statistics)
save.pairwise.table(Pk.pairwise.summary, "Pk_pairwise_summary.txt")


