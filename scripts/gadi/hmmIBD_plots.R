library(igraph)

.reorder.metadata <-
  function(IBD, metadata, sample.col = "Sample") {
    samples <- unique(c(IBD[, "sample1"], IBD[, "sample2"]))
    sample.order <- match(metadata[, sample.col], samples)
    actual.length <- length(na.omit(sample.order))
    
    metadata <- metadata[order(sample.order), ]
    metadata[1:actual.length, ]
  }


.get.labels <-
  function(metadata,
           label.col,
           label.palette = NULL) {
    labels <- metadata[, label.col]
    labels[is.na(labels)] <- NaN
    
    labels <- factor(labels)
    
    if (is.null(label.palette)) {
      labels
    } else {
      factor(labels, levels = names(label.palette))
    }
  }


.generate.label.palette <-
  function(labels) {
    label.names <- levels(labels)
    label.palette <-  #c("#482677FF", "#94D840FF", "#3F4788FF", "#29AF7FFF", "#32648EFF") # for Clusters
       viridis(length(label.names), option = ifelse(length(label.names) < 6, "D", "H")) 
    names(label.palette) <- label.names
    
    label.palette
  }


.generate.label.colours <-
  function(labels, label.palette)
    label.palette[labels]


.create.edgelist <-
  function(IBD)
    IBD[, c("sample1", "sample2", "fract_sites_IBD")]


.create.vertices <-
  function(metadata, sample.col = "Sample") {
    vertices <- data.frame(metadata[, sample.col])
    names(vertices) <- sample.col
    
    vertices
  }


.internal.plot.IBD <- function(IBD.graph,
                               unlabelled = TRUE,
                               coords = layout_nicely,
                               percent.cutoff = NA) {
  if (unlabelled) {
    plot(
      IBD.graph,
      layout = coords,
      vertex.size = 4,
      vertex.label = NA,
      main = paste0("IBD >=", percent.cutoff, "%")
    )
  } else {
    plot(
      IBD.graph,
      layout = coords,
      vertex.size = 4,
      vertex.label.cex = 0.3,
      main = paste0("IBD >=", percent.cutoff, "%")
    )
  }
}


.plot.IBD <-
  function(file,
           IBD.graph,
           unlabelled = TRUE,
           label.palette = NULL,
           labels = NULL,
           legend.title = NA,
           coords = layout_nicely,
           percent.cutoff = NA) {
    pdf(file)
    par(mar = rep.int(1, 4) + 0.1)
    
    .internal.plot.IBD(
      IBD.graph,
      unlabelled = unlabelled,
      coords = layout_nicely,
      percent.cutoff = percent.cutoff
    )
    
    if (!is.null(label.palette))
      legend(
        "bottomleft",
        legend = levels(labels),
        fill = label.palette,
        title = legend.title,
        cex = 0.6
      )
    
    dev.off()
  }


.plot.IBDs <- function(prefix,
                       IBD.graph,
                       label.palette,
                       labels,
                       legend.title,
                       coords,
                       percent.cutoff) {
  labelled.file <-
    paste0(prefix, "_labelled_IBD", percent.cutoff, ".pdf")
  .plot.IBD(
    labelled.file,
    IBD.graph,
    unlabelled = FALSE,
    label.palette = label.palette,
    labels = labels,
    legend.title = legend.title,
    coords = coords,
    percent.cutoff = percent.cutoff
  )
  
  unlabelled.file <-
    paste0(prefix, "_unlabelled_IBD", percent.cutoff, ".pdf")
  .plot.IBD(
    unlabelled.file,
    IBD.graph,
    unlabelled = TRUE,
    label.palette = label.palette,
    labels = labels,
    legend.title = legend.title,
    coords = coords,
    percent.cutoff = percent.cutoff
  )
}


plot.IBD <- function(file,
                     IBD.cutoffs,
                     IBD,
                     metadata,
                     label.col,
                     unlabelled = TRUE,
                     label.palette = NULL,
                     legend.title = NA,
                     sample.col = "Sample") {
  edgelist <- .create.edgelist(IBD)
  
  metadata <- .reorder.metadata(IBD, metadata)
  
  vertices <- .create.vertices(metadata)
  
  sqrt.n.plots <- ceiling(sqrt(length(IBD.cutoffs)))
  size <- 7 * sqrt.n.plots
  
  legend.position <- sqrt.n.plots ^ 2 - sqrt.n.plots + 1
  legend.plot <- FALSE
  
  pdf(file, width = size, height = size)
  par(mar = rep.int(1, 4) + 0.1,
      mfrow = c(sqrt.n.plots, sqrt.n.plots))
  
  if (is.null(label.palette)) {
    labels <- .get.labels(metadata, label.col)
    label.palette <- .generate.label.palette(labels)
    
  } else {
    labels <-
      .get.labels(metadata, label.col, label.palette = label.palette)
  }
  
  label.colours <-
    .generate.label.colours(labels, label.palette)
  
  for (i in seq_along(IBD.cutoffs)) {
    d <- edgelist[edgelist[, "fract_sites_IBD"] >= IBD.cutoffs[i], ]
    
    IBD.graph <-
      graph_from_data_frame(d, directed = FALSE, vertices = vertices)
    
    coords <- layout_(IBD.graph, nicely())
    percent.cutoff <- round(IBD.cutoffs[i] * 100, digits = 2)
    
    IBD.graph <-
      set_vertex_attr(IBD.graph, "color", value = label.colours)
    
    .internal.plot.IBD(IBD.graph,
                       unlabelled,
                       coords,
                       percent.cutoff)
    
    if (i == legend.position) {
      legend(
        "bottomleft",
        legend = levels(labels),
        fill = label.palette,
        title = legend.title,
        cex = 0.6
      )
      
      legend.plot <- TRUE
    }
  }
  
  if (!legend.plot) {
    while (i != legend.position) {
      plot.new()
      i <- i + 1
    }
    
    legend(
      "bottomleft",
      legend = levels(labels),
      fill = label.palette,
      title = legend.title,
      cex = 0.6
    )
  }
  
  dev.off()
}


plot.IBDs <- function(prefixes,
                      IBD.cutoffs,
                      IBD,
                      metadata,
                      label.cols,
                      label.palettes = NULL,
                      legend.titles = NA,
                      sample.col = "Sample") {
  edgelist <- .create.edgelist(IBD)
  
  metadata <- .reorder.metadata(IBD, metadata)
  
  vertices <- .create.vertices(metadata)
  
  for (IBD.cutoff in IBD.cutoffs) {
    d <- edgelist[edgelist[, "fract_sites_IBD"] >= IBD.cutoff, ]
    
    IBD.graph <-
      graph_from_data_frame(d, directed = FALSE, vertices = vertices)
    
    coords <- layout_(IBD.graph, nicely())
    percent.cutoff <- round(IBD.cutoff * 100, digits = 2)
    
    for (label.col in label.cols) {
      label.palette <- label.palettes[[label.col]]
      
      if (is.null(label.palette)) {
        labels <- .get.labels(metadata, label.col)
        label.palette <- .generate.label.palette(labels)
        
      } else {
        labels <-
          .get.labels(metadata, label.col, label.palette = label.palette)
      }
      
      label.colours <-
        .generate.label.colours(labels, label.palette)
      
      IBD.graph <-
        set_vertex_attr(IBD.graph, "color", value = label.colours)
      
      legend.title <- legend.titles[[label.col]]
      
      prefix <- prefixes[[label.col]]
      
      .plot.IBDs(prefix,
                 IBD.graph,
                 label.palette,
                 labels,
                 legend.title,
                 coords,
                 percent.cutoff)
    }
  }
}


read.pairwise.matrix <- function(file) {
  as.matrix(read.delim(file, check.names = FALSE))
}


IBD.cutoffs <- c(0.125, 0.25, 0.5)

# include +-5% IBD
relaxed.IBD.cutoffs <- c(IBD.cutoffs, 1)
relaxed.IBD.cutoffs <- relaxed.IBD.cutoffs * 0.95


# Plotting
library(igraph)
library(tidyverse)
#library(MetamapsDB)
library(viridis)
library(janitor)

# Metadata

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


## Read in cluster data (based on ADMIXTURE)
metadata <- read_tsv("/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/outputs/05_Analyses/Population_genetics/outputs/admix_clusters.tsv")


Pk.IBD.file <- "Pk.hmm_fract.txt"
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
  
## Fraction IBD within and between Clusters
IBD_fract_plot <- IBD_meta_combined %>%
  unite("Clusters", c("Cluster.x", "Cluster.y"), sep = "-") %>% 
  ggplot(aes(x = Clusters, y = fract_sites_IBD, colour = Clusters)) +
  geom_boxplot() +
  theme(legend.position = "none", 
    legend.title = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  coord_flip() +
  scale_color_viridis_d() +
  ylab("Fraction of IBD sites")

ggsave("fract_IBD_clusters.png", dpi = 300, height = 7, width = 7, IBD_fract_plot)

# Base plots
Pk.label.cols <- c("Cluster")
Pk.legend.titles <- list("Cluster" = "Cluster")
Pk.prefixes <- metadata %>% 
  select(Cluster) %>%
  unique() %>% 
  as.list()

plot.IBDs(
  Pk.prefixes,
  IBD.cutoffs,
  Pk.IBD,
  as.data.frame(metadata), # doesn't like it as a tibble!
  Pk.label.cols,
  legend.titles = Pk.legend.titles
)


Pk.median.file <- "Pk_pairwise_IBD_median.txt"
Pk.median <- read.pairwise.matrix(Pk.median.file)

Pk.IBD.cutoffs <- c(fivenum(Pk.median), relaxed.IBD.cutoffs)

IBDs.plot.file <- "Pk_major_IBDs.pdf"

plot.IBD(
  IBDs.plot.file,
  Pk.IBD.cutoffs,
  Pk.IBD,
  as.data.frame(metadata),
  "Cluster",
  legend.title = "Cluster"
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

clonal_samples("Mn", 0.95) %>% 
  rbind(
    clonal_samples("Mf", 0.95)
  ) %>%
  rbind(
    clonal_samples("Peninsular", 0.95)
  ) %>% 
  write_tsv("IBD_clonal_transmission.tsv")

# Investigate samples of high relatedness wihtin clusters





# Base plots - District

# Mn
Mn.IBD <- Pk.IBD %>%
    as.data.frame() %>%
    left_join(
      metadata %>% mutate(sample1 = Sample) %>% select(sample1, Cluster),
      by = "sample1"
    ) %>% 
    left_join(
      metadata %>% mutate(sample2 = Sample) %>% select(sample2, Cluster),
      by = "sample2" 
    ) %>% 
    filter(Cluster.x == "Mn" & Cluster.y == "Mn") %>% 
    dplyr::select(-starts_with("Cluster"))
  

IBD.cutoffs <- c(0.36, 0.365, 0.37, 0.375, 0.38, 0.385, 0.39, 0.395, 0.4)

outbreak.label.cols <- c("district")
outbreak.legend.titles <- list("district" = "district")

metadata_dsitrict <- metadata %>% 
  mutate(sample1 = Sample) %>%
  mutate(sample1 = str_remove(sample1, "_DK.*")) %>% 
  mutate(sample1 = str_replace(sample1, "_", "-")) %>% 
  left_join(
    metadata_full %>% 
      select(Sample, district) %>% 
      rename(sample1 = Sample)
  ) %>%
  select(-sample1)

outbreak.prefixes <- metadata_dsitrict %>% 
  select(district) %>%
  unique() %>% 
  as.list()

plot.IBDs(
  outbreak.prefixes,
  IBD.cutoffs,
  Mn.IBD,
  as.data.frame(metadata_dsitrict),
  outbreak.label.cols,
  legend.titles = outbreak.legend.titles
)


clonal_samples("Mn", .375) %>%
  write_tsv("Mn_clusters.tsv")

# Clusters/outbreaks
# Kapit outbreak

clonal_samples("Mn", .375) %>%
  #group_by(district1) %>% 
  #summarise(n())
  filter((district1 == "Kapit") & (district2 == "Kapit" | district2 == "Sarikei" | district2 == "Betong"))

# NOT our samples - no dates

#Mf

Mf.IBD <- Pk.IBD %>%
    as.data.frame() %>%
    left_join(
      metadata %>% mutate(sample1 = Sample) %>% select(sample1, Cluster),
      by = "sample1"
    ) %>% 
    left_join(
      metadata %>% mutate(sample2 = Sample) %>% select(sample2, Cluster),
      by = "sample2" 
    ) %>% 
    filter(Cluster.x == "Mf" & Cluster.y == "Mf") %>% 
    dplyr::select(-starts_with("Cluster"))
  

IBD.cutoffs <- c(0.05, 0.055, .06, .065, .07, 0.075)

outbreak.label.cols <- c("district")
outbreak.legend.titles <- list("district" = "district")

metadata_dsitrict <- metadata %>% 
  mutate(sample1 = Sample) %>%
  mutate(sample1 = str_remove(sample1, "_DK.*")) %>% 
  mutate(sample1 = str_replace(sample1, "_", "-")) %>% 
  left_join(
    metadata_full %>% 
      select(Sample, district) %>% 
      rename(sample1 = Sample)
  ) %>%
  select(-sample1)

outbreak.prefixes <- metadata_dsitrict %>% 
  select(district) %>%
  unique() %>% 
  as.list()

plot.IBDs(
  outbreak.prefixes,
  IBD.cutoffs,
  Mf.IBD,
  as.data.frame(metadata_dsitrict),
  outbreak.label.cols,
  legend.titles = outbreak.legend.titles
)


clonal_samples("Mf", .05) %>% head()
  write_tsv("Mf_clusters.tsv")



# Kudat outbreak
Kudat_outbreak <- clonal_samples("Mf", .05) %>%
  #group_by(district1) %>% 
  #summarise(n())
  filter((district1 == "Kudat" | district1 == "Tuaran" | district1 == "Kota Marudu" | district1 == "Beluran") & ((district2 == "Kudat" | district2 == "Tuaran" | district2 == "Kota Marudu" | district2 == "Beluran"))) %>% 
  filter(!grepl("REDACTED", sample1) & !grepl("REDACTED", sample2)) %>% 
  filter(!grepl("REDACTED", sample1) & !grepl("REDACTED", sample2)) %>%
  filter(!grepl("REDACTED", sample1) & !grepl("REDACTED", sample2)) %>%
  filter(!grepl("REDACTED", sample1) & !grepl("REDACTED", sample2)) %>%
  filter(!grepl("REDACTED", sample1) & !grepl("REDACTED", sample2)) %>%
  filter(!grepl("REDACTED", sample1) & !grepl("REDACTED", sample2)) %>%
  filter(!grepl("REDACTED", sample1) & !grepl("REDACTED", sample2)) %>%
  filter(!grepl("REDACTED", sample1) & !grepl("REDACTED", sample2)) %>%
  filter(!grepl("REDACTED", sample1) & !grepl("REDACTED", sample2)) 
   # need to make sure there are samples joining MANY samples to comfirm its the large cluster

cluster_plot <- Kudat_outbreak %>% 
  select(sample1, enroldate1, district1, village1) %>% 
  rename(sample = sample1, enroldate = enroldate1, district = district1, village = village1) %>%
  rbind(
  Kudat_outbreak %>% 
    select(sample2, enroldate2, district2, village2) %>%
    rename(sample = sample2, enroldate = enroldate2, district = district2, village = village2) 
  ) %>% 
  ggplot(aes(x = enroldate, fill = district)) +
  geom_histogram()  +
  scale_fill_manual(values = c("#888888", "#D55E00", "#88CCEE", "#117733"))

  
ggsave("Mf_plots/mf_IBD_cluster_kudat.png", dpi = 300, cluster_plot)
