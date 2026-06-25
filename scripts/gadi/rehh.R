# Packages
library(rehh)
library(tidyverse)
library(data.table)
library(R.utils)

# WHOLE GENOME iHS

## need to read in chr files independently and build up an object resulting from the scan_hh outputs
### first create an empty df to be appended in the loop and a list of chr names
### the loop then creates the vcf file name (hap_file), and reads in the file to hh, which is then parsed into the scan_hh function to produce an object i that is appended to the growing wgscan object
### the resulting wgscan produces the same results as the scan_hh object above
### this can then be parsed into ihh2ihs

chromsomes_files <- c("ordered_PKNH_01_v2", "ordered_PKNH_02_v2", "ordered_PKNH_03_v2", "ordered_PKNH_04_v2", "ordered_PKNH_05_v2",
    "ordered_PKNH_06_v2", "ordered_PKNH_07_v2", "ordered_PKNH_08_v2", "ordered_PKNH_09_v2", "ordered_PKNH_10_v2", "ordered_PKNH_11_v2",
    "ordered_PKNH_12_v2", "ordered_PKNH_13_v2", "ordered_PKNH_14_v2")

wgscan <- tibble("CHR" = NA, "POSITION" = NA, "FREQ_A" = NA, "FREQ_D" = NA, "NHAPLO_A" = NA, "NHAPLO_D" = NA, "IHH_A" = NA, "IHH_D" = NA, "IES" = NA, "INES" = NA)

for(i in chromsomes_files) {
    # haplotype file name for each chromosome
    hap_file = paste(i, ".vcf.gz", sep = "") # filename pattern
    # create internal representation
    hh <- data2haplohh(hap_file = hap_file, 
                    chr.name = i, 
                    polarize_vcf = FALSE, # if the AA key is absent
                    min_perc_geno.mrk = 100, # discard markers genotyped on < 100% of haplotypes
                    min_maf = 0.05, # discard markers with a minor allele frequency of < 0.05)
                    vcf_reader = "data.table") # use this package
    scan <- scan_hh(hh) # perform scan on a single chromosome (calculate iHH/iES values)
    wgscan <- wgscan %>% # append the wgscan object every round with additional data
    rbind(scan)
}

wgs_ihs <- wgscan %>%
    na.omit() %>%
    ihh2ihs(., min_maf = 0.05, # default
            freqbin = 0.025 # default
            )

### Output is a list with two elements:
### ihs: a data frame with iHS and corresponding p-value piHS (p-value in a −log10 scale assuming the iHS are normally distributed under the neutral hypothesis)
### frequency.class: a data frame summarizing the derived allele frequency bins used for standardization and mean and standard deviation of the un-standardized values
wgs_ihs$ihs
wgs_ihs$frequency.class


# Scan for candidate regions 
## First update iHS table so that its compatible
wgs_ihs_2 <- wgs_ihs
wgs_ihs_2$ihs <- wgs_ihs_2$ihs %>%
    mutate(CHR = str_remove(CHR, "ordered_PKNH_")) %>%
    mutate(CHR = str_remove(CHR, "_v2")) %>%
    mutate(CHR = as.numeric(CHR)) %>%
    as.data.frame()

candidate_regions <- calc_candidate_regions(wgs_ihs_2,
                                 threshold = 4,
                                 pval = TRUE,
                                 window_size = 10000,
                                 overlap = 1000,
                                 min_n_extr_mrk = 3) %>%
                                 add_column(Stat = "iHS") 

## Save candidate regions for all stats
candidate_regions %>%
  write_tsv("candidate_regions_iHS.tsv")

selection_plots <- function(DATA, MODEL){
  plot_data <- DATA %>%
    na.omit() %>%
    mutate(CHR = str_remove(CHR, "ordered_PKNH_")) %>%
    mutate(CHR = str_remove(CHR, "_v2")) %>%
    mutate(CHR = as.numeric(CHR)) %>%
    arrange(CHR, POSITION) %>%
    mutate(POS = 1:nrow(.)) %>%
    mutate(CHR = as.factor(CHR))

  x_axis <- plot_data %>%
    group_by(CHR) %>%
    summarise(POS = median(POS)) %>%
    mutate(CHR = str_remove(CHR, "^0"))

  hline <- plot_data %>% 
    summarise(enframe(quantile(!!sym(MODEL), c(0.001, 0.5, 0.999)), "quantile", MODEL))

  selection_plot <- plot_data %>%
    ggplot(aes(x = POS, y = !!sym(MODEL), colour = CHR)) +
    geom_point() + 
    theme(legend.position = "none", axis.text=element_text(size=14), axis.title=element_text(size=16)) +
    scale_colour_manual(values = c("#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF")) +
    scale_x_continuous(breaks = x_axis$POS, labels = x_axis$CHR) +
    xlab("Chromosome") +
    ylab("iHS") +
    geom_hline(yintercept = c(as.numeric(hline[1,2]), as.numeric(hline[3,2])), linetype="dotted")

  ggsave(paste0(MODEL, ".png"), dpi = 300, width = 16, selection_plot)

  selection_plot <- plot_data %>%
      ggplot(aes(x = POS, y = LOGPVALUE, colour = CHR)) +
      geom_jitter() + 
      theme(legend.position = "none", axis.text=element_text(size=14), axis.title=element_text(size=16)) +
      scale_x_continuous(breaks = x_axis$POS, labels = x_axis$CHR) +
      scale_colour_manual(values = c("#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF")) +
      xlab("Chromosomes") +
      ylab(expression(iHS ~ ~-log[10](italic(p))))

  ggsave(paste0(MODEL, "_pvalue.png"), dpi = 300, width = 16, selection_plot)
}

selection_plots(wgs_ihs$ihs, "IHS")
