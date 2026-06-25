library(tidyverse)

baf_df <- read_tsv("moimix/BAF_dataframe.tsv") %>% # Generated with previous script
  mutate(chromosome = str_remove(chromosome, "ordered_PKNH_")) %>% 
  mutate(chromosome = str_remove(chromosome, "_v2")) %>%
  mutate(chromosome = as.factor(as.numeric(chromosome))) %>%
  rename_with(~str_remove(., "-"))

## Plot BAF for different Fws samples with ggplot
plot_baf <- function(DATA, SAMPLE){
DATA %>%
    ggplot(aes_string(x = "variant.id", y = SAMPLE, colour = "chromosome")) + # use aes string so that we can paste the str in from the function arguments - passing the sample in not as a string does not work
    geom_point() +
    scale_colour_manual(values = c(rep_len(c("#404788FF", "#1F968BFF"), length(unique(baf_df$chromosome))))) +
    xlab("Chromosome") +
    ylab("NRAF") + 
    theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(nrow = 1))
}
# top 20
sample_list <- read_tsv("moimix/fws_MOI.tsv") %>%
  filter(Proportion < 0.95) %>%
  filter(!grepl("ERR", sample)) %>%
  arrange(Proportion) %>%
  slice(1:20) %>%
  select(sample) %>%
  mutate(sample = str_remove(sample, "-")) %>% 
  pull()

for (i in sample_list){
  baf_plot <- plot_baf(baf_df, paste0(i))
  ggsave(paste0("moimix/",i,"_NRAF.png"), width = 14, dpi = 300, baf_plot)
}


