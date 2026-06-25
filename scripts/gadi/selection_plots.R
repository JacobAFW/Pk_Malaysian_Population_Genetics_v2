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
    theme(legend.position = "none") +
    scale_colour_manual(values = c("#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF")) +
    scale_x_continuous(breaks = x_axis$POS, labels = x_axis$CHR) +
    xlab("Position") +
    ylab(MODEL) +
    geom_hline(yintercept = c(as.numeric(hline[1,2]), as.numeric(hline[3,2])), linetype="dotted")

  ggsave(paste0(MODEL, ".png"), dpi = 300, width = 16, selection_plot)

  selection_plot <- plot_data %>%
      ggplot(aes(x = POS, y = LOGPVALUE, colour = CHR)) +
      geom_jitter() + 
      theme(legend.position = "none") +
      scale_x_continuous(breaks = x_axis$POS, labels = x_axis$CHR) +
      scale_colour_manual(values = c("#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF", "#39568CFF", "#29AF7FFF")) +
      xlab("Chromosomes") +
      ylab(expression(MODEL ~ ~-log[10](italic(p))))

  ggsave(paste0(MODEL, "_pvalue.png"), dpi = 300, width = 16, selection_plot)
}
