# Stage B-7 plot: Fst manhattan (Mn vs Mf). Lifted from Analyses.Rmd L1483-1499.
suppressPackageStartupMessages({ library(tidyverse); library(here) })
source(here::here("scripts/_setup.R"))

Fst_matrix <- read_table(file.path(DATA_PROC, "fst", "Pk.fst"),
                         col_names = TRUE, show_col_types = FALSE) |>
  mutate(CHR = str_remove(SNP, "ordered_PKNH_"),
         CHR = str_remove(CHR, "_v2.*"),
         CHR = str_remove(CHR, "^0"))

plot_data <- Fst_matrix |>
  filter(FST != "nan", FST != "-nan") |>
  mutate(FST = as.numeric(FST)) |>
  filter(FST > 0, !is.na(FST)) |>
  arrange(CHR, POS) |>
  mutate(ROW = row_number())

x_axis <- plot_data |>
  group_by(CHR) |>
  summarise(ROW = median(ROW))

Fst_plot <- ggplot(plot_data, aes(x = ROW, y = FST, colour = CHR)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_colour_viridis_d("Chr") +
  scale_x_continuous(breaks = x_axis$ROW, labels = x_axis$CHR) +
  xlab("Chromosome") + ylab("Fst") +
  theme_pk() +
  theme(legend.position = "none")

save_fig(Fst_plot, "stage_b/15_fst_manhattan_mf_vs_mn", width = 12, height = 4)

mean_fst <- mean(plot_data$FST, na.rm = TRUE)
write_tsv(tibble(mean_fst = mean_fst),
          file.path(RES_DIR, "stage_b", "15_fst_mean_mf_vs_mn.tsv"))
message(sprintf("[15b] regen mean_fst (positive only) = %.6f  vs truth 0.15411", mean_fst))
message(sprintf("[15b] regen mean_fst (all valid, incl 0 & neg) = %.6f",
                mean(Fst_matrix$FST |> as.numeric(), na.rm = TRUE)))
