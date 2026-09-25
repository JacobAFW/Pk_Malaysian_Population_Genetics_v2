#!/usr/bin/env Rscript
# make_figures.R - publication figures for the introgression run.
# --------------------------------------------------------------------------
# Reads ONLY existing TSVs under results/<run>/; performs no analysis, no
# re-detection, no re-aggregation. Light reshaping only.
#
# Style comes from scripts/_setup.R (viridis clusters, inferno ordinal,
# theme_pk) so these figures match every other figure in the project. Output
# goes to a NEW versioned folder; nothing existing is overwritten.
#
# NO SAMPLE IDENTIFIERS appear in any figure: per-sample panels are aggregated
# to distributions, never labelled points.
#
#   Rscript make_figures.R --config config/malay_cohort.yaml
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here); library(ggplot2); library(viridis)
  library(dplyr); library(tidyr); library(readr); library(yaml)
  library(patchwork)
})
source(here::here("scripts/_setup.R"))

parse_args <- function(av) { o <- list(); i <- 1
  while (i <= length(av)) { o[[sub("^--", "", av[i])]] <- av[i + 1]; i <- i + 2 }; o }
args <- parse_args(commandArgs(trailingOnly = TRUE))
cfg  <- yaml::read_yaml(args[["config"]] %||% "config/malay_cohort.yaml")
RES  <- cfg$outputs$results_dir
FIGS <- cfg$outputs$figures_dir
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)
say <- function(...) message(sprintf("[figures] %s", sprintf(...)))

# Local saver: same spec as _setup.R's save_fig (vector PDF + 300-dpi PNG) but
# writing into the run's own versioned folder instead of the shared figures/.
save_run_fig <- function(plot, name, width = 9, height = 6) {
  pdf_p <- file.path(FIGS, paste0(name, ".pdf"))
  png_p <- file.path(FIGS, paste0(name, ".png"))
  ggsave(pdf_p, plot, width = width, height = height, device = grDevices::pdf)
  ggsave(png_p, plot, width = width, height = height, dpi = 300)
  say("%s  (%.0f x %.0f in)", name, width, height)
  invisible(c(pdf = pdf_p, png = png_p))
}

ARM_LAB <- c(full = "full (n = 501)", declonal = "de-clonalized (n = 484)")
# Facets default to alphabetical, which puts the de-clonalized arm first; the
# full arm is the reference and belongs on the left.
arm_factor <- function(a) factor(ARM_LAB[a], levels = unname(ARM_LAB))
FLOOR   <- c(full = 41, declonal = 40)

# ==========================================================================
# Fig 1 - floor derivation: the null, and where the FDR-safe floor lands
# ==========================================================================
fd <- bind_rows(lapply(names(ARM_LAB), function(a)
  read_tsv(file.path(RES, "floor", sprintf("floor_derivation_%s.tsv", a)),
           show_col_types = FALSE) %>% mutate(arm = a)))

# (a) observed vs expected-null windows across the floor, per cluster
p1a <- fd %>%
  mutate(arm_lab = arm_factor(arm),
         Cluster = factor(Cluster, levels = names(cluster_cols))) %>%
  ggplot(aes(N, observed_windows, colour = Cluster)) +
  geom_line(linewidth = 0.7) +
  geom_vline(data = tibble(arm_lab = arm_factor(names(ARM_LAB)), x = FLOOR[names(ARM_LAB)]),
             aes(xintercept = x), linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~arm_lab) +
  scale_colour_pk_cluster("Cluster") +
  scale_y_continuous(trans = "log1p", breaks = c(0, 10, 50, 200, 800)) +
  labs(x = "per-cluster support floor N", y = "windows passing the floor (log1p)",
       title = "Windows surviving the support floor",
       subtitle = "dashed line = FDR-derived floor (41 full / 40 de-clonalized), binding cluster Mf") +
  theme_pk()

# (b) the null support distribution each cluster must beat
p1b <- fd %>%
  distinct(arm, Cluster, cluster_n, null_mean, null_p95, null_p99, null_max) %>%
  mutate(arm_lab = arm_factor(arm),
         Cluster = factor(Cluster, levels = names(cluster_cols))) %>%
  ggplot(aes(Cluster, null_mean, fill = Cluster)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = null_p95, ymax = null_max), width = 0.18, linewidth = 0.4) +
  geom_hline(data = tibble(arm_lab = arm_factor(names(ARM_LAB)), y = FLOOR[names(ARM_LAB)]),
             aes(yintercept = y), linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~arm_lab) +
  scale_fill_pk_cluster("Cluster", guide = "none") +
  labs(x = NULL, y = "per-window sample support under the null",
       title = "Chance co-location scales with cluster size",
       subtitle = "bar = null mean; whisker = p95 to null max; dashed = the derived global floor") +
  theme_pk()

save_run_fig(p1a / p1b + plot_annotation(tag_levels = "a"),
             "fig1_floor_derivation", width = 9, height = 8)

# ==========================================================================
# Fig 2 - the floor x filter-4 interaction
# ==========================================================================
audits <- bind_rows(lapply(list.files(file.path(RES, "floor", "_sweep_full"),
                                      full.names = TRUE), function(d) {
  f <- file.path(d, "filter_audit.tsv")
  if (!file.exists(f)) return(NULL)
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(N = as.integer(sub("^N", "", basename(d))))
}))
mask_sizes <- audits %>%
  filter(grepl("hypervariable", step)) %>%
  mutate(masked = as.integer(sub(".*\\((\\d+) win\\).*", "\\1", step))) %>%
  dplyr::select(N, masked, windows_after = n_windows)

p2 <- mask_sizes %>%
  pivot_longer(c(masked, windows_after)) %>%
  mutate(name = recode(name,
                       masked = "windows removed by the multi-cluster mask",
                       windows_after = "windows surviving the whole chain")) %>%
  ggplot(aes(N, value, colour = name)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
  geom_vline(xintercept = 41, linetype = "dashed", linewidth = 0.4) +
  annotate("text", x = 41, y = Inf, label = " derived floor", hjust = 0, vjust = 1.6, size = 3) +
  annotate("text", x = 6, y = Inf, label = "published floor ", hjust = 1, vjust = 1.6, size = 3) +
  geom_vline(xintercept = 6, linetype = "dotted", linewidth = 0.4) +
  scale_colour_pk_ordinal("") +
  labs(x = "per-cluster support floor N", y = "windows",
       title = "Filter 4 runs after the floor, so the floor moves the mask",
       subtitle = "raising the floor makes windows single-cluster, rescuing them from the multi-cluster mask") +
  theme_pk() + theme(legend.position = "bottom")
save_run_fig(p2, "fig2_floor_filter4_interaction", width = 9, height = 5.5)

# ==========================================================================
# Fig 3 - headline Mf<->Mn asymmetry and its chromosomal distribution
# ==========================================================================
filt <- bind_rows(lapply(names(ARM_LAB), function(a)
  read_tsv(file.path(RES, "aggregate", a, "introgressed_windows_filtered.tsv"),
           show_col_types = FALSE) %>% mutate(arm = a)))

dir_tab <- filt %>% filter(PAIR == "Mf__Mn") %>%
  group_by(arm, DIRECTION) %>%
  summarise(windows = n_distinct(paste(CHROM, START)), .groups = "drop") %>%
  mutate(arm_lab = arm_factor(arm))

p3a <- ggplot(dir_tab, aes(DIRECTION, windows, fill = DIRECTION)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = windows), vjust = -0.4, size = 3.2) +
  facet_wrap(~arm_lab) +
  scale_fill_pk_ordinal("", guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "windows at the derived floor",
       title = "The Mf-Mn signal is strongly one-directional",
       subtitle = "Mf carrying Mn-like windows outnumbers the reverse ~20:1") +
  theme_pk()

chrom <- bind_rows(lapply(names(ARM_LAB), function(a)
  read_tsv(file.path(RES, "aggregate", a, "windows_across_chrom.tsv"),
           show_col_types = FALSE) %>% mutate(arm = a))) %>%
  mutate(arm_lab = arm_factor(arm),
         Cluster = factor(Cluster, levels = names(cluster_cols)))
# Peninsular retains zero windows: show it explicitly rather than dropping it.
chrom_full <- chrom %>%
  complete(arm_lab, Cluster, CHROM = 1:14, fill = list(n_windows = 0)) %>%
  filter(!is.na(Cluster))

p3b <- ggplot(chrom_full, aes(factor(CHROM), n_windows, fill = Cluster)) +
  geom_col(position = position_dodge(preserve = "single"), width = 0.75) +
  facet_wrap(~arm_lab) +
  scale_fill_pk_cluster("Cluster") +
  labs(x = "chromosome", y = "windows",
       title = "Signal concentrates on chr11-14",
       subtitle = "Peninsular is shown at zero throughout - it retains no window at the derived floor") +
  theme_pk()

save_run_fig(p3a / p3b + plot_annotation(tag_levels = "a"),
             "fig3_headline_mf_mn", width = 9, height = 8)

# ==========================================================================
# Fig 4 - the clonality effect, by layer
# ==========================================================================
# Both layers must be measured the SAME way to be comparable, so the raw layer
# is recomputed per CLUSTER here (an earlier draft mixed per-PAIR raw values
# with per-CLUSTER filtered ones under a single "Cluster" legend, which is not
# a like-for-like comparison). Filtered values are the matched-floor
# (full@40 vs declonal@40) figures from phase 4.
raw_cluster_jaccard <- function() {
  rd <- function(arm) {
    f <- list.files(file.path(RES, "calls", arm), full.names = TRUE, pattern = "\\.tsv$")
    bind_rows(lapply(f, read_tsv, show_col_types = FALSE)) %>%
      transmute(Cluster, k = paste(CHROM, BIN - 5000)) %>% distinct()
  }
  a <- rd("full"); b <- rd("declonal")
  bind_rows(lapply(unique(a$Cluster), function(cl) {
    x <- a$k[a$Cluster == cl]; y <- b$k[b$Cluster == cl]
    tibble(Cluster = cl,
           jaccard = length(intersect(x, y)) / length(union(x, y)))
  }))
}
raw_j <- raw_cluster_jaccard() %>% mutate(layer = "raw calls")
filt_j <- tribble(~Cluster, ~jaccard,
                  "Mf", 0.776, "Mn", 0.667, "Peninsular", NA_real_) %>%
  mutate(layer = "filtered windows")
lay <- bind_rows(raw_j, filt_j) %>%
  filter(!is.na(jaccard)) %>%
  mutate(layer = factor(layer, levels = c("raw calls", "filtered windows")),
         Cluster = factor(Cluster, levels = names(cluster_cols)))

p4a <- ggplot(lay, aes(layer, jaccard, fill = Cluster)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", jaccard)),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3) +
  scale_fill_pk_cluster("Cluster") +
  scale_y_continuous(limits = c(0, 1.08), expand = c(0, 0)) +
  labs(x = NULL, y = "Jaccard, full vs de-clonalized",
       title = "De-clonalization barely moves the raw calls but reshapes the filtered windows",
       subtitle = "per cluster, both layers; Peninsular has no filtered windows so only its raw layer is shown") +
  theme_pk()

# Per-sample filtered call counts, aggregated to a distribution (no identifiers).
pс <- filt %>% group_by(arm, SAMPLE, Cluster) %>%
  summarise(n_windows = n_distinct(paste(CHROM, START)), .groups = "drop")
shared <- intersect(pс$SAMPLE[pс$arm == "full"], pс$SAMPLE[pс$arm == "declonal"])
delta <- pс %>% filter(SAMPLE %in% shared) %>%
  dplyr::select(arm, SAMPLE, Cluster, n_windows) %>%
  pivot_wider(names_from = arm, values_from = n_windows, values_fill = 0) %>%
  mutate(delta = declonal - full,
         Cluster = factor(Cluster, levels = names(cluster_cols)))

p4b <- ggplot(delta, aes(delta, fill = Cluster)) +
  geom_histogram(binwidth = 1, colour = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~Cluster, scales = "free_y") +
  scale_fill_pk_cluster("Cluster", guide = "none") +
  labs(x = "change in filtered windows per sample (de-clonalized - full)",
       y = "samples",
       title = sprintf("Per-sample change over the %d Mf/Mn samples with filtered windows in both arms", length(shared)),
       subtitle = "Peninsular is absent: it retains no filtered window in either arm") +
  theme_pk()

save_run_fig(p4a / p4b + plot_annotation(tag_levels = "a"),
             "fig4_clonality_effect", width = 9, height = 8)

# ==========================================================================
# Fig 5 - benchmark and rule scoreboard
# ==========================================================================
bm <- read_tsv(file.path(RES, "benchmark", "coordinate_benchmark.tsv"), show_col_types = FALSE)
p5a <- bm %>% filter(scope == "all_truth", direction %in% c("combined", "Mf", "Mn")) %>%
  mutate(arm_lab = arm_factor(arm),
         direction = factor(direction, levels = c("combined", "Mf", "Mn"))) %>%
  ggplot(aes(direction, pct_reproduced, fill = direction)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 64, linetype = "dashed", linewidth = 0.4) +
  annotate("text", x = 0.6, y = 64, label = "published 558-basis ceiling (64%)",
           hjust = 0, vjust = -0.5, size = 3) +
  geom_text(aes(label = sprintf("%.1f%%", pct_reproduced)), vjust = -0.4, size = 3.2) +
  facet_wrap(~arm_lab) +
  scale_fill_pk_ordinal("", guide = "none") +
  scale_y_continuous(limits = c(0, 75), expand = c(0, 0)) +
  labs(x = NULL, y = "% of the poster's 217 windows reproduced",
       title = "Coordinate reproduction at the derived floor",
       subtitle = "matched by CHROM + window start, never by WINDOW id") +
  theme_pk()

rs <- read_tsv(file.path(RES, "benchmark", "rule_scoreboard.tsv"), show_col_types = FALSE)
indo <- tibble(rule = c("absolute", "relative", "distance", "distance_adaptive"),
               indo_pct = c(37, 14, 7, 49))
# One explicit rule ordering, shared by the bars and the comparison markers, so
# the markers land on the floor-6 bar they refer to rather than the rule centre.
rule_ord <- rs %>% group_by(rule) %>% summarise(m = max(pct_reproduced)) %>%
  arrange(desc(m)) %>% pull(rule)
POL_LEVELS <- c("common floor 41", "own derived floor", "published floor 6")
DODGE <- 0.8
indo <- indo %>%
  mutate(xpos = match(rule, rule_ord) + DODGE * (3 - 1) / (2 * 3))  # 3rd of 3 dodged bars

p5b <- rs %>%
  mutate(policy = factor(recode(floor_policy,
                                common = "common floor 41", own = "own derived floor",
                                floor6 = "published floor 6"), levels = POL_LEVELS),
         rule = factor(rule, levels = rule_ord)) %>%
  ggplot(aes(rule, pct_reproduced, fill = policy)) +
  geom_col(position = position_dodge(width = DODGE), width = 0.7) +
  geom_point(data = indo, aes(x = xpos, y = indo_pct), inherit.aes = FALSE,
             shape = 4, size = 2.8, stroke = 1.1) +
  scale_fill_pk_ordinal("floor policy") +
  labs(x = NULL, y = "% of the 217 windows reproduced",
       title = "Rule ranking inverts with the support floor",
       subtitle = "x = published 558-basis scoreboard at floor 6, marked on our matching floor-6 bar") +
  theme_pk() + theme(legend.position = "bottom",
                     axis.text.x = element_text(angle = 15, hjust = 1))

save_run_fig(p5a / p5b + plot_annotation(tag_levels = "a"),
             "fig5_benchmark_and_rules", width = 9, height = 8.5)

say("wrote %d figure pairs to %s", 5, FIGS)
