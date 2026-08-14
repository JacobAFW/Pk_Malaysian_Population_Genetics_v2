# Stage B-8: rehh iHS per cluster (Mf/Mn/Peninsular).
# Ports scripts/gadi/rehh.R to run for all 3 clusters in one R session, writing
# candidate regions + iHS/pvalue Manhattan plots via save_fig().
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(rehh); library(R.utils); library(data.table)
})
source(here::here("scripts/_setup.R"))

CLUSTERS  <- list(Mf = "Mf", Mn = "Mn", Peninsular = "Pen")   # display name -> dir name
CONTIGS   <- sprintf("ordered_PKNH_%02d_v2", 1:14)
SEL_DIR   <- file.path(DATA_PROC, "selection")
OUT_TSV   <- file.path(SEL_DIR, "candidate_regions_iHS_clusters.tsv")
dir.create(file.path(FIG_DIR, "stage_b"), showWarnings = FALSE, recursive = TRUE)

selection_plot <- function(ihs_df, model, y_var, ylab_str, name) {
  pd <- ihs_df |>
    tidyr::drop_na(all_of(y_var)) |>
    mutate(CHR = str_remove(CHR, "ordered_PKNH_"),
           CHR = as.integer(str_remove(CHR, "_v2"))) |>
    arrange(CHR, POSITION) |>
    mutate(ROW = row_number(), CHR = factor(CHR))
  x_axis <- pd |> group_by(CHR) |> summarise(ROW = median(ROW))
  cols <- rep(c("#39568CFF", "#29AF7FFF"), 7)
  p <- ggplot(pd, aes(x = ROW, y = .data[[y_var]], colour = CHR)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_colour_manual(values = cols) +
    scale_x_continuous(breaks = x_axis$ROW, labels = x_axis$CHR) +
    xlab("Chromosome") + ylab(ylab_str) +
    theme_pk() + theme(legend.position = "none") +
    ggtitle(paste0(model, " — ", ylab_str))
  save_fig(p, file.path("stage_b", name), width = 12, height = 4)
}

all_regions <- list()

for (nm in names(CLUSTERS)) {
  dirtag <- CLUSTERS[[nm]]
  CDIR <- file.path(SEL_DIR, dirtag)
  message(sprintf("[16b] === %s (%s) ===", nm, CDIR))

  wgscan_parts <- list()
  for (i in CONTIGS) {
    hap_file <- file.path(CDIR, paste0(i, ".vcf.gz"))
    if (!file.exists(hap_file) || file.info(hap_file)$size < 1000) {
      message(sprintf("[16b]   skip %s (missing/empty)", i)); next
    }
    hh <- tryCatch(
      data2haplohh(hap_file = hap_file, chr.name = i,
                   polarize_vcf = FALSE, min_perc_geno.mrk = 100,
                   min_maf = 0.05, vcf_reader = "data.table",
                   verbose = FALSE),
      error = function(e) { message(sprintf("[16b]   err %s: %s", i, conditionMessage(e))); NULL }
    )
    if (is.null(hh) || nhap(hh) == 0 || nmrk(hh) == 0) {
      message(sprintf("[16b]   empty hh for %s (nhap=%d nmrk=%d)", i,
                      if(is.null(hh)) NA else nhap(hh), if(is.null(hh)) NA else nmrk(hh)))
      next
    }
    scan <- scan_hh(hh, discard_integration_at_border = FALSE)
    wgscan_parts[[i]] <- scan
    message(sprintf("[16b]   %s: %d markers, %d haplotypes", i, nmrk(hh), nhap(hh)))
  }
  wgscan <- bind_rows(wgscan_parts)
  if (nrow(wgscan) == 0) { message(sprintf("[16b] %s: no data, skip", nm)); next }
  wgscan <- as.data.frame(wgscan) |> tidyr::drop_na()

  wgs_ihs <- ihh2ihs(wgscan, min_maf = 0.05, freqbin = 0.025)
  ihs_out <- wgs_ihs$ihs

  # Save the raw per-marker iHS for record
  write_tsv(as_tibble(ihs_out), file.path(SEL_DIR, sprintf("iHS_%s.tsv", dirtag)))

  # Candidate regions (windowed extreme-iHS clusters, pvalue-based threshold 4)
  wgs_ihs_2 <- wgs_ihs
  wgs_ihs_2$ihs <- as.data.frame(wgs_ihs_2$ihs) |>
    mutate(CHR = str_remove(CHR, "ordered_PKNH_"),
           CHR = as.numeric(str_remove(CHR, "_v2")))
  regs <- calc_candidate_regions(wgs_ihs_2, threshold = 4, pval = TRUE,
                                 window_size = 10000, overlap = 1000,
                                 min_n_extr_mrk = 3)
  if (nrow(regs) > 0) regs$Stat <- paste0("iHS_", dirtag)
  all_regions[[nm]] <- regs
  message(sprintf("[16b] %s: %d candidate regions", nm, nrow(regs)))

  # Manhattan plots
  selection_plot(as.data.frame(ihs_out), nm, "IHS",       "iHS",                       sprintf("16_ihs_%s", dirtag))
  selection_plot(as.data.frame(ihs_out), nm, "LOGPVALUE", "iHS -log10(p)",             sprintf("16_ihs_%s_pvalue", dirtag))
}

combined <- bind_rows(all_regions)
write_tsv(combined, OUT_TSV)
message(sprintf("[16b] wrote %s (total %d regions across clusters)",
                OUT_TSV, nrow(combined)))
cat("[16b] cluster breakdown:\n")
print(combined |> count(Stat))
