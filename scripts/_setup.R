# _setup.R — shared config, palette, theme, figure-saver for the Pk pop-gen project.
# Source at the top of each numbered script:  source(here::here("scripts/_setup.R"))
# Env-independent (no heavy pkgs); safe to source before the analysis libraries.
#
# STYLE SOURCE OF TRUTH. Every figure in 01/02 inherits palette + theme from here, so a
# restyle is a one-file change. Spec: viridis family throughout, theme_bw(11) with the grid
# and border stripped for non-map plots, a dark-panel theme for maps, png+pdf+svg output.

suppressMessages({
  library(here)
  library(ggplot2)
  library(viridis)
})

# ---- paths (single source of truth; here() resolves the project root via .here)
PROJ      <- here::here()
DATA_RAW  <- here::here("data", "raw")        # imported, read-only
DATA_PROC <- here::here("data", "processed")  # script-generated derivations
FIG_DIR   <- here::here("figures")            # vector + png figure outputs
RES_DIR   <- here::here("results")            # tables / stats
for (d in c(DATA_PROC, FIG_DIR, RES_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ---- palette -----------------------------------------------------------------
# Genomic clusters: viridis discrete. These three hues are the ones already used in the
# Sabah admixture map, kept as ONE named mapping so every figure agrees.
cluster_cols <- c(
  Peninsular = "#440154FF",   # viridis 0.00
  Mf         = "#73D055FF",   # viridis 0.85
  Mn         = "#39568CFF"    # viridis 0.30
)

# Ordinal fills (ADMIXTURE K, introgression / degree tertiles): inferno, trimmed so neither
# the near-black nor the near-white end of the ramp is used.
ORD_OPTION <- "inferno"; ORD_BEGIN <- 0.2; ORD_END <- 0.9

# Continuous rasters (ALC vector density, tree cover, elevation): viridis, matching the poster.
CONT_OPTION <- "viridis"

scale_colour_pk_cluster <- function(name = "Genomic cluster", ...)
  scale_colour_manual(values = cluster_cols, name = name, ...)
scale_color_pk_cluster <- scale_colour_pk_cluster
scale_fill_pk_cluster <- function(name = "Genomic cluster", ...)
  scale_fill_manual(values = cluster_cols, name = name, ...)

scale_fill_pk_ordinal <- function(...)
  scale_fill_viridis_d(option = ORD_OPTION, begin = ORD_BEGIN, end = ORD_END, ...)
scale_colour_pk_ordinal <- function(...)
  scale_colour_viridis_d(option = ORD_OPTION, begin = ORD_BEGIN, end = ORD_END, ...)

scale_fill_pk_continuous <- function(...) scale_fill_viridis_c(option = CONT_OPTION, ...)

# Discrete viridis for open-ended categories (e.g. clonal cluster IDs, K components).
pk_viridis_d <- function(n) viridis::viridis(n)

# ---- themes ------------------------------------------------------------------
# Non-map plots: bars, scatter, Manhattan, regression.
theme_pk <- function(base_size = 11, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid      = element_blank(),
      panel.border    = element_blank(),
      axis.line       = element_line(colour = "black", linewidth = 0.4),
      axis.ticks      = element_line(colour = "black", linewidth = 0.4),
      text            = element_text(colour = "black", family = base_family),
      axis.text       = element_text(colour = "black"),
      legend.position = "right",
      legend.key      = element_blank(),
      strip.background = element_blank(),
      strip.text      = element_text(colour = "black", face = "bold"),
      plot.title      = element_text(hjust = 0.5)
    )
}

# Maps: dark panel so viridis fills and light sample points read against it.
# NOTE: `theme_linedraw() + theme_dark()` in the original chunks is effectively just
# theme_dark(), because adding a COMPLETE theme replaces the previous one wholesale. The
# linedraw part is reinstated explicitly below (black border + ticks), which is the evident
# intent. Change here to restyle every map at once.
theme_pk_map <- function(base_size = 11, base_family = "sans") {
  theme_dark(base_size = base_size, base_family = base_family) +
    theme(
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid      = element_blank(),
      axis.ticks      = element_line(colour = "black", linewidth = 0.4),
      text            = element_text(colour = "black", family = base_family),
      axis.text       = element_text(colour = "black"),
      legend.key      = element_blank(),
      strip.background = element_blank(),
      strip.text      = element_text(colour = "black", face = "bold"),
      plot.title      = element_text(hjust = 0.5),           # centred title
      plot.margin     = unit(c(0.4, 0.4, 0.4, 0.4), "cm")
    )
}

# ---- basemap -----------------------------------------------------------------
# malariaAtlas::getShp() no longer works: the geoserver WFS endpoint it calls
# (https://malariaatlas.org/geoserver/Explorer/ows) returns HTTP 404 — the service was
# retired, so this fails for everyone, not just offline. Cache the polygon on first
# successful build and fall back to rnaturalearth so the maps stay reproducible.
pk_malaysia_shp <- function(cache = file.path(DATA_PROC, "malaysia_admin_sf.rds"),
                            quiet = FALSE) {
  if (file.exists(cache)) return(readRDS(cache))
  shp <- NULL
  if (requireNamespace("malariaAtlas", quietly = TRUE)) {
    # NB getShp() wraps its own call in try() and RETURNS the try-error object instead of
    # signalling it, so tryCatch(error=) never fires -- check the returned class as well.
    shp <- tryCatch(
      malariaAtlas::getShp(country = "Malaysia",
                           admin_level = c("admin0", "admin1", "admin2")),
      error = function(e) NULL)
    if (inherits(shp, "try-error") ||
        !(inherits(shp, "sf") || inherits(shp, "Spatial"))) shp <- NULL
  }
  src <- "malariaAtlas::getShp"
  if (is.null(shp)) {
    src <- "rnaturalearth::ne_states  [SUBSTITUTE: malariaAtlas geoserver returns HTTP 404]"
    shp <- rnaturalearth::ne_states(country = "malaysia", returnclass = "sf")
  }
  if (!quiet) message("pk_malaysia_shp(): basemap from ", src)
  attr(shp, "pk_source") <- src
  saveRDS(shp, cache)
  shp
}

# The Anopheles leucosphyrus complex layer. malariaAtlas 1.0.0 dropped the `dataset_id`
# argument getRaster() was originally called with (it now takes `surface`), and the
# Explorer endpoint is gone, so the cached data frame is the only working source.
pk_anoph_raster_df <- function(cache = file.path(DATA_PROC, "anoph_raster_df.rds")) {
  if (file.exists(cache)) return(readRDS(cache))
  stop("anoph_raster_df cache missing at ", cache,
       " and malariaAtlas::getRaster() can no longer fetch it (API retired). ",
       "Restore the cache or re-fetch from the current MAP API.")
}

# ---- Q-matrix palettes -------------------------------------------------------
# A K-column Q matrix has UNLABELLED ancestry components; their column order is whatever the
# fit produced. Handing CreatePalette() the cluster colours in name order therefore assigns
# hues positionally and can silently mislabel the map. These helpers work out, from the
# truth cluster labels, which component each column actually is, and only use the named
# cluster colours when the correspondence is clean.
pk_q_cluster_order <- function(Q, sample_ids,
                               labels_file = file.path(DATA_PROC, "cluster_maf",
                                                       "cluster_labels.tsv"),
                               min_mean_q = 0.5) {
  lab <- utils::read.delim(labels_file, stringsAsFactors = FALSE)
  names(lab)[1:2] <- c("Sample", "Cluster")
  cl <- lab$Cluster[match(sample_ids, lab$Sample)]
  best <- character(ncol(Q)); score <- numeric(ncol(Q))
  for (j in seq_len(ncol(Q))) {
    mu <- tapply(Q[, j], cl, mean, na.rm = TRUE)
    best[j] <- names(mu)[which.max(mu)]; score[j] <- max(mu, na.rm = TRUE)
  }
  list(order = best, score = score,
       ok = !anyDuplicated(best) && all(score >= min_mean_q) && all(best %in% names(cluster_cols)))
}

pk_q_palette <- function(Q, sample_ids, n = 10, quiet = FALSE) {
  o <- try(pk_q_cluster_order(Q, sample_ids), silent = TRUE)
  if (!inherits(o, "try-error") && isTRUE(o$ok)) {
    if (!quiet) message("pk_q_palette(): Q columns anchored to ",
                        paste(o$order, collapse = " / "))
    cols <- unname(cluster_cols[o$order])
  } else {
    sc <- if (inherits(o, "try-error")) "unavailable" else paste(round(o$score, 2), collapse = ", ")
    warning("pk_q_palette(): Q components do NOT correspond 1:1 to the named clusters ",
            "(max mean Q per column: ", sc, "). Using neutral viridis component colours; ",
            "read this map as ancestry components, NOT as Mf/Mn/Peninsular.", call. = FALSE)
    cols <- viridis::viridis(ncol(Q))
  }
  tess3r::CreatePalette(cols, n)
}

# ---- figure saver: png (300 dpi) + pdf + svg ---------------------------------
save_fig <- function(plot, name, width = 8, height = 6, dpi = 300) {
  out <- c(
    png = file.path(FIG_DIR, paste0(name, ".png")),
    pdf = file.path(FIG_DIR, paste0(name, ".pdf")),
    svg = file.path(FIG_DIR, paste0(name, ".svg"))
  )
  ggsave(out[["png"]], plot, width = width, height = height, dpi = dpi)
  # base pdf(), not cairo_pdf(): cairo needs X11, which is not installed here, and
  # cairo_pdf() then warns "failed to load cairo DLL" and silently writes nothing.
  ggsave(out[["pdf"]], plot, width = width, height = height, device = grDevices::pdf)
  ggsave(out[["svg"]], plot, width = width, height = height)   # editable handoff
  invisible(out)
}

# Default theme for anything that does not set one explicitly.
theme_set(theme_pk())

message("Pk pop-gen setup loaded · root: ", PROJ, " · style: viridis / theme_pk")
