# _setup.R — shared config, palette, theme, figure-saver for the Pk pop-gen project.
# Source at the top of each numbered script:  source(here::here("scripts/_setup.R"))
# Env-independent (no heavy pkgs); safe to source before the analysis libraries.

suppressMessages({
  library(here)
  library(ggplot2)
})

# ---- paths (single source of truth; here() resolves the project root via .here)
PROJ      <- here::here()
DATA_RAW  <- here::here("data", "raw")        # imported, read-only
DATA_PROC <- here::here("data", "processed")  # script-generated derivations
FIG_DIR   <- here::here("figures")            # vector + png figure outputs
RES_DIR   <- here::here("results")            # tables / stats
for (d in c(DATA_PROC, FIG_DIR, RES_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ---- accessible palette (Okabe-Ito; colour-blind safe)
okabe_ito <- c(
  black  = "#000000", orange = "#E69F00", skyblue = "#56B4E9", green = "#009E73",
  yellow = "#F0E442", blue   = "#0072B2", vermillion = "#D55E00", purple = "#CC79A7"
)

# Genomic-cluster colours. INFORMED CAVEAT: assigned from Okabe-Ito for
# accessibility — NOT verified against the poster's exact hues. Adjust here if a
# specific mapping is required; every figure inherits from this one definition.
cluster_cols <- c(
  Peninsular = unname(okabe_ito["green"]),
  Mf         = unname(okabe_ito["blue"]),
  Mn         = unname(okabe_ito["vermillion"])
)

# continuous fills (rasters/densities): inferno via viridis option
seq_option <- "inferno"

# ---- shared theme (clean, publication-leaning)
theme_pk <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      axis.title       = element_text(face = "bold"),
      legend.position  = "right",
      plot.title       = element_text(face = "bold"),
      strip.text       = element_text(face = "bold")
    )
}

# ---- figure saver: write BOTH vector (.svg) and raster (.png) into figures/
save_fig <- function(plot, name, width = 8, height = 6, dpi = 300) {
  png <- file.path(FIG_DIR, paste0(name, ".png"))
  svg <- file.path(FIG_DIR, paste0(name, ".svg"))
  ggsave(png, plot, width = width, height = height, dpi = dpi)
  ggsave(svg, plot, width = width, height = height)        # editable handoff
  invisible(c(png = png, svg = svg))
}

message("Pk pop-gen setup loaded · root: ", PROJ)
