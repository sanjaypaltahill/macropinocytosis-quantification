#!/usr/bin/env Rscript
# =============================================================================
#  plot_macropinocytosis.R
#  Macropinocytosis Uptake Plot — Multi-Condition Summary
# =============================================================================
#
#  PURPOSE
#  -------
#  Reads per-cell CSVs produced by analyze_macropinocytosis.py from an
#  arbitrary number of experiment folders, then plots a bar chart with
#  individual data points and SEM error bars.
#
#  USAGE
#  -----
#  1.  Edit the "USER CONFIGURATION" section below (folders, labels, colours).
#  2.  Run:
#        Rscript plot_macropinocytosis.R
#      or source it inside RStudio.
#
#  WHAT IT READS
#  -------------
#  For each folder, it globs all  <folder>/results/*_per_cell.csv  files
#  and stacks them.  The "integrated_intensity" column is used as the
#  primary metric (total dye uptake per cell, accounts for cell-size
#  variation).  Change METRIC below to use "mean_intensity" instead.
#
#  OUTPUTS
#  -------
#  A publication-ready PNG (and optionally PDF) saved to OUTPUT_DIR.
#  A tidy CSV of the per-condition summary stats (mean ± SEM, n).
#
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# 0.  PACKAGES
# ─────────────────────────────────────────────────────────────────────────────
required_packages <- c("ggplot2", "dplyr", "readr", "ggbeeswarm", "scales")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}


# ═════════════════════════════════════════════════════════════════════════════
#  USER CONFIGURATION  ← edit everything in this block
# ═════════════════════════════════════════════════════════════════════════════

# ── Folders containing a "results/" sub-directory with *_per_cell.csv files ─
#    Add or remove entries freely; order = left-to-right on the plot.
FOLDERS <- c(
  "/Users/spaltahill/Library/CloudStorage/GoogleDrive-sanjay01@stanford.edu/My Drive/Macropinocytosis Project/SF188 PELP1 pInducer Assay/0 ug doxy",
  "/Users/spaltahill/Library/CloudStorage/GoogleDrive-sanjay01@stanford.edu/My Drive/Macropinocytosis Project/SF188 PELP1 pInducer Assay/1 ug doxy",
  "/Users/spaltahill/Library/CloudStorage/GoogleDrive-sanjay01@stanford.edu/My Drive/Macropinocytosis Project/SF188 PELP1 pInducer Assay/2 ug doxy"
)

# ── Display labels — one per folder, same order ──────────────────────────────
LABELS <- c(
  "0 ug",
  "1 ug",
  "2 ug"
)

# ── Bar/point colours — one hex (or R colour name) per folder ────────────────
COLOURS <- c(
  "#4878CF",
  "#4878CF",
  "#4878CF"
)

# ── Metric to plot ────────────────────────────────────────────────────────────
#    "integrated_intensity"  ← recommended (total uptake per cell)
#    "mean_intensity"        ← average pixel intensity per cell
METRIC <- "mean_intensity"

# ── Axis Labels ─────────────────────────────────────────────────────────────
X_LABEL <- "[Doxycycline]"
Y_LABEL <- "Mean Intensity (AU)"

# ── Plot title (set to "" for no title) ──────────────────────────────────────
PLOT_TITLE <- "SF188 PELP1 pInducer Assay"

# ── Output ────────────────────────────────────────────────────────────────────
OUTPUT_DIR <- dirname(FOLDERS[1])          # folder where output files are written
OUTPUT_STEM <- "macropinocytosis_plot"   # filename without extension
SAVE_PNG    <- TRUE
SAVE_PDF    <- TRUE
PNG_WIDTH   <- 5
PNG_HEIGHT  <- 4
BASE_FONT_SIZE <- 14

# ── Plot aesthetics (optional fine-tuning) ────────────────────────────────────
BAR_ALPHA      <- 0.75   # bar transparency  (0 = transparent, 1 = solid)
POINT_SIZE     <- 1.5    # jittered individual cell dot size
POINT_ALPHA    <- 0.45   # dot transparency
ERRORBAR_WIDTH <- 0.18   # width of SEM error-bar caps
BASE_FONT_SIZE <- 12

# ═════════════════════════════════════════════════════════════════════════════


# ─────────────────────────────────────────────────────────────────────────────
# 1.  VALIDATION
# ─────────────────────────────────────────────────────────────────────────────
stopifnot(
  "FOLDERS, LABELS, and COLOURS must all have the same length" =
    length(FOLDERS) == length(LABELS) && length(FOLDERS) == length(COLOURS)
)

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  message("Created output directory: ", OUTPUT_DIR)
}


# ─────────────────────────────────────────────────────────────────────────────
# 2.  LOAD DATA
# ─────────────────────────────────────────────────────────────────────────────
load_condition <- function(folder, label) {
  results_dir <- file.path(folder, "results")
  if (!dir.exists(results_dir)) {
    warning("No 'results' sub-folder found in: ", folder, " — skipping.")
    return(NULL)
  }
  csv_files <- list.files(
    path       = results_dir,
    pattern    = "_per_cell\\.csv$",
    full.names = TRUE
  )
  if (length(csv_files) == 0) {
    warning("No *_per_cell.csv files found in: ", results_dir, " — skipping.")
    return(NULL)
  }
  df <- lapply(csv_files, readr::read_csv, show_col_types = FALSE) |>
    dplyr::bind_rows()
  df$condition <- label
  df
}

raw_list <- mapply(
  load_condition,
  folder = FOLDERS,
  label  = LABELS,
  SIMPLIFY = FALSE
)

all_data <- dplyr::bind_rows(raw_list)

if (nrow(all_data) == 0) {
  stop("No data loaded. Check that FOLDERS point to valid experiment directories.")
}

# Enforce factor level order = order given by the user
all_data$condition <- factor(all_data$condition, levels = LABELS)

# Map condition → colour for ggplot
colour_map <- setNames(COLOURS, LABELS)

message(
  sprintf(
    "Loaded %d cells across %d condition(s).",
    nrow(all_data), length(unique(all_data$condition))
  )
)


# ─────────────────────────────────────────────────────────────────────────────
# 3.  SUMMARY STATISTICS  (mean ± SEM per condition)
# ─────────────────────────────────────────────────────────────────────────────
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

summary_df <- all_data |>
  dplyr::group_by(condition) |>
  dplyr::summarise(
    n       = dplyr::n(),
    mean    = mean(.data[[METRIC]], na.rm = TRUE),
    sd      = sd(.data[[METRIC]],   na.rm = TRUE),
    sem     = sem(.data[[METRIC]]),
    .groups = "drop"
  )

message("\n── Per-condition summary ──────────────────────────────────")
print(as.data.frame(summary_df))
message("───────────────────────────────────────────────────────────\n")

# Write summary CSV
summary_csv <- file.path(OUTPUT_DIR, paste0(OUTPUT_STEM, "_summary_stats.csv"))
readr::write_csv(summary_df, summary_csv)
message("Summary stats saved to: ", summary_csv)


# ─────────────────────────────────────────────────────────────────────────────
# 4.  PLOT
# ─────────────────────────────────────────────────────────────────────────────
p <- ggplot2::ggplot(
  summary_df,
  ggplot2::aes(x = condition, y = mean, colour = condition)
) +
  
  # Mean points
  ggplot2::geom_point(
    size = 1
  ) +
  
  # SEM error bars
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean - sem, ymax = mean + sem),
    width     = 0.10,
    linewidth = 1.0
  ) +
  
  # Colour scales
  ggplot2::scale_colour_manual(values = colour_map, guide = "none") +
  
  # Tighten spacing between conditions
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.08, 0.08))) +
  
  # Axis formatting
  ggplot2::scale_y_continuous(
    limits = c(0, NA),
    labels = scales::label_comma(),
    expand = ggplot2::expansion(mult = c(0, 0.05))
  ) +
  
  # Labels
  ggplot2::labs(
    title = PLOT_TITLE,
    x     =  X_LABEL, 
    y     = Y_LABEL
  ) +
  
  # Theme
  ggplot2::theme_classic(base_size = BASE_FONT_SIZE) +
  ggplot2::theme(
    
    # Title
    plot.title = ggplot2::element_text(
      face   = "bold",
      size   = 20,
      hjust  = 0.5,
      margin = ggplot2::margin(b = 14)
    ),
    
    # X labels
    axis.text.x = ggplot2::element_text(
      colour = "black",
      face   = "bold",
      size   = 16,
      margin = ggplot2::margin(t = 8)
    ),
    
    # Y tick labels
    axis.text.y = ggplot2::element_text(
      colour = "black",
      face   = "bold",
      size   = 14
    ),
    
    # Axis titles
    axis.title.x = ggplot2::element_text(
      face = "bold",
      size = 18
    ),
    
    axis.title.y = ggplot2::element_text(
      face   = "bold",
      size   = 18,
      margin = ggplot2::margin(r = 14)
    ),
    
    # Axes
    axis.line = ggplot2::element_line(
      colour   = "black",
      linewidth = 1.0
    ),
    
    axis.ticks = ggplot2::element_line(
      colour   = "black",
      linewidth = 0.9
    ),
    
    axis.ticks.length = grid::unit(0.22, "cm"),
    
    # Clean Nature-style background
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor   = ggplot2::element_blank(),
    
    # Reduce outer whitespace
    plot.margin = ggplot2::margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )
# ─────────────────────────────────────────────────────────────────────────────
# 5.  SAVE
# ─────────────────────────────────────────────────────────────────────────────
if (SAVE_PNG) {
  png_path <- file.path(OUTPUT_DIR, paste0(OUTPUT_STEM, ".png"))
  ggplot2::ggsave(
    filename = png_path,
    plot     = p,
    width    = PNG_WIDTH,
    height   = PNG_HEIGHT,
    dpi      = PNG_DPI,
    bg       = "white"
  )
  message("PNG saved to: ", png_path)
}

if (SAVE_PDF) {
  pdf_path <- file.path(OUTPUT_DIR, paste0(OUTPUT_STEM, ".pdf"))
  ggplot2::ggsave(
    filename = pdf_path,
    plot     = p,
    width    = PNG_WIDTH,
    height   = PNG_HEIGHT,
    device   = cairo_pdf
  )
  message("PDF saved to: ", pdf_path)
}

message("\nDone.")

