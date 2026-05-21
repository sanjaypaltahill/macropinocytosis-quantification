# Macropinocytosis Image Analysis Pipeline

A two-part toolkit for quantifying macropinocytosis from fluorescence microscopy images. A Python script (`analysis.py`) segments cells and measures dye uptake from multi-channel TIFF images; an R script (`graph-maker.R`) aggregates results across conditions and produces publication-ready figures.

---

## Table of Contents

1. [Overview](#overview)
2. [How It Works](#how-it-works)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [File Naming Conventions](#file-naming-conventions)
6. [Usage](#usage)
   - [Python Analysis Script](#python-analysis-script)
   - [R Plotting Script](#r-plotting-script)
7. [Output Files](#output-files)
8. [Tunable Parameters](#tunable-parameters)
9. [Metrics Explained](#metrics-explained)
10. [Tips and Troubleshooting](#tips-and-troubleshooting)

---

## Overview

Macropinocytosis is a form of bulk endocytosis in which cells engulf extracellular fluid and its contents into large vesicles. This pipeline quantifies macropinocytic uptake at the single-cell level from 63× widefield or confocal fluorescence microscopy images. It is designed to handle multi-condition experiments (e.g., dose–response, knockdown vs. control) and produces both per-cell measurements and condition-level summary statistics.

**Typical workflow:**

```
Raw TIFF images (3 channels per field of view)
        │
        ▼
analysis.py   ←── segments cells, measures dye intensity
        │
        ▼
  results/*_per_cell.csv      ←── one row per cell, per image
        │
        ▼
graph-maker.R       ←── bar/dot plot with SEM across conditions
        │
        ▼
  macropinocytosis_plot.png/.pdf
```

---

## How It Works

### Cell Segmentation (`analysis.py`)

The pipeline uses a three-channel imaging scheme:

| Channel | Default suffix | Role |
|---------|---------------|------|
| Actin (red) | `_ch00` | Cell boundary marker |
| Dye / FITC / TMR (green) | `_ch01` | Macropinocytosis readout |
| DAPI (blue) | `_ch02` | Nuclear stain |

Segmentation proceeds in five steps:

1. **Nucleus detection** — The DAPI channel is Gaussian-smoothed and thresholded with Otsu's method. Morphological closing fills holes in nuclear blobs. Each connected component is treated as one nucleus (optional watershed splitting is available for densely packed cells).

2. **Nucleus filtering** — Very small objects (debris, fragments) are discarded based on a minimum area threshold. Nuclei touching the image border can optionally be excluded.

3. **Cell territory definition** — The actin channel is smoothed and thresholded at a user-defined percentile. Pixels above this threshold are treated as "cell territory" — the space in which nuclei are allowed to expand.

4. **Cell body expansion** — Each nucleus is expanded outward up to a maximum radius, but expansion is constrained to the actin-positive territory. This prevents cell regions from bleeding into background or neighboring cells.

5. **Intensity measurement** — Raw dye pixel values are extracted for each cell mask. No background subtraction is applied, so measurements represent total signal within the cell boundary.

### Aggregation and Plotting (`graph-maker.R`)

The R script reads all `*_per_cell.csv` files from the `results/` subfolder of each specified experiment directory, stacks them, and computes per-condition mean ± SEM. It produces a dot plot (mean ± SEM) using `ggplot2` in a clean, Nature-style theme.

---

## Requirements

### Python (analysis)

- Python ≥ 3.9
- `tifffile`
- `imagecodecs`
- `numpy`
- `pandas`
- `matplotlib`
- `scikit-image`
- `scipy`
- `openpyxl`

### R (plotting)

- R ≥ 4.0
- `ggplot2`
- `dplyr`
- `readr`
- `ggbeeswarm`
- `scales`

Missing R packages are installed automatically when the script is run.

---

## Installation

```bash
# Clone or download the repository, then install Python dependencies:
pip install tifffile imagecodecs numpy pandas matplotlib scikit-image scipy openpyxl
```

No additional installation is needed for the R script — it self-installs missing packages on first run.

---

## File Naming Conventions

Each field of view must be represented by **three separate single-channel TIFF files** sharing a common base name and differing only in their channel suffix:

```
<basename>_ch00.tif    →  actin  (cell boundary)
<basename>_ch01.tif    →  dye    (uptake signal)
<basename>_ch02.tif    →  DAPI   (nuclei)
```

**Example:**

```
experiment01_pos001_ch00.tif
experiment01_pos001_ch01.tif
experiment01_pos001_ch02.tif
experiment01_pos002_ch00.tif
...
```

The script anchors on the actin suffix (`_ch00` by default) and locates the corresponding dye and DAPI siblings automatically. If your files use different suffixes, override them at the command line (see [Usage](#usage) below).

> **Alternative mode:** If your files do not share a common base name, set `TRIPLET_MODE = "sorted_order"` inside the script. Files will then be grouped as consecutive sorted triples (files 0–2 → image 1, files 3–5 → image 2, etc.).

---

## Usage

### Python Analysis Script

**Basic usage** (default channel mapping: `_ch00` = actin, `_ch01` = dye, `_ch02` = DAPI):

```bash
python analysis.py /path/to/image/folder
```

**Custom channel mapping:**

```bash
python analysis.py /path/to/image/folder \
    --actin _ch02 \
    --dapi  _ch00 \
    --dye   _ch01
```

**Arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `folder` | *(required)* | Path to the folder containing channel TIFF files |
| `--actin` | `_ch00` | Filename suffix for the actin channel |
| `--dapi` | `_ch02` | Filename suffix for the DAPI channel |
| `--dye` | `_ch01` | Filename suffix for the dye/FITC channel |

Run the script once per experiment folder (i.e., once per condition). Results are written to a `results/` subdirectory inside the input folder.

---

### R Plotting Script

Edit the **USER CONFIGURATION** block near the top of `graph-maker.R`:

```r
# Paths to experiment folders — one per condition
FOLDERS <- c(
  "/path/to/condition_1",
  "/path/to/condition_2",
  "/path/to/condition_3"
)

# Display labels (left to right on the plot)
LABELS <- c("Control", "Treatment A", "Treatment B")

# Bar/point colours (one hex code per condition)
COLOURS <- c("#4878CF", "#E84646", "#56A956")

# Metric: "integrated_intensity" (recommended) or "mean_intensity"
METRIC <- "integrated_intensity"

# Axis labels and plot title
X_LABEL    <- "Condition"
Y_LABEL    <- "Integrated Intensity (AU)"
PLOT_TITLE <- "Macropinocytosis Uptake"
```

Then run:

```bash
Rscript graph-maker.R
```

Or source the script inside RStudio.

---

## Output Files

### Per-image outputs (in `<folder>/results/`)

| File | Description |
|------|-------------|
| `<basename>_per_cell.csv` | One row per detected cell; columns: `image`, `cell_id`, `mean_intensity`, `integrated_intensity`, `integrated_intensity_per_area`, `pixel_count` |
| `side_by_side/<basename>_analysis.png` | Three-panel QC figure: original composite with outlines, cell mask overlay, dye heat-map |

### Folder-level summary (in `<folder>/results/`)

| File | Description |
|------|-------------|
| `summary_table.csv` | Per-image summary: `n_cells`, `mean_integrated_intensity`, `total_integrated_intensity`, plus a final `TOTALS` row |
| `summary_table.xlsx` | Same as above in Excel format |

### Plotting outputs (written to the parent of the first `FOLDERS` entry)

| File | Description |
|------|-------------|
| `macropinocytosis_plot.png` | Publication-ready bar/dot plot (mean ± SEM) |
| `macropinocytosis_plot.pdf` | Vector PDF version of the same plot |
| `macropinocytosis_plot_summary_stats.csv` | Per-condition mean, SD, SEM, and n |

---

## Tunable Parameters

The following constants can be edited directly at the top of `analysis.py`:

| Parameter | Default | Effect |
|-----------|---------|--------|
| `DAPI_BLUR_SIGMA` | `4` | Gaussian smoothing of DAPI before thresholding; increase to merge fragmented nuclei |
| `DAPI_CLOSE_RADIUS` | `8` | Morphological closing radius; increase to fill larger holes in nuclear masks |
| `MIN_NUCLEUS_AREA` | `5000` | Minimum nucleus size in px²; increase to remove larger debris |
| `EXCLUDE_EDGE_NUCLEI` | `False` | If `True`, discard any nucleus touching the image border |
| `WATERSHED_MIN_DIST` | `20` | Minimum peak separation for watershed splitting; set to `0` to disable splitting |
| `CELL_EXPANSION_PX` | `150` | Maximum outward expansion from nucleus edge (pixels) |
| `ACTIN_PERCENTILE` | `85` | Actin threshold percentile; lower values allow more permissive cell body expansion |

---

## Metrics Explained

Three complementary intensity metrics are reported for each cell:

**`mean_intensity`**
The average pixel value within the cell mask. Equivalent to uptake density; comparable across cells of different sizes, but does not capture the total amount of material internalized.

**`integrated_intensity`**
The sum of all pixel values within the cell mask. Represents total dye uptake per cell. This is the **recommended primary metric** when cell size varies between conditions (e.g., after knockdown or treatment), as it reflects total endocytic load rather than concentration.

**`integrated_intensity_per_area`**
Integrated intensity divided by cell area in pixels. Mathematically equivalent to `mean_intensity`; included for convenience when both uptake density and total uptake are reported side by side.

> **Which metric to use?** Use `integrated_intensity` as the default. Switch to `mean_intensity` / `integrated_intensity_per_area` only when you specifically want to control for cell size and report uptake concentration.

---

## Tips and Troubleshooting

**Too many / too few nuclei detected**
Adjust `DAPI_BLUR_SIGMA` and `MIN_NUCLEUS_AREA`. Inspect the QC panels in `side_by_side/` to evaluate segmentation quality visually.

**Cells being merged together**
Enable or reduce `WATERSHED_MIN_DIST` to split closely touching nuclei. Alternatively, lower `CELL_EXPANSION_PX` to reduce over-expansion of cell bodies.

**Cell bodies expanding into background**
Increase `ACTIN_PERCENTILE` (e.g., from 85 to 90) to raise the actin threshold and restrict expansion to higher-confidence cell territory.

**Files not being found**
Confirm that file suffixes match the `--actin`, `--dye`, and `--dapi` arguments. The script prints a warning for any expected sibling file it cannot locate.

**Z-stack TIFFs**
The loader automatically max-projects along the Z axis. No pre-processing is required.

**sgRNA / condition notes**
Per the experimental design of this project: use `sgRNA2` for SNAP experiments and `sgRNA1` for PELP1/AMBRA experiments.

---

## License

MIT License. See `LICENSE` for details.
