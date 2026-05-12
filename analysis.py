"""
Macropinocytosis 63x Microscopy Image Analysis Pipeline
=====================================================
Strategy
--------
  1.  Detect nuclei from the DAPI (blue) channel using Otsu threshold +
      connected-component labelling.  No watershed splitting — each connected
      DAPI blob is treated as exactly one cell nucleus.
  2.  Filter out debris (too small) and partial cells cut off by the image edge.
  3.  Expand each nucleus outward, constrained to actin (red) signal, to fill
      the full cell body up to the red membrane boundary.
  4.  Measure mean FITC/TMR (green) intensity inside each complete cell region.
  5.  Save a side-by-side PNG (original composite | cell masks | FITC heat-map)
      for every image, plus a summary CSV/XLSX across the whole folder.

Channel layout (standard RGB TIF from fluorescence microscope)
  index 0  →  Red   = actin  (cell boundary marker)
  index 1  →  Green = FITC / TMR  (signal of interest)
  index 2  →  Blue  = DAPI  (nuclear stain, used for segmentation)

USAGE
-----
    python analyze_macropinocytosis.py /path/to/folder

REQUIREMENTS
------------
    pip install tifffile imagecodecs numpy pandas matplotlib scikit-image scipy openpyxl

OUTPUTS  (written to <folder>/results/)
---------------------------------------
    side_by_side/<name>_analysis.png  – 3-panel figure per image
    <name>_per_cell.csv               – per-cell data for every image
    summary_table.csv / .xlsx         – n_cells + mean FITC intensity per image
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import ndimage
from skimage.color import label2rgb
from skimage.filters import gaussian, threshold_otsu
from skimage.measure import label as sklabel, regionprops
from skimage.morphology import closing, disk
from skimage.segmentation import expand_labels, find_boundaries


# ═══════════════════════════ TUNABLE PARAMETERS ═══════════════════════════════

# ── Channel indices (0-based) after loading as (C, H, W) ────────────────────
ACTIN_CH  = 0    # red   – defines cell body boundary
FITC_CH   = 1    # green – fluorescence signal to measure
DAPI_CH   = 2    # blue  – nuclei, used for cell counting & segmentation

# ── DAPI nucleus detection ───────────────────────────────────────────────────
DAPI_BLUR_SIGMA     = 4      # gaussian smoothing before threshold (pixels)
DAPI_CLOSE_RADIUS   = 8      # morphological closing to fill holes in nuclei
MIN_NUCLEUS_AREA    = 5000   # minimum nucleus area in px² (removes debris)
EXCLUDE_EDGE_NUCLEI = False   # if True, skip nuclei touching the image border

# ── Cell body expansion ───────────────────────────────────────────────────────
CELL_EXPANSION_PX   = 150    # max expansion radius from nucleus edge (pixels)
ACTIN_PERCENTILE    = 50     # actin pixels above this percentile define cell
                             # territory.  Lower = more permissive expansion.

# ── Output ────────────────────────────────────────────────────────────────────
RESULTS_DIR_NAME = "results"

# ══════════════════════════════════════════════════════════════════════════════


# ─── image loading ────────────────────────────────────────────────────────────

def load_tif(path: Path):
    """
    Load a TIF and return (actin, fitc, dapi) 2-D float32 arrays
    plus an (H, W, 3) RGB array for display.

    Handles:
      (H, W, C)    interleaved  – most colour camera / ImageJ RGB exports
      (C, H, W)    channel-first – OMERO / multi-channel stacks
      (Z, C, H, W) Z-stack      – max-projected along Z first
    """
    raw = tifffile.imread(str(path))
    img = raw.astype(np.float32)

    if img.ndim == 4:                  # Z-stack → max-project
        img = img.max(axis=0)

    if img.ndim != 3:
        raise ValueError(f"Unexpected image shape {raw.shape} in {path.name}")

    h, w = img.shape[0], img.shape[1]
    last = img.shape[2]
    if last <= 4 and last < h and last < w:   # (H, W, C) → (C, H, W)
        img = np.moveaxis(img, -1, 0)

    if img.shape[0] < 3:
        raise ValueError(
            f"{path.name}: only {img.shape[0]} channel(s) — need R/G/B.")

    actin = img[ACTIN_CH]
    fitc  = img[FITC_CH]
    dapi  = img[DAPI_CH]
    rgb   = np.stack([img[0], img[1], img[2]], axis=-1)
    return actin, fitc, dapi, rgb


# ─── helpers ──────────────────────────────────────────────────────────────────

def normalise(arr: np.ndarray, pmin: float = 1.0, pmax: float = 99.5) -> np.ndarray:
    lo, hi = np.percentile(arr, pmin), np.percentile(arr, pmax)
    if hi <= lo:
        return np.zeros_like(arr, dtype=np.float32)
    return np.clip((arr - lo) / (hi - lo), 0.0, 1.0).astype(np.float32)


def normalise_rgb(rgb: np.ndarray) -> np.ndarray:
    out = np.zeros_like(rgb, dtype=np.float32)
    for c in range(rgb.shape[2]):
        out[..., c] = normalise(rgb[..., c])
    return out


# ─── segmentation ─────────────────────────────────────────────────────────────

def segment_cells(actin: np.ndarray, dapi: np.ndarray):
    """
    Returns (nuclear_labels, cell_labels) — 2-D int32 arrays,
    0 = background, 1..N = individual cells.

    Each DAPI connected component = one cell nucleus (no splitting).
    Nuclei are then expanded into actin-positive territory to define
    the full cell body.
    """
    H, W = dapi.shape

    # ── Step 1: threshold DAPI ───────────────────────────────────────────────
    dapi_smooth = gaussian(dapi, sigma=DAPI_BLUR_SIGMA)
    thresh      = threshold_otsu(dapi_smooth)
    nuclear_bin = dapi_smooth > thresh
    nuclear_bin = closing(nuclear_bin, disk(DAPI_CLOSE_RADIUS))

    # ── Step 2: connected-component labelling (no watershed splitting) ───────
    raw_labels, _ = sklabel(nuclear_bin, return_num=True)

    # ── Step 3: filter nuclei ────────────────────────────────────────────────
    nuclear_labels = np.zeros_like(raw_labels, dtype=np.int32)
    new_id = 1
    for p in regionprops(raw_labels):
        r0, c0, r1, c1 = p.bbox
        touches_edge = (r0 == 0 or c0 == 0 or r1 == H or c1 == W)
        if p.area < MIN_NUCLEUS_AREA:
            continue
        if EXCLUDE_EDGE_NUCLEI and touches_edge:
            continue
        nuclear_labels[raw_labels == p.label] = new_id
        new_id += 1

    # ── Step 4: define cell territory from actin signal ──────────────────────
    actin_smooth  = gaussian(actin, sigma=2)
    actin_thresh  = np.percentile(actin_smooth, ACTIN_PERCENTILE)
    cell_territory = (actin_smooth > actin_thresh) | (nuclear_labels > 0)

    # ── Step 5: expand nuclei outward, clipped to cell territory ─────────────
    cell_labels = expand_labels(nuclear_labels, distance=CELL_EXPANSION_PX)
    cell_labels[~cell_territory] = 0

    return nuclear_labels, cell_labels.astype(np.int32)


# ─── measurement ──────────────────────────────────────────────────────────────

def measure_intensity(fitc: np.ndarray, cell_labels: np.ndarray) -> pd.DataFrame:
    records = []
    for cid in range(1, int(cell_labels.max()) + 1):
        px = fitc[cell_labels == cid]
        if px.size == 0:
            continue
        records.append({
            "cell_id":        cid,
            "mean_intensity": round(float(px.mean()), 4),
            "pixel_count":    int(px.size),
        })
    return pd.DataFrame(records)


# ─── visualisation ────────────────────────────────────────────────────────────

def make_side_by_side(
    rgb:            np.ndarray,
    fitc:           np.ndarray,
    nuclear_labels: np.ndarray,
    cell_labels:    np.ndarray,
    cell_df:        pd.DataFrame,
    image_name:     str,
    out_path:       Path,
) -> None:
    """
    Three-panel figure:
      1. Contrast-stretched RGB composite with nucleus (blue) and
         cell-boundary (yellow) outlines drawn on top.
      2. Cell regions, each uniquely coloured.
      3. FITC intensity heat-map restricted to cell regions.
    """
    composite = normalise_rgb(rgb)

    nuc_bounds  = find_boundaries(nuclear_labels, mode="outer")
    cell_bounds = find_boundaries(cell_labels,    mode="outer")
    comp_ann = composite.copy()
    comp_ann[cell_bounds] = [1.0, 0.85, 0.0]   # yellow = cell edge
    comp_ann[nuc_bounds]  = [0.2, 0.65, 1.0]   # blue   = nucleus edge

    mask_rgb    = label2rgb(cell_labels, bg_label=0, bg_color=(0.08, 0.08, 0.08))

    fitc_n      = normalise(fitc)
    fitc_masked = fitc_n.copy()
    fitc_masked[cell_labels == 0] = 0.0
    fitc_rgb    = plt.cm.hot(fitc_masked)[..., :3]

    n_cells = len(cell_df)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), facecolor="#0e0e0e")
    fig.suptitle(
        f"{image_name}   |   {n_cells} cells",
        color="white", fontsize=13, fontweight="bold", y=1.01,
    )
    panels = [
        (comp_ann,  "Original composite  (blue=nucleus | yellow=cell boundary)"),
        (mask_rgb,  f"Cell masks  ({n_cells} cells)"),
        (fitc_rgb,  "FITC intensity inside cells"),
    ]
    for ax, (im, title) in zip(axes, panels):
        ax.imshow(im, interpolation="nearest")
        ax.set_title(title, color="#cccccc", fontsize=10, pad=6)
        ax.axis("off")

    sm = plt.cm.ScalarMappable(cmap="hot", norm=mcolors.Normalize(0, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes[2], fraction=0.035, pad=0.02)
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white", fontsize=8)
    cbar.set_label("norm. intensity", color="white", fontsize=8)

    plt.tight_layout()
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight", facecolor="#0e0e0e")
    plt.close(fig)
    print(f"    Saved -> {out_path.name}")


# ─── main pipeline ────────────────────────────────────────────────────────────

def run_pipeline(folder: str) -> None:
    folder_path = Path(folder).resolve()
    if not folder_path.exists():
        sys.exit(f"[ERROR] Folder not found: {folder_path}")

    tif_files = sorted(
        list(folder_path.glob("*.tif")) + list(folder_path.glob("*.tiff"))
    )
    if not tif_files:
        sys.exit(f"[ERROR] No .tif/.tiff files found in: {folder_path}")

    print(f"\n{'='*62}")
    print(f"  Macropinocytosis Analysis Pipeline")
    print(f"  Folder : {folder_path}")
    print(f"  Images : {len(tif_files)}")
    print(f"  Method : DAPI nuclei → actin-bounded expansion")
    print(f"{'='*62}\n")

    results_dir      = folder_path / RESULTS_DIR_NAME
    side_by_side_dir = results_dir / "side_by_side"
    side_by_side_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []

    for idx, tif_path in enumerate(tif_files, start=1):
        print(f"[{idx}/{len(tif_files)}] {tif_path.name}")

        try:
            actin, fitc, dapi, rgb = load_tif(tif_path)
        except Exception as exc:
            print(f"    WARNING: Could not load — {exc}  (skipping)\n")
            continue

        print("    Segmenting ...")
        nuclear_labels, cell_labels = segment_cells(actin, dapi)
        n_cells = int(cell_labels.max())
        print(f"    {n_cells} cell(s) detected.")

        if n_cells == 0:
            print("    WARNING: No cells found.\n")
            summary_rows.append({
                "image":              tif_path.name,
                "n_cells":            0,
                "mean_intensity_all": float("nan"),
            })
            continue

        cell_df  = measure_intensity(fitc, cell_labels)
        mean_all = float(cell_df["mean_intensity"].mean())
        print(f"    Mean FITC intensity across all cells: {mean_all:.4f}")

        # side-by-side panel
        panel_path = side_by_side_dir / (tif_path.stem + "_analysis.png")
        make_side_by_side(rgb, fitc, nuclear_labels, cell_labels,
                          cell_df, tif_path.name, panel_path)

        # per-cell CSV
        per_cell_csv = results_dir / (tif_path.stem + "_per_cell.csv")
        cell_df.insert(0, "image", tif_path.name)
        cell_df.to_csv(str(per_cell_csv), index=False)

        summary_rows.append({
            "image":              tif_path.name,
            "n_cells":            n_cells,
            "mean_intensity_all": round(mean_all, 4),
        })
        print()

    # ── summary tables ───────────────────────────────────────────────────────
    if not summary_rows:
        print("[ERROR] No images were successfully processed.")
        return

    summary_df = pd.DataFrame(summary_rows)
    csv_path   = results_dir / "summary_table.csv"
    xlsx_path  = results_dir / "summary_table.xlsx"
    summary_df.to_csv(str(csv_path),  index=False)
    summary_df.to_excel(str(xlsx_path), index=False)

    print("Summary tables written:")
    print(f"  {csv_path}")
    print(f"  {xlsx_path}")
    print()
    print(summary_df.to_string(index=False))
    print(f"\nDone.  All outputs in: {results_dir}\n")


# ─── entry point ──────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Macropinocytosis pipeline: DAPI segmentation + "
                    "actin-bounded cell expansion + FITC intensity.")
    parser.add_argument(
        "folder",
        help="Path to folder containing .tif / .tiff microscopy images.")
    args = parser.parse_args()
    run_pipeline(args.folder)
