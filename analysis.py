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
  4.  Measure raw dye intensity inside each cell region and report:
        - mean_intensity                 : mean pixel value within the cell mask
        - integrated_intensity           : sum of all pixel values within the cell mask
                                           (total uptake per cell; preferred metric when
                                           cell size varies between conditions)
        - integrated_intensity_per_area  : integrated_intensity / pixel_count
                                           (uptake density; corrects for cell-size bias)
        - pixel_count                    : cell area in px²
  5.  Save a side-by-side PNG (original composite | cell masks | dye heat-map)
      for every image, plus a per-cell CSV for every image and a summary
      CSV/XLSX across the whole folder.

  The folder-level summary reports, per image:
      n_cells, mean_integrated_intensity (mean across cells in that image),
      total_integrated_intensity (sum across all cells in that image)
  and a final TOTALS row across all images in the folder:
      total n_cells, grand mean integrated intensity per cell,
      grand total integrated intensity.

Channel file naming convention
-------------------------------
  Each "image" is represented by THREE separate single-channel TIF files that
  share a common base name and differ only in their channel suffix:

      <basename>_ch00.tif   →  actin  (cell boundary marker)   [default]
      <basename>_ch01.tif   →  dye / FITC / TMR  (signal)      [default]
      <basename>_ch02.tif   →  DAPI  (nuclear stain)           [default]

  Override which suffix maps to which role at the command line:

      python analyze_macropinocytosis.py /path/to/folder \
          --actin _ch00 --dye _ch01 --dapi _ch02

  The script discovers triplets by scanning for all files whose stem ends with
  the actin suffix, then locating the matching dye / dapi siblings.
  A metadata subfolder (or any file that does not match) is silently skipped.

  If your files have no shared base name, set TRIPLET_MODE = "sorted_order"
  and the script will group sorted TIFs as consecutive triples instead.

USAGE
-----
    python analyze_macropinocytosis.py /path/to/folder
    python analyze_macropinocytosis.py /path/to/folder --actin _ch02 --dapi _ch00 --dye _ch01

REQUIREMENTS
------------
    pip install tifffile imagecodecs numpy pandas matplotlib scikit-image scipy openpyxl

OUTPUTS  (written to <folder>/results/)
---------------------------------------
    side_by_side/<basename>_analysis.png  – 3-panel figure per image group
    <basename>_per_cell.csv               – per-cell data for every image group
    summary_table.csv / .xlsx             – per-image + TOTALS row summary
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from skimage.color import label2rgb
from skimage.filters import gaussian, threshold_otsu
from skimage.measure import label as sklabel, regionprops
from skimage.morphology import closing, disk
from skimage.segmentation import expand_labels, find_boundaries, watershed
from scipy.ndimage import distance_transform_edt


# ═══════════════════════════ TUNABLE PARAMETERS ═══════════════════════════════

# ── How channel TIFs are grouped ─────────────────────────────────────────────
# "suffix_match"  : look for files ending in the actin suffix, then find the
#                   dye / dapi siblings with the same base name.
# "sorted_order"  : group every sorted TIF as consecutive triples
#                   (files 0-2 → image 1, files 3-5 → image 2, …).
TRIPLET_MODE = "suffix_match"   # "suffix_match" | "sorted_order"

# Default channel suffixes — override at runtime with --actin / --dapi / --dye
DEFAULT_ACTIN_SUFFIX = "_ch02"
DEFAULT_DYE_SUFFIX   = "_ch01"
DEFAULT_DAPI_SUFFIX  = "_ch00"

# ── DAPI nucleus detection ───────────────────────────────────────────────────
DAPI_BLUR_SIGMA     = 4      # gaussian smoothing before threshold (pixels)
DAPI_CLOSE_RADIUS   = 8      # morphological closing to fill holes in nuclei
MIN_NUCLEUS_AREA    = 5000   # minimum nucleus area in px² (removes debris)
EXCLUDE_EDGE_NUCLEI = False  # if True, skip nuclei touching the image border
WATERSHED_MIN_DIST  = 20     # minimum distance in px between nucleus peaks for
                             # watershed splitting; increase to split less aggressively,
                             # decrease to split more.  Set to 0 to disable.

# ── Cell body expansion ───────────────────────────────────────────────────────
CELL_EXPANSION_PX   = 150    # max expansion radius from nucleus edge (pixels)
ACTIN_PERCENTILE    = 90     # actin pixels above this percentile define cell
                             # territory.  Lower = more permissive expansion.

# ── Output ────────────────────────────────────────────────────────────────────
RESULTS_DIR_NAME = "results"

# ══════════════════════════════════════════════════════════════════════════════


# ─── file discovery ───────────────────────────────────────────────────────────

def find_triplets(
    folder: Path,
    actin_suffix: str,
    dapi_suffix: str,
    dye_suffix: str,
) -> list[tuple[Path, Path, Path, str]]:
    """
    Return a list of (actin_path, dye_path, dapi_path, basename) tuples,
    one entry per logical image, sorted by basename.

    In suffix_match mode the actin suffix is used as the anchor; the dye
    and dapi siblings are located by replacing that suffix in the filename.
    Skips any files / subdirectories that do not match the expected pattern.
    """
    tif_exts = {".tif", ".tiff"}

    if TRIPLET_MODE == "sorted_order":
        all_tifs = sorted(
            p for p in folder.iterdir()
            if p.is_file() and p.suffix.lower() in tif_exts
        )
        if len(all_tifs) % 3 != 0:
            print(
                f"  WARNING: {len(all_tifs)} TIF files found — "
                "not a multiple of 3.  Trailing file(s) will be ignored."
            )
        triplets = []
        for i in range(0, len(all_tifs) - 2, 3):
            f_actin, f_dye, f_dapi = all_tifs[i], all_tifs[i + 1], all_tifs[i + 2]
            basename = f_actin.stem
            triplets.append((f_actin, f_dye, f_dapi, basename))
        return triplets

    # ── suffix_match (default) — anchor on the actin suffix ─────────────────
    actin_files = sorted(
        p for p in folder.iterdir()
        if p.is_file()
        and p.suffix.lower() in tif_exts
        and p.stem.lower().endswith(actin_suffix.lower())
    )

    triplets = []
    missing  = []

    for f_actin in actin_files:
        stem_base = re.sub(
            re.escape(actin_suffix) + r"$", "", f_actin.stem, flags=re.IGNORECASE
        )
        ext = f_actin.suffix

        def find_sibling(suffix: str) -> Path:
            p = f_actin.parent / (stem_base + suffix + ext)
            if not p.exists():
                alt = f_actin.parent / (stem_base + suffix + ext.upper())
                p = alt if alt.exists() else p
            return p

        f_dye  = find_sibling(dye_suffix)
        f_dapi = find_sibling(dapi_suffix)

        ok = True
        for p, label in [(f_dye, "dye (" + dye_suffix + ")"),
                         (f_dapi, "dapi (" + dapi_suffix + ")")]:
            if not p.exists():
                missing.append(
                    f"  MISSING {label}: {p.name}  (expected sibling of {f_actin.name})"
                )
                ok = False
        if ok:
            triplets.append((f_actin, f_dye, f_dapi, stem_base))

    if missing:
        print("\n".join(missing))

    return triplets


# ─── image loading ────────────────────────────────────────────────────────────

def load_channel(path: Path) -> np.ndarray:
    """
    Load a single-channel TIF and return a 2-D float32 array (H, W).

    Handles:
      (H, W)       grayscale — most single-channel exports
      (H, W, 1)    extra trailing dimension
      (1, H, W)    channel-first with one channel
      (Z, H, W)    Z-stack — max-projected along Z
      (Z, 1, H, W) Z-stack channel-first — max-projected
    """
    raw = tifffile.imread(str(path)).astype(np.float32)

    if raw.ndim == 2:
        return raw

    if raw.ndim == 3:
        if raw.shape[0] == 1:
            return raw[0]
        if raw.shape[2] == 1:
            return raw[..., 0]
        return raw.max(axis=0)   # assume (Z, H, W) Z-stack

    if raw.ndim == 4:
        projected = raw.max(axis=0)
        if projected.shape[0] == 1:
            return projected[0]
        if projected.shape[2] == 1:
            return projected[..., 0]
        return projected.max(axis=0)

    raise ValueError(
        f"Cannot interpret shape {raw.shape} as a single-channel image: {path.name}"
    )


def load_triplet(f_actin: Path, f_dye: Path, f_dapi: Path):
    """
    Load the three channel TIFs and return:
      actin  (float32, H×W)   – cell boundary marker
      dye    (float32, H×W)   – signal of interest
      dapi   (float32, H×W)   – nuclear stain
      rgb    (float32, H×W×3) – composite for display (actin/dye/dapi → R/G/B)
    """
    actin = load_channel(f_actin)
    dye   = load_channel(f_dye)
    dapi  = load_channel(f_dapi)

    shapes = {actin.shape, dye.shape, dapi.shape}
    if len(shapes) > 1:
        raise ValueError(
            f"Channel size mismatch:\n"
            f"  actin {actin.shape}  {f_actin.name}\n"
            f"  dye   {dye.shape}    {f_dye.name}\n"
            f"  dapi  {dapi.shape}   {f_dapi.name}"
        )

    rgb = np.stack([actin, dye, dapi], axis=-1)   # (H, W, 3) — R=actin G=dye B=dapi
    return actin, dye, dapi, rgb


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

    # ── Step 2: connected-component labelling + watershed splitting ─────────
    raw_labels = sklabel(nuclear_bin)
    if WATERSHED_MIN_DIST > 0:
        distance   = distance_transform_edt(nuclear_bin)
        from skimage.feature import peak_local_max
        peak_coords = peak_local_max(
            distance,
            min_distance=WATERSHED_MIN_DIST,
            labels=nuclear_bin,
        )
        markers = np.zeros_like(nuclear_bin, dtype=np.int32)
        for i, (r, c) in enumerate(peak_coords, start=1):
            markers[r, c] = i
        raw_labels = watershed(-distance, markers, mask=nuclear_bin)

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
    actin_smooth   = gaussian(actin, sigma=2)
    actin_thresh   = np.percentile(actin_smooth, ACTIN_PERCENTILE)
    cell_territory = (actin_smooth > actin_thresh) | (nuclear_labels > 0)

    # ── Step 5: expand nuclei outward, clipped to cell territory ─────────────
    cell_labels = expand_labels(nuclear_labels, distance=CELL_EXPANSION_PX)
    cell_labels[~cell_territory] = 0

    return nuclear_labels, cell_labels.astype(np.int32)


# ─── measurement ──────────────────────────────────────────────────────────────

def measure_intensity(
    dye: np.ndarray,
    cell_labels: np.ndarray,
) -> pd.DataFrame:
    """
    For each cell region, measure raw dye intensity (no background subtraction):
      - mean_intensity                : mean pixel value within the cell mask
      - integrated_intensity          : sum of all pixel values within the cell mask;
                                        preferred metric when cell size varies between
                                        conditions, as it captures total uptake per cell
      - integrated_intensity_per_area : integrated_intensity / pixel_count;
                                        uptake density — corrects for cell-size bias,
                                        equivalent to mean_intensity but reported
                                        explicitly alongside integrated for convenience
      - pixel_count                   : number of pixels in the cell region (area in px²)
    """
    records = []
    for cid in range(1, int(cell_labels.max()) + 1):
        px = dye[cell_labels == cid]
        if px.size == 0:
            continue
        integrated  = float(px.sum())
        pixel_count = int(px.size)
        records.append({
            "cell_id":                         cid,
            "mean_intensity":                  round(float(px.mean()),         4),
            "integrated_intensity":            round(integrated,                4),
            "integrated_intensity_per_area":   round(integrated / pixel_count,  4),
            "pixel_count":                     pixel_count,
        })
    return pd.DataFrame(records)


# ─── visualisation ────────────────────────────────────────────────────────────

def make_side_by_side(
    rgb:            np.ndarray,
    dye:            np.ndarray,
    nuclear_labels: np.ndarray,
    cell_labels:    np.ndarray,
    cell_df:        pd.DataFrame,
    image_name:     str,
    dye_suffix:     str,
    out_path:       Path,
) -> None:
    """
    Three-panel figure:
      1. Contrast-stretched RGB composite with nucleus (blue) and
         cell-boundary (yellow) outlines drawn on top.
      2. Cell regions, each uniquely coloured.
      3. Dye intensity heat-map restricted to cell regions.
    """
    composite = normalise_rgb(rgb)

    nuc_bounds  = find_boundaries(nuclear_labels, mode="outer")
    cell_bounds = find_boundaries(cell_labels,    mode="outer")
    comp_ann = composite.copy()
    comp_ann[cell_bounds] = [1.0, 0.85, 0.0]   # yellow = cell edge
    comp_ann[nuc_bounds]  = [0.2, 0.65, 1.0]   # blue   = nucleus edge

    mask_rgb = label2rgb(cell_labels, bg_label=0, bg_color=(0.08, 0.08, 0.08))

    dye_n      = normalise(dye)
    dye_masked = dye_n.copy()
    dye_masked[cell_labels == 0] = 0.0
    dye_rgb    = plt.cm.hot(dye_masked)[..., :3]

    dye_label  = dye_suffix.lstrip("_")   # e.g. "_ch01" → "ch01"
    n_cells    = len(cell_df)
    mean_integ = cell_df["integrated_intensity"].mean()

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), facecolor="#0e0e0e")
    fig.suptitle(
        f"{image_name}   |   {n_cells} cells   |   "
        f"mean integrated intensity: {mean_integ:.1f}",
        color="white", fontsize=13, fontweight="bold", y=1.01,
    )
    panels = [
        (comp_ann,  "Original composite  (blue=nucleus | yellow=cell boundary)"),
        (mask_rgb,  f"Cell masks  ({n_cells} cells)"),
        (dye_rgb,   f"Dye ({dye_label}) intensity inside cells"),
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

def run_pipeline(
    folder: str,
    actin_suffix: str,
    dapi_suffix: str,
    dye_suffix: str,
) -> None:
    folder_path = Path(folder).resolve()
    if not folder_path.exists():
        sys.exit(f"[ERROR] Folder not found: {folder_path}")

    print(f"\n{'='*62}")
    print(f"  Macropinocytosis Analysis Pipeline")
    print(f"  Folder      : {folder_path}")
    print(f"  Triplet mode: {TRIPLET_MODE}")
    print(f"  Channels    : actin={actin_suffix}  dye={dye_suffix}  dapi={dapi_suffix}")
    print(f"  Method      : DAPI nuclei → actin-bounded expansion")
    print(f"  Intensity   : raw pixel values, no background subtraction")
    print(f"  Primary metric: integrated intensity per cell")
    print(f"{'='*62}\n")

    triplets = find_triplets(folder_path, actin_suffix, dapi_suffix, dye_suffix)
    if not triplets:
        sys.exit(
            "[ERROR] No complete channel triplets found.\n"
            f"  Looking for actin suffix '{actin_suffix}' as anchor,\n"
            f"  with siblings '{dye_suffix}' (dye) and '{dapi_suffix}' (dapi)\n"
            f"  in: {folder_path}\n"
            f"  Tip: use --actin / --dye / --dapi to change which suffix maps to which role.\n"
            f"  Or set TRIPLET_MODE = 'sorted_order' if files have no shared base name."
        )

    print(f"  Found {len(triplets)} image group(s).\n")

    results_dir      = folder_path / RESULTS_DIR_NAME
    side_by_side_dir = results_dir / "side_by_side"
    side_by_side_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []

    for idx, (f_actin, f_dye, f_dapi, basename) in enumerate(triplets, start=1):
        print(f"[{idx}/{len(triplets)}] {basename}")
        print(f"    actin ({actin_suffix}) : {f_actin.name}")
        print(f"    dye   ({dye_suffix}) : {f_dye.name}")
        print(f"    dapi  ({dapi_suffix}) : {f_dapi.name}")

        try:
            actin, dye, dapi, rgb = load_triplet(f_actin, f_dye, f_dapi)
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
                "image":                      basename,
                "n_cells":                    0,
                "mean_integrated_intensity":  float("nan"),
                "total_integrated_intensity": float("nan"),
            })
            continue

        cell_df     = measure_intensity(dye, cell_labels)
        mean_integ  = float(cell_df["integrated_intensity"].mean())
        total_integ = float(cell_df["integrated_intensity"].sum())
        print(f"    Mean integrated intensity per cell : {mean_integ:.4f}")
        print(f"    Total integrated intensity (image) : {total_integ:.4f}")

        panel_path = side_by_side_dir / (basename + "_analysis.png")
        make_side_by_side(
            rgb, dye, nuclear_labels, cell_labels,
            cell_df, basename, dye_suffix, panel_path,
        )

        per_cell_csv = results_dir / (basename + "_per_cell.csv")
        cell_df.insert(0, "image", basename)
        cell_df.to_csv(str(per_cell_csv), index=False)

        summary_rows.append({
            "image":                      basename,
            "n_cells":                    n_cells,
            "mean_integrated_intensity":  round(mean_integ,  4),
            "total_integrated_intensity": round(total_integ, 4),
        })
        print()

    # ── summary tables ───────────────────────────────────────────────────────
    if not summary_rows:
        print("[ERROR] No images were successfully processed.")
        return

    summary_df = pd.DataFrame(summary_rows)

    # Folder-level totals row — grand mean is weighted by cell count per image
    valid             = summary_df.dropna(subset=["mean_integrated_intensity"])
    total_n_cells     = int(valid["n_cells"].sum())
    grand_mean_integ  = (
        float(valid["mean_integrated_intensity"].mul(valid["n_cells"]).sum() / total_n_cells)
        if total_n_cells > 0 else float("nan")
    )
    grand_total_integ = float(valid["total_integrated_intensity"].sum())

    totals_row = pd.DataFrame([{
        "image":                      "TOTALS",
        "n_cells":                    total_n_cells,
        "mean_integrated_intensity":  round(grand_mean_integ,  4),
        "total_integrated_intensity": round(grand_total_integ, 4),
    }])
    summary_with_totals = pd.concat([summary_df, totals_row], ignore_index=True)

    csv_path  = results_dir / "summary_table.csv"
    xlsx_path = results_dir / "summary_table.xlsx"
    summary_with_totals.to_csv(str(csv_path),  index=False)
    summary_with_totals.to_excel(str(xlsx_path), index=False)

    print("Summary tables written:")
    print(f"  {csv_path}")
    print(f"  {xlsx_path}")
    print()
    print(summary_with_totals.to_string(index=False))
    print(f"\nDone.  All outputs in: {results_dir}\n")


# ─── entry point ──────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Macropinocytosis pipeline: loads per-channel TIF files, "
            "segments cells from DAPI with actin-bounded expansion, "
            "and measures raw dye intensity (mean and integrated per cell)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "folder",
        help="Path to folder containing channel TIF files.",
    )
    parser.add_argument(
        "--actin", metavar="SUFFIX",
        default=DEFAULT_ACTIN_SUFFIX,
        help="Filename suffix (without extension) for the actin channel.",
    )
    parser.add_argument(
        "--dapi", metavar="SUFFIX",
        default=DEFAULT_DAPI_SUFFIX,
        help="Filename suffix (without extension) for the DAPI channel.",
    )
    parser.add_argument(
        "--dye", metavar="SUFFIX",
        default=DEFAULT_DYE_SUFFIX,
        help="Filename suffix (without extension) for the dye/FITC channel.",
    )
    args = parser.parse_args()
    run_pipeline(
        args.folder,
        actin_suffix=args.actin,
        dapi_suffix=args.dapi,
        dye_suffix=args.dye,
    )