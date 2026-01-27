# Kidney Orion Image Processing Scripts

This repository contains Python scripts for image processing and quantification
related to kidney tissue analysis. The workflow involves segmentation of
different tissue components (tubules, interstitial regions), followed by
quantification of cellular features within these segments and assignment of
object types.

## Scripts Overview

### `00-tubule-segmentation.py`

This script focuses on segmenting tubule structures within kidney images. It
performs the following key operations:

- **ROI Parsing and Mask Generation**: Reads ROI data (e.g., from CSV files) and
  converts them into image masks. It handles various ROI shapes like polygons
  and ellipses.
- **Glomerulus Segmentation**: Identifies and segments glomeruli (both normal
  and sclerosed) based on ROI data.
- **Vessel Clearance**: Excludes blood vessels from the tubule segmentation to
  avoid misclassification.
- **Tubule and Object Type Mask Generation**: Generates masks for tubules and
  assigns object types (tubule, blood vessel, glomerulus, sclerosed glomerulus)
  to segmented regions.
- **Derived Mask Generation**: Creates tubule ring masks (dilated tubules) and
  tubule epithelial masks (eroded tubules).

### `01-interstitial-segmentation.py`

This script is responsible for segmenting the interstitial regions of the kidney
tissue. It leverages outputs from the tubule segmentation step:

- **Interstitial Mask Generation**: Initial interstitial areas were acquired by
  applying intensity threshold using the Col III channel and excluding tubule,
  ring, and glomerulus regions. The initial interstitial areas are then refined
  by intersecting them with dilated glomeruli masks to produce the final
  interstitial masks.

### `04-run-quantification-o2.py`

This is the primary script for running quantification on the segmented images.
It orchestrates calls to an external `CommandSingleCellExtraction.py` script
from [this module](https://github.com/labsyspharm/quantification):

- **Image and Mask Loading**: Prepares image and mask paths for quantification.
- **Quantification Execution**: Runs single-object extraction for Orion images
  and tubule type masks.
- **CSV Processing**: Includes utility functions (`minify_csv_gated`,
  `round_csv`) for post-processing the generated quantification CSV files.

### `05-assign-object-type.py`

This script is used to refine the object type assignments in the quantification
data:

- **Type Assignment**: Assigns specific object types (e.g., "IN" for
  interstitial, "TU" for tubule, "BV" for blood vessel, "GL" for glomerulus,
  "GS" for sclerosed glomerulus) to segments based on the previously generated
  tubule type masks.
- **Quantification Data Update**: Modifies existing quantification CSV files
  with the assigned object types.

## Workflow

The typical workflow for processing images with these scripts would involve:

1. Run `00-tubule-segmentation.py` to generate tubule, ring, and type masks.
2. Run `01-interstitial-segmentation.py` to generate interstitial masks,
   utilizing the outputs from step 1.
3. Execute `04-run-quantification-o2.py` to perform single-object quantification
   using the generated masks.
4. Run `05-assign-object-type.py` to refine object type assignments in the
   quantification data.

## Usage

Each script is designed to be run independently or as part of a larger pipeline.
Refer to the individual script files for specific function parameters and
execution details.

## Dependencies

The scripts rely on several Python libraries, including:

- `numpy`
- `pandas`
- `scipy`
- `scikit-image`
- `opencv-python`
- `matplotlib`
- `shapely`
- `palom`
- `dask`
- `joblib`
- `tifffile`
- `zarr`
- `fire`

Please ensure these dependencies are installed in your environment before
running the scripts.
