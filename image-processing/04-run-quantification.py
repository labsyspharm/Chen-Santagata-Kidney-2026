import pathlib
import subprocess
import joblib


filepaths = (
    open(
        r"Z:\yc00000\computation\YC-20241105-kidney-segmentation\script-processing\files-ashlar.csv"
    )
    .read()
    .replace('"', "")
    .strip()
    .split("\n")
)
filepaths = {pathlib.Path(ff).name.split("-")[0]: ff for ff in filepaths}
out_mask_dir = pathlib.Path(
    r"Z:\yc00000\computation\YC-20241105-kidney-segmentation\tubule-mask"
)
masks = {
    ii: next(out_mask_dir.glob(f"{ii}-tubule_mask.ome.tif")) for ii in filepaths.keys()
}
out_ring_mask_dir = pathlib.Path(
    r"Z:\yc00000\computation\YC-20241105-kidney-segmentation\tubule-ring-mask"
)
ring_masks = {
    ii: next(out_ring_mask_dir.glob(f"{ii}-tubule_ring_mask_30.ome.tif"))
    for ii in filepaths.keys()
}
out_interstitial_mask_dir = pathlib.Path(
    r"Z:\yc00000\computation\YC-20241105-kidney-segmentation\interstitial-mask"
)
interstitial_masks = {
    ii: next(out_interstitial_mask_dir.glob(f"{ii}-interstitial_mask.ome.tif"))
    for ii in filepaths.keys()
}

slides = list(filepaths.keys())
subset = slides[:]


def parallel_wrap(ii, mt, mr, mi):
    subprocess.run(
        [
            "python",
            r"script-processing\quantification\CommandSingleCellExtraction.py",
            "--image",
            ii,
            "--output",
            "quantification",
            "--channel_names",
            r"script-processing\markers.csv",
            "--intensity_props",
            "intensity_mean_top_p",
            "--masks",
            mt,
            mr,
            mi,
        ],
        cwd=r"Z:\yc00000\computation\YC-20241105-kidney-segmentation",
    )
    return


joblib.Parallel(n_jobs=4, verbose=1)(
    joblib.delayed(parallel_wrap)(
        filepaths[nn], masks[nn], ring_masks[nn], interstitial_masks[nn]
    )
    for nn in subset
)


def round_csv(csv_path):
    import pandas as pd

    df = pd.read_csv(csv_path, index_col="CellID")
    cols = df.loc[:, "Hoechst":"MinorAxisLength"].columns
    df[cols] = df[cols].round().astype("int")
    df.loc[:, "Eccentricity":] = df.loc[:, "Eccentricity":].round(decimals=4)
    df.to_csv(csv_path)
    return


def run_round_csv():
    import pathlib
    import tqdm

    csvs = sorted(
        pathlib.Path(
            r"Z:\yc00000\computation\YC-20241105-kidney-segmentation\quantification"
        ).glob("*.csv")
    )

    len(csvs)
    for pp in tqdm.tqdm(csvs):
        round_csv(pp)


"""
python script-processing\quantification\CommandSingleCellExtraction.py
    --masks
    tubule-mask\LSP22525-tubule_mask.ome.tif
    tubule-ring-mask\LSP22525-tubule_ring_mask_30.ome.tif
    interstitial-mask\LSP22525-interstitial_mask.ome.tif
    --image
    "X:\cycif-production\152-Kidney_Imaging_ProspectiveCasesRound2-2\Ashlar\LSP22525-ashlar.ome.tif"
    --output
    quantification
    --channel_names
    script-processing\markers.csv
    --intensity_props
    intensity_mean_top_p

"""
