import pathlib
import subprocess
import joblib


filepaths = (
    open(
        r"/n/scratch/users/y/yc00000/152-152-Kidney_Imaging_ProspectiveCasesRound2-2/YC-20241105-kidney-segmentation/script-processing/files-ashlar-o2.csv"
    )
    .read()
    .replace('"', "")
    .strip()
    .split("\n")
)
filepaths = {pathlib.Path(ff).name.split("-")[0]: ff for ff in filepaths}

slides = list(filepaths.keys())
subset = slides[:]

mask_dirs = [
    r"/n/scratch/users/y/yc00000/152-152-Kidney_Imaging_ProspectiveCasesRound2-2/YC-20241105-kidney-segmentation/tubule-mask",
    r"/n/scratch/users/y/yc00000/152-152-Kidney_Imaging_ProspectiveCasesRound2-2/YC-20241105-kidney-segmentation/tubule-ring-mask",
    r"/n/scratch/users/y/yc00000/152-152-Kidney_Imaging_ProspectiveCasesRound2-2/YC-20241105-kidney-segmentation/tubule-epithelial-mask",
    r"/n/scratch/users/y/yc00000/152-152-Kidney_Imaging_ProspectiveCasesRound2-2/YC-20241105-kidney-segmentation/interstitial-mask",
]
mask_dirs = [pathlib.Path(mm) for mm in mask_dirs]

masks = {
    ii: [sorted(dd.glob(f"{ii}*mask*.ome.tif")) for dd in mask_dirs]
    for ii in filepaths.keys()
}


def run_image(slide_id):
    img_path = (
        pathlib.Path(
            r"/n/scratch/users/y/yc00000/152-152-Kidney_Imaging_ProspectiveCasesRound2-2/Ashlar"
        )
        / f"{slide_id}-ashlar.ome.tif"
    )
    _masks = masks[slide_id]
    masks_flat = []
    for mm in _masks:
        masks_flat.extend(mm)
    masks_flat = [str(pp) for pp in masks_flat]
    assert pathlib.Path(img_path).exists()
    assert len(masks_flat) == 9
    cmd = [
        "/home/yc00000/mambaforge/envs/ashlar/bin/python",
        r"script-processing/quantification/CommandSingleCellExtraction.py",
        "--image",
        str(img_path),
        "--output",
        "quantification",
        "--channel_names",
        r"script-processing/markers.csv",
        "--intensity_props",
        "intensity_mean_top_p",
        "--masks",
    ]
    cmd += masks_flat
    result = subprocess.run(
        cmd,
        cwd=r"/n/scratch/users/y/yc00000/152-152-Kidney_Imaging_ProspectiveCasesRound2-2/YC-20241105-kidney-segmentation",
    )
    return result.returncode


def run_gated_image(slide_id):
    img_path = (
        pathlib.Path(
            r"/n/scratch/users/y/yc00000/152-152-Kidney_Imaging_ProspectiveCasesRound2-2/Ashlar_gated"
        )
        / f"{slide_id}-ashlar-gated.ome.tif"
    )
    _masks = masks[slide_id]
    masks_flat = []
    for mm in _masks:
        masks_flat.extend(mm)
    masks_flat = [str(pp) for pp in masks_flat]
    assert pathlib.Path(img_path).exists()
    assert len(masks_flat) == 9
    cmd = [
        "/home/yc00000/mambaforge/envs/ashlar/bin/python",
        r"script-processing/quantification/CommandSingleCellExtraction.py",
        "--image",
        str(img_path),
        "--output",
        "quantification-gated",
        "--channel_names",
        r"script-processing/markers.csv",
        "--intensity_props",
        "intensity_positive_mean",
        "intensity_positive_median",
        "intensity_positive_count",
        "--masks",
    ]
    cmd += masks_flat
    result = subprocess.run(
        cmd,
        cwd=r"/n/scratch/users/y/yc00000/152-152-Kidney_Imaging_ProspectiveCasesRound2-2/YC-20241105-kidney-segmentation",
    )
    return result.returncode


def run_type_mask(slide_id):
    img_path = (
        pathlib.Path(
            r"/n/scratch/users/y/yc00000/152-152-Kidney_Imaging_ProspectiveCasesRound2-2/YC-20241105-kidney-segmentation"
        )
        / "tubule-type-mask"
        / f"{slide_id}-tubule_type_mask.ome.tif"
    )
    mask_path = (
        pathlib.Path(
            r"/n/scratch/users/y/yc00000/152-152-Kidney_Imaging_ProspectiveCasesRound2-2/YC-20241105-kidney-segmentation"
        )
        / "tubule-mask"
        / f"{slide_id}-tubule_mask.ome.tif"
    )
    assert img_path.exists()
    assert mask_path.exists()
    cmd = [
        "/home/yc00000/mambaforge/envs/ashlar/bin/python",
        r"script-processing/quantification/CommandSingleCellExtraction.py",
        "--image",
        str(img_path),
        "--output",
        "quantification-type",
        "--channel_names",
        r"script-processing/markers-object-type.csv",
        "--masks",
        str(mask_path),
    ]
    result = subprocess.run(
        cmd,
        cwd=r"/n/scratch/users/y/yc00000/152-152-Kidney_Imaging_ProspectiveCasesRound2-2/YC-20241105-kidney-segmentation",
    )
    return result.returncode


# joblib.Parallel(n_jobs=4, verbose=1)(
#     joblib.delayed(parallel_wrap)(filepaths[nn], masks[nn]) for nn in subset
# )


if __name__ == "__main__":
    import sys
    import fire
    import numpy as np

    def fire_wrap(slide_id):
        r0 = run_image(slide_id)
        r1 = run_gated_image(slide_id)
        r2 = run_type_mask(slide_id)
        print([r0, r1, r2])
        if np.all(np.array([r0, r1, r2]) == 0):
            return 0
        return 1

    sys.exit(fire.Fire(fire_wrap))


def minify_csv_gated(csv_path):
    import pandas as pd

    df = pd.read_csv(csv_path, index_col="CellID")
    cols = df.loc[:, "Hoechst":"MinorAxisLength"].columns
    df[cols] = df[cols].round().astype("int")
    df.loc[:, "Eccentricity":] = df.loc[:, "Eccentricity":].round(decimals=4)

    df.rename(
        columns=lambda x: x.replace("intensity_positive_mean", "pMean")
        .replace("intensity_positive_median", "pMedian")
        .replace("intensity_positive_count", "pCount"),
        inplace=True,
    )
    df.to_csv(csv_path)


def round_csv(csv_path):
    import pandas as pd

    df = pd.read_csv(csv_path, index_col="CellID")
    cols = df.loc[:, "Hoechst":"MinorAxisLength"].columns
    df[cols] = df[cols].round().astype("int")
    df.loc[:, "Eccentricity":] = df.loc[:, "Eccentricity":].round(decimals=4)

    ordered_cols = "Hoechst,AF,C3,PODXL,Albumin,Lambda,Fibrinogen,C1q,Nephrin,IgG,IgM,Synaptopodin,KIM1,Kappa,AQP1,IgA,ColIII,ColIV,SMA,Hoechst_intensity_mean_top_p-0,AF_intensity_mean_top_p-0,C3_intensity_mean_top_p-0,PODXL_intensity_mean_top_p-0,Albumin_intensity_mean_top_p-0,Lambda_intensity_mean_top_p-0,Fibrinogen_intensity_mean_top_p-0,C1q_intensity_mean_top_p-0,Nephrin_intensity_mean_top_p-0,IgG_intensity_mean_top_p-0,IgM_intensity_mean_top_p-0,Synaptopodin_intensity_mean_top_p-0,KIM1_intensity_mean_top_p-0,Kappa_intensity_mean_top_p-0,AQP1_intensity_mean_top_p-0,IgA_intensity_mean_top_p-0,ColIII_intensity_mean_top_p-0,ColIV_intensity_mean_top_p-0,SMA_intensity_mean_top_p-0,Hoechst_intensity_mean_top_p-1,AF_intensity_mean_top_p-1,C3_intensity_mean_top_p-1,PODXL_intensity_mean_top_p-1,Albumin_intensity_mean_top_p-1,Lambda_intensity_mean_top_p-1,Fibrinogen_intensity_mean_top_p-1,C1q_intensity_mean_top_p-1,Nephrin_intensity_mean_top_p-1,IgG_intensity_mean_top_p-1,IgM_intensity_mean_top_p-1,Synaptopodin_intensity_mean_top_p-1,KIM1_intensity_mean_top_p-1,Kappa_intensity_mean_top_p-1,AQP1_intensity_mean_top_p-1,IgA_intensity_mean_top_p-1,ColIII_intensity_mean_top_p-1,ColIV_intensity_mean_top_p-1,SMA_intensity_mean_top_p-1,Hoechst_intensity_mean_top_p-2,AF_intensity_mean_top_p-2,C3_intensity_mean_top_p-2,PODXL_intensity_mean_top_p-2,Albumin_intensity_mean_top_p-2,Lambda_intensity_mean_top_p-2,Fibrinogen_intensity_mean_top_p-2,C1q_intensity_mean_top_p-2,Nephrin_intensity_mean_top_p-2,IgG_intensity_mean_top_p-2,IgM_intensity_mean_top_p-2,Synaptopodin_intensity_mean_top_p-2,KIM1_intensity_mean_top_p-2,Kappa_intensity_mean_top_p-2,AQP1_intensity_mean_top_p-2,IgA_intensity_mean_top_p-2,ColIII_intensity_mean_top_p-2,ColIV_intensity_mean_top_p-2,SMA_intensity_mean_top_p-2,Hoechst_intensity_mean_top_p-3,AF_intensity_mean_top_p-3,C3_intensity_mean_top_p-3,PODXL_intensity_mean_top_p-3,Albumin_intensity_mean_top_p-3,Lambda_intensity_mean_top_p-3,Fibrinogen_intensity_mean_top_p-3,C1q_intensity_mean_top_p-3,Nephrin_intensity_mean_top_p-3,IgG_intensity_mean_top_p-3,IgM_intensity_mean_top_p-3,Synaptopodin_intensity_mean_top_p-3,KIM1_intensity_mean_top_p-3,Kappa_intensity_mean_top_p-3,AQP1_intensity_mean_top_p-3,IgA_intensity_mean_top_p-3,ColIII_intensity_mean_top_p-3,ColIV_intensity_mean_top_p-3,SMA_intensity_mean_top_p-3,X_centroid,Y_centroid,Area,MajorAxisLength,MinorAxisLength,Eccentricity,Solidity,Extent,Orientation".split(
        ","
    )

    df.loc[:, ordered_cols].to_csv(csv_path)
    return


def run_round_csv():
    import pathlib
    import tqdm

    csvs = sorted(
        pathlib.Path(
            r"Z:\yc00000\computation\YC-20241105-kidney-segmentation\quantification"
        ).glob("*ashlar*.csv")
    )

    len(csvs)
    for pp in tqdm.tqdm(csvs):
        round_csv(pp)


def assign_object_type():
    import pathlib
    import tqdm
    import pandas as pd
    import numpy as np

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
    slides = list(filepaths.keys())

    # quant_dir = pathlib.Path(
    #     
    # )
    quant_dir = pathlib.Path(
        r"Z:\yc00000\computation\YC-20241105-kidney-segmentation\quantification"
        # r"Z:\yc00000\computation\YC-20241105-kidney-segmentation\quantification-gated"
    )
    quant_type_dir = pathlib.Path(
        r"Z:\yc00000\computation\YC-20241105-kidney-segmentation\quantification-type"
    )

    csvs = sorted(quant_dir.glob("*interstitial*.csv"))
    for pp in tqdm.tqdm(csvs):
        df = pd.read_csv(pp, index_col="CellID")
        df["Type"] = "IN"
        df.to_csv(pp)

    # record object type
    # {1: tubule, 2: blood vessel, 3: glomerulus, 4: sclerosed glomerulus}
    for ss in tqdm.tqdm(slides[174:175]):
        type_path = next(quant_type_dir.glob(f"{ss}*type*.csv"))
        tub_path = sorted(quant_dir.glob(f"{ss}*tubule*.csv"))
        tub_path = list(filter(lambda x: "type" not in x.name, tub_path))

        tdf = pd.read_csv(type_path, usecols=["CellID", "Type"], index_col="CellID")
        object_type = np.full(tdf["Type"].shape, "TU")
        object_type[tdf["Type"] == 2] = "BV"
        object_type[tdf["Type"] == 3] = "GL"
        object_type[tdf["Type"] == 4] = "GS"
        tseries = pd.Series(object_type, index=tdf.index)

        for pp in tub_path:
            df = pd.read_csv(pp, index_col="CellID")
            df["Type"] = tseries[df.index]
            df.to_csv(pp)


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
