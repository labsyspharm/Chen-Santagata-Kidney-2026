import tifffile
import zarr
import dask.array as da
import skimage.feature
import pathlib
import tqdm.contrib
import pandas as pd
import numpy as np


before_path = r"X:\cycif-production\152-Kidney_Imaging\24-01_JY_ArgoPaletteFiles_Kidney_V4\Kidney_V4\Kidney_v4_515_C3_2023-12-21-03-32\Kidney_v4_515_C3_2023-12-21-03-32.pysed.ome.tifp"
after_path = r"X:\cycif-production\152-Kidney_Imaging\24-01_JY_ArgoPaletteFiles_Kidney_V4\250609_Re-processed_Artemis41\Kidney_v4_515_C3_2023-12-21-03-32\HMS_41_Reprocess_000002\Kidney_v4_515_C3_2023-12-21-03-32_HMS_41_Reprocess_000002.pysed.ome.tif"
out_dir = r"X:\cycif-production\152-Kidney_Imaging\24-01_JY_ArgoPaletteFiles_Kidney_V4\YC-plot-extraction"
channel = 2


def process_one_scan(before_path, after_path, channel, out_dir, min_distance=50):
    before = da.from_zarr(
        zarr.open(tifffile.imread(before_path, aszarr=True)), name=False
    ).reshape(-1, 19, 2000, 2000)

    after = da.from_zarr(
        zarr.open(tifffile.imread(after_path, aszarr=True)), name=False
    ).reshape(-1, 19, 2000, 2000)

    values_before = []
    values_after = []
    for tb, ta in tqdm.contrib.tzip(before, after):
        locations = skimage.feature.peak_local_max(
            tb[channel].compute(), min_distance=min_distance, num_peaks=10
        )
        values_before.append([cc.compute()[tuple(locations.T)] for cc in tb])
        values_after.append([cc.compute()[tuple(locations.T)] for cc in ta])

    def to_frame(intensities, norm_channel):
        df = pd.DataFrame(
            np.log1p(np.hstack(intensities).astype("float"))
            / np.log1p(np.hstack(intensities)[norm_channel]),
            index=[f"ch-{cc:02}" for cc in range(19)],
        )
        df.index.name = "Channel"
        return df.melt(var_name="Location", value_name="Intensity", ignore_index=False)

    out_dir = pathlib.Path(out_dir)
    to_frame(values_before, channel).to_csv(
        out_dir / f"channel-{channel:02}-0_before_extraction.csv"
    )
    to_frame(values_after, channel).to_csv(
        out_dir / f"channel-{channel:02}-1_after_extraction.csv"
    )


processed_dir = pathlib.Path(
    r"X:\cycif-production\152-Kidney_Imaging\24-01_JY_ArgoPaletteFiles_Kidney_V4\250609_Re-processed_Artemis41"
)
ori_dir = pathlib.Path(
    r"X:\cycif-production\152-Kidney_Imaging\24-01_JY_ArgoPaletteFiles_Kidney_V4\Kidney_V4"
)

scans = """
Kidney_v4_Hoechst_2024-01-05-12-11
Kidney_v4_AF_2024-01-05-11-59
Kidney_v4_515_C3_2023-12-21-03-32
Kidney_v4_555L_PODXL_2023-12-21-03-39
Kidney_v4_535_Albumin_2023-12-21-03-45
Kidney_v4_580L_Lambda_2023-12-21-03-51
Kidney_v4_660L_Fibrinogen_2024-01-05-10-03
Kidney_v4_572_C1q_2023-12-21-03-59
Kidney_v4_602_Nephrin_2023-12-22-10-07
Kidney_v4_624_IgG_2023-12-22-10-15
Kidney_v4_662_IgM_2023-12-22-10-31
Kidney_v4_686_Synaptopodi_2023-12-22-10-42
Kidney_v4_730_KIM1_2023-12-21-12-32
Kidney_v4_760_Kappa_2023-12-21-10-32
Kidney_v4_784_AQP1_2024-01-05-09-44
Kidney_v4_810_IgA_2023-12-21-01-32
LSP20339_Argo845_ColIII_2024-01-25-11-01
KidneyV4_874_ColIV_2024-04-05-01-15
KidneyV4_890_SMABio_001
""".strip().split("\n")

bpaths = [next((ori_dir / ss).glob("*.pysed.ome.tifp")) for ss in scans]
apaths = [
    next((processed_dir / ss).glob("HMS_41_Reprocess*/*.pysed.ome.tif")) for ss in scans
]
channels = list(range(len(bpaths)))
min_distances = [50] * len(bpaths)
min_distances[-1] = 5

out_dir = (
    r"X:\cycif-production\152-Kidney_Imaging\24-01_JY_ArgoPaletteFiles_Kidney_V4\test"
)

import joblib


joblib.Parallel(n_jobs=4, verbose=1)(
    joblib.delayed(process_one_scan)(bb, aa, cc, out_dir, mm)
    for bb, aa, cc, mm in zip(bpaths[:], apaths[:], channels[:], min_distances[:])
)


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


curr = pathlib.Path(
    r"X:\cycif-production\152-Kidney_Imaging\24-01_JY_ArgoPaletteFiles_Kidney_V4\YC-plot-extraction"
)

before_csvs = sorted(curr.glob("*_before*.csv"))
after_csvs = sorted(curr.glob("*_after*.csv"))


yticklabels = "01-DNA-Hoechst,02-AF-AF,03-C3-FITC,04-PODXL-Argo555L,05-Albumin-Argo535,06-Lambda-Argo580L,07-Fibrinogen-Argo660L,08-C1q-Argo572,09-Nephrin-Argo602,10-IgG-Argo624,11-IgM-Argo662,12-Synaptopodin-Argo686,13-KIM1-Argo730,14-Kappa-Argo760,15-AQP1-Argo784,16-IgA-Argo810,17-ColIII-Argo845,18-ColIV-Argo874,19-SMA-Argo890".strip().split(
    ","
)
xticklabels = [f"Channel {cc + 1}" for cc in range(19)]

_mask = pd.read_csv(
    r"X:\cycif-production\152-Kidney_Imaging\24-01_JY_ArgoPaletteFiles_Kidney_V4\250609_Re-processed_Artemis41\Kidney_v4_515_C3_2023-12-21-03-32\HMS_41_Reprocess_000002\HMS_Kidney_v4_reprocess.csv",
    index_col=0,
).values == 0
fig, axs = plt.subplots(1, 2)
for pp, ax in zip([before_csvs, after_csvs], axs):
    df_median = pd.concat(
        [
            pd.read_csv(ppp, usecols=["Channel", "Intensity"])
            .groupby("Channel")
            .median()
            .assign(Donor=ii)
            for ii, ppp in enumerate(pp)
        ]
    )
    sns.heatmap(
        df_median.reset_index().pivot(index="Donor", columns="Channel").iloc[2:, 2:],
        ax=ax,
        vmin=0,
        vmax=0.8,
        cmap="cividis",
        linewidths=0.1,
        xticklabels=xticklabels[2:],
        yticklabels=yticklabels[2:],
        square=True,
        mask=_mask[2:, 2:]
    )
    ax.set_facecolor('k')
