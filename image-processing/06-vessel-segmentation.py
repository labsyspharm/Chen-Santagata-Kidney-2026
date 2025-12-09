import napari
import pathlib
import palom
import numpy as np
import skimage.morphology
import scipy.ndimage as ndi
import pandas as pd


def mapping_indexer(df, column_name=None, value_for_missing_key=0):
    if column_name is None:
        assert df.shape[1] == 1
        column_name = df.columns[0]
    indexer = np.full(df.index.max() + 1, fill_value=value_for_missing_key)
    indexer[df.index.values] = df[column_name].values
    return indexer


def map_mask_labels(mask_pyramid, df_mapper):
    def recolor(mask, indexer):
        return indexer[mask].astype(indexer.dtype)

    mapped_masks = []
    for kk in df_mapper:
        idxr = mapping_indexer(df_mapper[[kk]], column_name=kk).astype(
            df_mapper[kk].dtype
        )
        mapped_masks.append(
            [
                ll[0]
                .map_blocks(contour, thickness=2, dtype=ll.dtype)
                .map_blocks(recolor, indexer=idxr, dtype=idxr.dtype)
                for ll in mask_pyramid
            ]
        )

    return mapped_masks


def contour(tmask, thickness=1):
    assert thickness in range(1, 10)
    contour = ndi.grey_dilation(
        tmask, footprint=skimage.morphology.disk(thickness)
    ) != ndi.grey_erosion(tmask, footprint=skimage.morphology.disk(thickness))
    return contour * tmask


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
subset = slides[1::150][:]

out_mask_dir = pathlib.Path(
    r"Z:\yc00000\computation\YC-20241105-kidney-segmentation\tubule-mask"
)
masks = {ii: next(out_mask_dir.glob(f"{ii}-tubule_mask.ome.tif")) for ii in subset}
quant_dir = pathlib.Path(
    r"Z:\yc00000\computation\YC-20241105-kidney-segmentation\quantification"
)

quant_tubule_paths = {
    ii: next(quant_dir.glob(f"{ii}-*-tubule_ring_mask_30.csv"))
    for ii in filepaths.keys()
}
# out_ring_mask_dir = pathlib.Path(
#     r"Z:\yc00000\computation\YC-20241105-kidney-segmentation\tubule-ring-mask"
# )
# ring_masks = {
#     ii: next(out_ring_mask_dir.glob(f"{ii}-tubule_ring_mask_30.ome.tif"))
#     for ii in subset
# }

v = napari.Viewer()
mask_readers = [palom.reader.OmePyramidReader(masks[nn]) for nn in subset]
# ring_mask_readers = [palom.reader.OmePyramidReader(ring_masks[nn]) for nn in subset]
img_readers = [palom.reader.OmePyramidReader(filepaths[nn]) for nn in subset]
quant_paths = [
    pd.read_csv(quant_tubule_paths[nn], index_col="CellID", engine="pyarrow")
    for nn in subset
]
# for mm, rr, mr in zip(mask_readers, img_readers, ring_mask_readers):
for mm, rr, qq in zip(mask_readers, img_readers, quant_paths):
    # for mm, rr in zip(mask_readers, img_readers):
    # for rr in img_readers:
    v.add_image(
        [np.log1p(pp[[1, 18]]) for pp in rr.pyramid],
        name=rr.path.name.split("-")[0],
        visible=False,
        contrast_limits=(np.log1p(200), np.log1p(5000)),
        channel_axis=0,
    )
    mask_pyramid = map_mask_labels(
        mm.pyramid, qq[["AF_intensity_mean_top_p-3", "AF_intensity_mean_top_p-2", "AF_intensity_mean_top_p-1"]].transform(np.log1p)
    )
    for mmm in mask_pyramid:
        v.add_image(
            mmm,
            name=rr.path.name.split("-")[0],
            visible=False,
            colormap="cividis",
            contrast_limits=np.log1p([100, 2000])
        )
    # v.add_labels(mm.pyramid, visible=False, name="mask")
    # v.add_labels(mr.pyramid, visible=False, name="ring mask")
