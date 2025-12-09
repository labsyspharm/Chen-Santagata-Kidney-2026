import pathlib
import re

import cv2
import dask.array as da
import joblib
import matplotlib.patches as mpatches
import numpy as np
import palom
import pandas as pd
import scipy.spatial.distance as sdistance
import shapely
import shapely.validation
import skimage.filters


def parse_roi_points(all_points):
    return np.array(re.findall(r"-?\d+\.?\d+", all_points), dtype=float).reshape(-1, 2)


def ellipse_points_to_patch(
    vertex_1, vertex_2, co_vertex_1, co_vertex_2, patch_kwargs={}
):
    """
    Parameters
    ----------
    vertex_1, vertex_2, co_vertex_1, co_vertex_2: array like, in the form of (x-coordinate, y-coordinate)

    """
    v_and_co_v = np.array([vertex_1, vertex_2, co_vertex_1, co_vertex_2])
    centers = v_and_co_v.mean(axis=0)

    d = sdistance.cdist(v_and_co_v, v_and_co_v, metric="euclidean")
    width = d[0, 1]
    height = d[2, 3]

    vector_2 = v_and_co_v[1] - v_and_co_v[0]
    vector_2 /= np.linalg.norm(vector_2)

    angle = np.degrees(np.arccos([1, 0] @ vector_2))

    ellipse_patch = mpatches.Ellipse(
        centers, width=width, height=height, angle=angle, **patch_kwargs
    )
    return ellipse_patch


def add_mpatch(roi, patch_kwargs=None):
    if patch_kwargs is None:
        patch_kwargs = {}
    roi = roi.copy()
    points = parse_roi_points(roi["all_points"])
    roi.loc["parsed_points"] = points

    roi_type = roi["type"]
    if roi_type in ["Point", "Line"]:
        roi_mpatch = mpatches.Polygon(points, closed=False, **patch_kwargs)
    elif roi_type in ["Rectangle", "Polygon", "Polyline"]:
        roi_mpatch = mpatches.Polygon(points, closed=True, **patch_kwargs)
    elif roi_type == "Ellipse":
        roi_mpatch = ellipse_points_to_patch(*points, patch_kwargs=patch_kwargs)
    else:
        raise ValueError

    roi.loc["mpatch"] = roi_mpatch
    return roi


def glomerulus_roi_path_to_mask(roi_path, mask_shape=None, buffer_size=0):
    _roi = pd.read_csv(roi_path)
    is_polygon = _roi["type"].fillna("").astype("str").apply(lambda x: x == "Polyline")
    assert is_polygon.sum() <= 255
    roi = _roi.loc[is_polygon].apply(add_mpatch, axis=1)

    polygons = [
        to_valid_polygon(patch_to_shapely_polygon(pp).buffer(buffer_size))
        for pp in roi["mpatch"].dropna()
    ]
    mask = shapely_polygons_to_mask(polygons, mask_shape, 1)

    indexer = np.zeros(mask.max() + 1, dtype="uint8")
    for ii, tt in enumerate(roi["Text"].fillna("GL")):
        tt = tt.upper()
        assert tt in ["GL", "GS"]
        if tt == "GL":
            indexer[ii + 1] = 1
        else:
            indexer[ii + 1] = 2

    _mask = indexer[mask]
    # return (glomeruli-mask, sclerosed-glomeruli-mask)
    return (np.where(mask == 1, mask, 0), np.where(mask == 2, mask, 0))


def roi_path_to_mask(roi_path, mask_shape=None, fill_value_min=None, buffer_size=0):
    _roi = pd.read_csv(roi_path)
    is_polygon = _roi["type"].fillna("").astype("str").apply(lambda x: x == "Polyline")
    assert len(_roi) <= 255
    roi = _roi.loc[is_polygon].apply(add_mpatch, axis=1)

    polygons = [
        to_valid_polygon(patch_to_shapely_polygon(pp).buffer(buffer_size))
        for pp in roi["mpatch"].dropna()
    ]
    return shapely_polygons_to_mask(polygons, mask_shape, fill_value_min)


def shapely_polygons_to_mask(polygons, mask_shape=None, fill_value_min=None):
    coords = [np.array(pp.exterior.coords) for pp in polygons]

    if mask_shape is None:
        mask_shape = np.ceil(np.vstack(coords).max(axis=0)).astype("int")[::-1]
    h, w = mask_shape

    if fill_value_min is None:
        fill_value_min = 1
    fill_value = int(fill_value_min)

    mask = np.zeros((int(h), int(w)), np.uint8)
    for pp in coords:
        cv2.fillPoly(mask, [pp.round().astype("int")], fill_value)
        fill_value += 1
    return mask


def to_valid_polygon(polygon: shapely.Polygon):
    polygon = shapely.validation.make_valid(polygon)

    if isinstance(polygon, shapely.Polygon):
        return polygon

    elif isinstance(polygon, (shapely.MultiPolygon, shapely.GeometryCollection)):

        def flatten(lst):
            flat_list = []
            for item in lst:
                if isinstance(item, list):
                    flat_list.extend(flatten(item))
                else:
                    flat_list.append(item)
            return flat_list

        def unpack_geom(geometry):
            if not hasattr(geometry, "geoms"):
                return [geometry]
            return [unpack_geom(gg) for gg in geometry.geoms]

        geometries = flatten(unpack_geom(polygon))
        polygons = list(filter(lambda x: isinstance(x, shapely.Polygon), geometries))
        if len(polygons) == 0:
            return None
        areas = [pp.area for pp in polygons]
        idx = np.argmax(areas)
        return polygons[idx]
    else:
        return None


def patch_to_shapely_polygon(patch: mpatches.Patch):
    polygon = shapely.Polygon(patch.get_verts())
    return to_valid_polygon(polygon)


def filter_label_area(label_img, area_min, area_max):
    if np.all(label_img == 0):
        return label_img
    ids, counts = da.compute(
        *da.unique(
            da.from_array(label_img, chunks=1000, name=False),
            return_counts=True,
        )
    )
    filtered = np.where((counts >= area_min) & (counts < area_max), ids, 0)
    indexer = np.arange(ids.max() + 1)
    indexer[ids] = filtered
    return indexer[label_img].astype(np.int32)


def segment_interstitial(
    img_path, tubule_mask_path, tubule_ring_mask_path, tubule_type_mask_path
):
    Reader = palom.reader.OmePyramidReader
    i_reader = Reader(img_path)
    m_reader = Reader(tubule_mask_path)
    r_reader = Reader(tubule_ring_mask_path)
    t_reader = Reader(tubule_type_mask_path)

    img = i_reader.pyramid[2][16].compute()
    log_img = np.log1p(np.clip(img, 200, None))
    threshold = skimage.filters.threshold_yen(log_img[::5, ::5])

    _mask = cv2.GaussianBlur(log_img, (0, 0), 2) > threshold
    _, _mask = cv2.connectedComponents(_mask.astype("uint8"), connectivity=8)
    _mask = filter_label_area(_mask, 50 * 50, np.inf).astype("bool")

    interstitial_mask = np.full(i_reader.pyramid[0][0].shape, False, dtype="bool")
    _mask = palom.img_util.repeat_2d(_mask, (2**2,) * 2)
    h, w = np.min([_mask.shape, interstitial_mask.shape], axis=0)
    interstitial_mask[:h, :w] = _mask[:h, :w]

    # _, interstitial_mask = cv2.connectedComponents(
    #     interstitial_mask.astype(np.uint8), connectivity=8
    # )

    tubule_mask = m_reader.pyramid[0][0] + r_reader.pyramid[0][0]
    glomerulus_id_min = tubule_mask[t_reader.pyramid[0][0] >= 3].min()

    glomerulus_mask = np.array(tubule_mask >= glomerulus_id_min)
    _tissue_mask = np.array((tubule_mask > 0) & (tubule_mask < glomerulus_id_min))
    tubule_mask = np.array(tubule_mask.astype("bool"))

    _interstitial_mask = np.full_like(interstitial_mask, False)
    _interstitial_mask[:] = interstitial_mask[:]

    interstitial_mask[tubule_mask] = False
    _interstitial_mask[glomerulus_mask] = False
    _interstitial_mask[_tissue_mask] = True

    return interstitial_mask, _interstitial_mask


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
out_type_mask_dir = pathlib.Path(
    r"Z:\yc00000\computation\YC-20241105-kidney-segmentation\tubule-type-mask"
)
type_masks = {
    ii: next(out_type_mask_dir.glob(f"{ii}-tubule_type_mask.ome.tif"))
    for ii in filepaths.keys()
}
rois = {
    ii: next(
        pathlib.Path(
            r"X:\cycif-production\152-Kidney_Imaging_ProspectiveCasesRound2-2\ROI export"
        ).glob(f"{ii}*-rois.csv")
    )
    for ii in filepaths.keys()
}
slides = list(filepaths.keys())
subset = slides[:]

out_interstitial_mask_dir = pathlib.Path(
    r"Z:\yc00000\computation\YC-20241105-kidney-segmentation\interstitial-mask"
)
out_interstitial_mask_dir.mkdir(exist_ok=True, parents=True)
out_paths = {
    ii: out_interstitial_mask_dir / f"{ii}-interstitial_mask.ome.tif"
    for ii in filepaths.keys()
}


def run_segment_interstitial():
    def parallel_wrap(
        img_path,
        tubule_mask_path,
        tubule_ring_mask_path,
        tubule_type_mask_path,
        roi_path,
        out_path,
    ):
        interstitial_mask, _interstitial_mask = segment_interstitial(
            img_path, tubule_mask_path, tubule_ring_mask_path, tubule_type_mask_path
        )

        expands = [150, 300, 450]
        expands = np.divide(expands, 0.325)

        for ee in expands:
            expanded_glomeruli_mask = roi_path_to_mask(
                roi_path, mask_shape=interstitial_mask.shape, buffer_size=ee
            )
            expanded_glomeruli_mask[interstitial_mask == 0] = 0
            _ee = int(np.ceil(ee))
            palom.pyramid.write_pyramid(
                [da.from_array(expanded_glomeruli_mask, chunks=1024, name=False)],
                out_path.parent / out_path.name.replace(".ome.tif", f"_{_ee}.ome.tif"),
                pixel_size=0.325,
                downscale_factor=2,
                is_mask=True,
                compression="zlib",
                save_RAM=True,
            )

            expanded_glomeruli_mask = roi_path_to_mask(
                roi_path, mask_shape=_interstitial_mask.shape, buffer_size=ee
            )
            expanded_glomeruli_mask[_interstitial_mask == 0] = 0
            _ee = int(np.ceil(ee))
            palom.pyramid.write_pyramid(
                [da.from_array(expanded_glomeruli_mask, chunks=1024, name=False)],
                out_path.parent
                / out_path.name.replace(".ome.tif", f"_{_ee}_for_area.ome.tif"),
                pixel_size=0.325,
                downscale_factor=2,
                is_mask=True,
                compression="zlib",
                save_RAM=True,
            )

    joblib.Parallel(n_jobs=6, verbose=1)(
        joblib.delayed(parallel_wrap)(
            filepaths[nn],
            masks[nn],
            ring_masks[nn],
            type_masks[nn],
            rois[nn],
            out_paths[nn],
        )
        for nn in subset
    )
