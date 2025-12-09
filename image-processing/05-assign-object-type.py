import pathlib
import re

import cv2
import joblib
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import scipy.spatial.distance as sdistance
import shapely
import shapely.validation
import tifffile
import zarr


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


def roi_path_to_mask(roi_path, mask_shape=None, fill_value_min=None, buffer_size=0):
    _roi = pd.read_csv(roi_path)
    is_polygon = _roi["type"].fillna("").astype("str").apply(lambda x: x == "Polyline")
    assert len(_roi) <= 255
    roi = _roi.loc[is_polygon].apply(add_mpatch, axis=1)

    polygons = [
        to_valid_polygon(patch_to_shapely_polygon(pp).buffer(buffer_size))
        for pp in roi["mpatch"].dropna()
    ]
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


def _assign_type(quant_path, type_dict):
    quant = pd.read_csv(quant_path, index_col="CellID")
    if "interstitial_mask" in str(quant_path):
        quant["Type"] = "IN"
    else:
        quant["Type"] = "TU"
        quant.loc[type_dict["GL"], "Type"] = "GL"
        quant.loc[type_dict["GS"], "Type"] = "GS"
    quant.to_csv(quant_path)


def run_one(
    roi_path,
    tubule_mask_path,
    quant_path_tubule,
    quant_path_ring,
    quant_path_interstitial,
):
    mask = roi_path_to_mask(roi_path, buffer_size=-20)

    _roi = pd.read_csv(roi_path)
    is_polygon = _roi["type"].fillna("").astype("str").apply(lambda x: x == "Polyline")
    assert len(_roi) <= 255
    roi = _roi.loc[is_polygon]

    assert mask.max() == len(roi)

    tubule_mask = zarr.open(tifffile.imread(tubule_mask_path, aszarr=True, level=0))
    mask = roi_path_to_mask(roi_path, buffer_size=-20, mask_shape=tubule_mask.shape)

    indexer = np.zeros(mask.max() + 1, dtype="uint8")
    for ii, tt in enumerate(roi["Text"].fillna("GL")):
        tt = tt.upper()
        assert tt in ["GL", "GS"]
        if tt == "GL":
            indexer[ii + 1] = 1
        else:
            indexer[ii + 1] = 2

    mask = indexer[mask]

    object_type = {
        "GL": np.unique(tubule_mask.vindex[np.where(mask == 1)]),
        "GS": np.unique(tubule_mask.vindex[np.where(mask == 2)]),
    }
    assert len(object_type["GL"]) + len(object_type["GS"]) == len(
        roi
    ), f"\n{roi_path}\n"

    for qq in [quant_path_tubule, quant_path_ring, quant_path_interstitial]:
        _assign_type(qq, object_type)
    return


rois = sorted(
    pathlib.Path(
        r"X:\cycif-production\152-Kidney_Imaging_ProspectiveCasesRound2-2\ROI export"
    ).glob("*-rois.csv")
)
filepaths = {ff.name.split("-")[0]: ff for ff in rois}

out_mask_dir = pathlib.Path(
    r"Z:\yc00000\computation\YC-20241105-kidney-segmentation\tubule-mask"
)
masks = {
    ii: next(out_mask_dir.glob(f"{ii}-tubule_mask.ome.tif")) for ii in filepaths.keys()
}

quant_dir = pathlib.Path(
    r"Z:\yc00000\computation\YC-20241105-kidney-segmentation\quantification"
)

quant_tubule_paths = {
    ii: next(quant_dir.glob(f"{ii}-*-tubule_mask.csv")) for ii in filepaths.keys()
}
quant_ring_paths = {
    ii: next(quant_dir.glob(f"{ii}-*-tubule_ring_mask_30.csv"))
    for ii in filepaths.keys()
}
quant_interstitial_paths = {
    ii: next(quant_dir.glob(f"{ii}-*-interstitial_mask.csv")) for ii in filepaths.keys()
}


slides = list(filepaths.keys())
subset = slides[:]

_ = joblib.Parallel(n_jobs=8, verbose=1)(
    joblib.delayed(run_one)(
        filepaths[nn],
        masks[nn],
        quant_tubule_paths[nn],
        quant_ring_paths[nn],
        quant_interstitial_paths[nn],
    )
    for nn in subset
)
