import itertools
import pathlib

import anndata as ad
import matplotlib.cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import shap
import sklearn.feature_selection
import sklearn.linear_model
import sklearn.metrics
import sklearn.model_selection
import sklearn.pipeline
import sklearn.preprocessing
import tqdm
import tqdm.contrib
from palom.cli.align_he import set_matplotlib_font

from cmap import Colormap


# sns.set_palette("deep")
sns.set_palette(np.roll(Colormap("tol:muted").to_matplotlib().colors, 1, axis=0))
sns.set_palette(
    Colormap("tol:muted_alt").to_matplotlib().colors[[0, 1, 3, 4, 6, 5, 2, 7, 8]]
)

_colors = """
#cc6677
#88ccee
#ddcc77
#1f6bb8
#117733
#994455
#808080
#44aa99
#999933
#aa4499
""".strip().split("\n")
sns.set_palette(_colors)

sc.settings.verbosity = 3


def save_all_figs(dpi=300, format="pdf", out_dir=None, prefix=None):
    figs = [plt.figure(i) for i in plt.get_fignums()]
    if prefix is not None:
        for f in figs:
            if f._suptitle:
                f.suptitle(f"{prefix} {f._suptitle.get_text()}")
            else:
                f.suptitle(prefix)
    names = [f._suptitle.get_text() if f._suptitle else "" for f in figs]
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    for f, n, nm in zip(figs, plt.get_fignums(), names):
        f.savefig(out_dir / f"{n}-{nm}.{format}", dpi=dpi, bbox_inches="tight")


# ---------------------------------------------------------------------------- #
#                            Analyze organized data                            #
# ---------------------------------------------------------------------------- #
data_paths = """
/Users/yuanchen/HMS Dropbox/Yu-An Chen/000 local remote sharing/20250327-kidney_uni_feature/anndata-slide/ALIGNED/titan-HE_norm.h5ad
/Users/yuanchen/HMS Dropbox/Yu-An Chen/000 local remote sharing/20250327-kidney_uni_feature/anndata-slide/ALIGNED/intensity-JRL.h5ad
/Users/yuanchen/HMS Dropbox/Yu-An Chen/000 local remote sharing/20250327-kidney_uni_feature/anndata-slide/ALIGNED/uni-pc4_1_3_compartments.h5ad
""".strip().split("\n")


# /Users/yuanchen/HMS Dropbox/Yu-An Chen/000 local remote sharing/20250327-kidney_uni_feature/anndata-slide/ALIGNED/uni-pc4_1_3_preproc_norm_compartments.h5ad
# /Users/yuanchen/HMS Dropbox/Yu-An Chen/000 local remote sharing/20250327-kidney_uni_feature/anndata-slide/ALIGNED/intensity.h5ad
# /Users/yuanchen/HMS Dropbox/Yu-An Chen/000 local remote sharing/20250327-kidney_uni_feature/anndata-slide/ALIGNED/titan-pc4_1_3_preproc_norm.h5ad
# /Users/yuanchen/HMS Dropbox/Yu-An Chen/000 local remote sharing/20250327-kidney_uni_feature/anndata-slide/ALIGNED/conch-pc4_1_3-conch-compartments.h5ad
# /Users/yuanchen/HMS Dropbox/Yu-An Chen/000 local remote sharing/20250327-kidney_uni_feature/anndata-slide/ALIGNED/conch-pc4_1_3_preproc_norm-conch-compartments.h5ad
# /Users/yuanchen/HMS Dropbox/Yu-An Chen/000 local remote sharing/20250327-kidney_uni_feature/anndata-slide/ALIGNED/titan-pc4_1_3.h5ad


adatas = [ad.read_h5ad(pp) for pp in data_paths]

feature_names = [
    "TITAN / H&E (1)",
    "Orion-compartment (2)",
    "UNI-compartment / PC image (3)",
    "(1) + (2)",
    "(1) + (3)",
    # "(2) + (3)",
    # "(1) + (2) + (3)",
]
features = [adatas[0].X, adatas[1].X, adatas[2].X]
features.append(np.hstack([features[0], features[1]]))
features.append(np.hstack([features[0], features[2]]))
features.append(np.hstack([features[1], features[2]]))
features.append(np.hstack([features[0], features[1], features[2]]))


# ---------------------------------------------------------------------------- #
#                          aggregate clinical targets                          #
# ---------------------------------------------------------------------------- #
df_clinical = {}


# continuous
cols = """
Mean_Log2_KL,KL_Zscore,Mean_Log2_C3C1q,C3C1q_Zscore,PercentCorticalColIII,PercentGlomColIV,mean_W4,GlobalGlomerulosclerosis,GlobalAndSegmentalGlomerulosclerosis,InterstitialFibrosisTubularAtrophy,sCr_2,eGFR,slideName
""".strip().split(",")

cols = """
slideName,GlobalAndSegmentalGlomerulosclerosis,InterstitialFibrosisTubularAtrophy
""".strip().split(",")

rename_map = {
    "GlobalGlomerulosclerosis": "GSG",
    "GlobalAndSegmentalGlomerulosclerosis": r"% FSGS",
    "InterstitialFibrosisTubularAtrophy": r"% IFTA",
    "ArterialArteriolorSclerosisScore": "AAS",
    "FootEffacementScore": "Foot Effacement",
    "ChronicChangesScore": "Chronic\nChanges",
    #
    "Podocytopathy": "Podocy-\ntopathy",
    "Sec FSGS": "Sec\nFSGS",
    "Collapsing GP": "Collapsing\nGP",
    "Diabetic GN": "Diabetic\nGN",
    "Crescentic GN": "Crescentic\nGN",
    #
    "eGFR > 60": "eGFR\n> 60",
    "Foot Processes Binary": "GD-FPE",
}

df = pd.read_excel(
    r"/Users/yuanchen/HMS Dropbox/Yu-An Chen/000 local remote sharing/20250327-kidney_uni_feature/clinical-data/PhenotypeReadable.xlsx",
    index_col="slideName",
    usecols=cols,
).rename(columns=rename_map)
df = df.loc[adatas[0].obs["LSPID"]].copy().astype("float")
df_clinical["continuous"] = df


# binary
df = pd.read_excel(
    r"/Users/yuanchen/HMS Dropbox/Yu-An Chen/000 local remote sharing/20250327-kidney_uni_feature/clinical-data/MergedDiagnosisV2.xlsx",
    index_col="caseID",
)
df = df.loc[adatas[0].obs.index].copy().astype("bool")
df = df.join(adatas[0].obs[["eGFR > 60", "Foot Processes Binary"]].astype("bool"))

use_cols = [
    "Podocytopathy",
    "Sec FSGS",
    "Collapsing GP",
    "MN",
    "IgA GN",
    "IC-GN",
    "Diabetic GN",
    "Crescentic GN",
    "TBM",
    "AIN",
    "CAIN",
    "ATI",
    "Vasc Dis",
    # "binary features"
    "eGFR > 60",
    "Foot Processes Binary",
]
target_cols = df.loc[:, df.sum() >= 7][use_cols].columns.tolist()
df_clinical["binary"] = df[target_cols].rename(columns=rename_map).copy()


# ordinal
cols = """
ArterialArteriolorSclerosisScore,ChronicChangesScore,slideName
""".strip().split(",")

df = pd.read_excel(
    r"/Users/yuanchen/HMS Dropbox/Yu-An Chen/000 local remote sharing/20250327-kidney_uni_feature/clinical-data/PhenotypeReadable.xlsx",
    index_col="slideName",
    usecols=cols,
).rename(columns=rename_map)
# NA in some columns
df = df.loc[adatas[0].obs["LSPID"]].copy().astype("float")
df_clinical["ordinal"] = df


# ---------------------------------------------------------------------------- #
#                               feature selection                              #
# ---------------------------------------------------------------------------- #
feature_sets = (
    [adatas[0].X]
    + [adatas[1].X[:, i * 20 : (i + 1) * 20] for i in range(4)]
    + [adatas[2].X[:, i * 1024 : (i + 1) * 1024] for i in range(4)]
)
feature_set_names = [
    "TITAN",
    "Orion-GLO",
    "Orion-INT",
    "Orion-TUB",
    "Orion-VES",
    "UNI-GLO",
    "UNI-TUB",
    "UNI-INT",
    "UNI-VES",
]

df_feature_pvalue = pd.DataFrame()
df_feature_count = pd.DataFrame()


for ft, fn, dt in zip(
    ["binary", "continuous", "ordinal"],
    [
        sklearn.feature_selection.f_classif,
        sklearn.feature_selection.f_regression,
        sklearn.feature_selection.f_classif,
    ],
    ["bool", "float", "int"],
):
    curr = df_clinical[ft]

    for cname in curr:
        col = curr[cname]
        valid = pd.notna(col)

        for ff, nn in zip(feature_sets, feature_set_names):
            kbest = sklearn.feature_selection.SelectKBest(fn)
            kbest.fit(ff[valid], col[valid].astype(dt))
            df_feature_pvalue.loc[cname.replace("-\n", "").replace("\n", " "), nn] = (
                np.nanpercentile(-np.log10(kbest.pvalues_), 95)
            )
            df_feature_count.loc[cname.replace("-\n", "").replace("\n", " "), nn] = (
                np.sum(-np.log10(kbest.pvalues_) >= 2) / kbest.pvalues_.size
            )


set_matplotlib_font(10)
for cmap in [
    # Colormap("cubehelix:cubehelix").to_matplotlib(),
    # Colormap("cmasher:chroma").to_matplotlib(),
    "jet",
    # "rainbow",
    # "Spectral_r",
    # "turbo",
    # "coolwarm",
    # "icefire",
]:
    fig, _ = plt.subplots(figsize=(9, 3.6))
    sns.heatmap(df_feature_count.T * 100, annot=True, cmap=cmap, fmt=".1f")
    fig.suptitle("Percentage of significant features")
    fig.tight_layout()

fig, _ = plt.subplots(figsize=(9, 3.6))
sns.heatmap(df_feature_pvalue.T, annot=True, cmap="coolwarm", center=2)
fig.suptitle("95-percentile of feature -log10(p-values)")
fig.tight_layout()
# ---------------------------------------------------------------------------- #
#                               continuous target                              #
# ---------------------------------------------------------------------------- #
df = df_clinical["continuous"]
target_cols = df.columns.tolist()[:]


scorings = [
    "r2",
    "neg_mean_squared_error",
    "explained_variance",
]
scoring_funcs = {name: sklearn.metrics.get_scorer(name) for name in scorings}

_df_summary = pd.DataFrame(
    index=pd.MultiIndex.from_product([target_cols, range(8 * 4)]),
    columns=pd.MultiIndex.from_product(
        [[f"test_{ss}" for ss in scorings], feature_names]
    ),
)

random_seed = np.random.randint(1, 100)
random_seed = 57


for y_column_name, (feature, feature_name) in tqdm.tqdm(
    itertools.product(target_cols[:], zip(features, feature_names))
):
    mask = pd.notna(df[y_column_name])

    X = feature[mask]
    y_ori = df[y_column_name][mask].values

    bins = np.percentile(y_ori, np.linspace(0, 100, 6))
    bins[-1] += 1
    y = np.digitize(y_ori, bins)

    clf = sklearn.linear_model.RidgeCV(alphas=np.logspace(-3, 5, 10))

    clf = sklearn.pipeline.make_pipeline(
        sklearn.preprocessing.StandardScaler(),
        clf,
    )

    cv = sklearn.model_selection.RepeatedStratifiedKFold(
        n_repeats=8, n_splits=4, random_state=random_seed
    )

    cv_scores = sklearn.model_selection.cross_validate(
        clf, X, y_ori, scoring=scorings, cv=list(cv.split(X, y)), n_jobs=8
    )
    for ss in scorings:
        _df_summary.loc[y_column_name, (f"test_{ss}", feature_name)] = cv_scores[
            f"test_{ss}"
        ]


sname = "r2"
splot = sns.catplot(
    _df_summary.loc[:, f"test_{sname}"]
    .astype("float")
    .melt(ignore_index=False, var_name="Feature", value_name=sname)
    .reset_index(names=["Task", "cv"]),
    hue="Feature",
    y=sname,
    x="Task",
    kind="bar",
    errorbar=("ci", 95),
)
splot.ax.set_title(f"Feature Comparison ({clf.steps[1][0].capitalize()})")

df_summary = _df_summary.loc[:, "test_r2"].astype("float")
df_summary_cont = df_summary.copy()

# ---------------------------------------------------------------------------- #
#                                binary columns                                #
# ---------------------------------------------------------------------------- #
df = df_clinical["binary"]
target_cols = df.columns.tolist()[:]


scorings = [
    # suggested by chatgpt
    "precision",
    "recall",
    "f1",
    "balanced_accuracy",
    "matthews_corrcoef",
    # comparison
    "roc_auc",
]
scoring_funcs = {name: sklearn.metrics.get_scorer(name) for name in scorings}

_df_summary = pd.DataFrame(
    index=pd.MultiIndex.from_product([target_cols, range(8 * 4)]),
    columns=pd.MultiIndex.from_product(
        [[f"test_{ss}" for ss in scorings], feature_names]
    ),
)

random_seed = np.random.randint(1, 100)
random_seed = 57

for y_column_name, (feature, feature_name) in tqdm.tqdm(
    itertools.product(target_cols, zip(features, feature_names))
):
    mask = pd.notna(df[y_column_name])

    X = feature[mask]
    y_ori = df[y_column_name][mask].values.astype("bool")

    y = y_ori.copy()

    clf = sklearn.linear_model.LogisticRegressionCV(
        Cs=10,
        max_iter=2000,
        scoring="roc_auc",
        solver="lbfgs",
        random_state=random_seed,
        class_weight="balanced",
        n_jobs=8,
    )
    clf = sklearn.pipeline.make_pipeline(
        sklearn.preprocessing.StandardScaler(),
        clf,
    )

    cv = sklearn.model_selection.RepeatedStratifiedKFold(
        n_repeats=8, n_splits=4, random_state=random_seed
    )

    cv_scores = sklearn.model_selection.cross_validate(
        clf, X, y_ori, scoring=scorings, cv=list(cv.split(X, y)), n_jobs=8
    )
    for ss in scorings:
        _df_summary.loc[y_column_name, (f"test_{ss}", feature_name)] = cv_scores[
            f"test_{ss}"
        ]

df_summary = _df_summary.loc[:, "test_roc_auc"].astype("float")
df_summary_binary = df_summary.copy()


sname = "roc_auc"
splot = sns.catplot(
    _df_summary.loc[:, f"test_{sname}"]
    .astype("float")
    .melt(ignore_index=False, var_name="Feature", value_name=sname)
    .reset_index(names=["Task", "cv"]),
    hue="Feature",
    y=sname,
    x="Task",
    kind="bar",
    errorbar=("ci", 95),
)
splot.ax.set_title(f"Feature Comparison ({clf.steps[1][0].capitalize()})")


# ---------------------------------------------------------------------------- #
#                                ordinal columns                               #
# ---------------------------------------------------------------------------- #
df = df_clinical["ordinal"]
target_cols = df.columns.tolist()[:]


scorings = [
    # "r2",
    # "neg_mean_squared_error",
    # "explained_variance",
    # suggested by chatgpt
    # "precision",
    # "recall",
    # "f1",
    "balanced_accuracy",
    "matthews_corrcoef",
    # comparison
    "roc_auc_ovo",
    "roc_auc_ovr",
]
_df_summary = pd.DataFrame(
    index=pd.MultiIndex.from_product([target_cols, range(8 * 4)]),
    columns=pd.MultiIndex.from_product(
        [[f"test_{ss}" for ss in scorings], feature_names]
    ),
)


random_seed = np.random.randint(1, 100)
random_seed = 57
# random_seed = 25


models = {}
for y_column_name, (feature, feature_name) in tqdm.tqdm(
    itertools.product(target_cols[:], zip(features, feature_names))
):
    mask = pd.notna(df[y_column_name])

    X = feature[mask]
    y_ori = df[y_column_name][mask].values

    y = y_ori.astype("int")

    clf = sklearn.linear_model.LogisticRegressionCV(
        Cs=10,
        max_iter=5000,
        # scoring="roc_auc_ovo",
        solver="lbfgs",
        random_state=random_seed,
        # class_weight="balanced",
        n_jobs=8,
    )
    clf = sklearn.pipeline.make_pipeline(
        sklearn.preprocessing.StandardScaler(),
        clf,
    )

    cv = sklearn.model_selection.RepeatedStratifiedKFold(
        n_repeats=8, n_splits=4, random_state=random_seed
    )

    cv_scores = sklearn.model_selection.cross_validate(
        clf,
        X,
        y,
        scoring=scorings,
        cv=list(cv.split(X, y)),
        n_jobs=8,
        return_estimator=True,
    )

    for ss in scorings:
        _df_summary.loc[y_column_name, (f"test_{ss}", feature_name)] = cv_scores[
            f"test_{ss}"
        ]

    if y_column_name not in models:
        models[y_column_name] = {}
    models[y_column_name][feature_name] = cv_scores["estimator"]

df_summary = _df_summary.loc[:, "test_roc_auc_ovo"].astype("float")
df_summary_ord = df_summary.copy()

sname = "roc_auc_ovo"
sns.catplot(
    _df_summary.loc[:, f"test_{sname}"]
    .astype("float")
    .melt(ignore_index=False, var_name="Feature", value_name=sname)
    .reset_index(names=["Task", "cv"]),
    hue="Feature",
    y=sname,
    x="Task",
    kind="bar",
    errorbar=("ci", 95),
)
plt.gca().tick_params(axis="x", labelrotation=90)

# ------------------------------ output as table ----------------------------- #
out_scoring_excel = rf"/Users/yuanchen/HMS Dropbox/Yu-An Chen/000 local remote sharing/20250327-kidney_uni_feature/scoring-seed_{random_seed}-20251117.xlsx"

df_summary = pd.concat([df_summary_binary, df_summary_cont, df_summary_ord])

gb = df_summary.reset_index(0, names="Diagnosis").groupby("Diagnosis", sort=False)
df_table = gb.mean().map(lambda x: f"{x:.3f} ± ").astype("str") + gb.std().map(
    lambda x: f"{x:.3f}"
).astype("str")
df_table.index = df_table.index.str.replace("-\n", "").str.replace("\n", " ")
df_table["Scoring"] = (
    ["ROC AUC"] * len(df_summary_binary.index.get_level_values(0).unique())
    + ["R2"] * len(df_summary_cont.index.get_level_values(0).unique())
    + ["ROC AUC"] * len(df_summary_ord.index.get_level_values(0).unique())
)
df_table.reset_index().set_index(["Diagnosis", "Scoring"]).to_excel(out_scoring_excel)
# ---------------------------------------------------------------------------- #
#                                 SHAP explore                                 #
# ---------------------------------------------------------------------------- #
shap_values = {}
for y_column_name, (feature, feature_name) in tqdm.tqdm(
    itertools.product(target_cols[:], zip(features, feature_names))
):
    mask = pd.notna(df[y_column_name])

    X = feature[mask]
    y_ori = df[y_column_name][mask].values

    y = y_ori.astype("int")

    _shap_values = np.full((8, *X.shape, len(np.unique(y))), fill_value=np.nan)

    for fold, (train, test) in enumerate(cv.split(X, y)):
        nth_repeat = np.floor(fold / 4).astype("int")

        clf = models[y_column_name][feature_name][fold]
        model, norm = clf.steps[1][1], clf.steps[0][1]
        explainer = shap.LinearExplainer(model, norm.transform(X[train]))

        _explanation = explainer(norm.transform(X[test]))
        _shap_values[nth_repeat, test] = _explanation.values

    if y_column_name not in shap_values:
        shap_values[y_column_name] = {}
    shap_values[y_column_name][feature_name] = _shap_values


# ---------------------------------------------------------------------------- #
#                               post-hoc analysis                              #
# ---------------------------------------------------------------------------- #
import scikit_posthocs
from statannotations.Annotator import Annotator

value_name = "Scoring"

set_matplotlib_font(8)
splot = sns.catplot(
    df_summary.melt(
        ignore_index=False, var_name="Feature", value_name=value_name
    ).reset_index(names=["Task", "cv"]),
    hue="Feature",
    y=value_name,
    x="Task",
    kind="bar",
    err_kws={"linewidth": 1},
    capsize=0.2,
    # errorbar=("ci", 95),
)
splot.ax.axhline(0.5, linestyle="--", color="#666")
splot.ax.set_title(f"Feature Comparison ({clf.steps[1][0].capitalize()})")
sns.move_legend(
    splot,
    "upper right",
    bbox_to_anchor=(1, 1),
    ncol=5,
    title=None,
    frameon=False,
)
splot.figure.set_size_inches(13.4, 3.35)
splot.figure.tight_layout()


splot = sns.catplot(
    df_summary.melt(
        ignore_index=False, var_name="Feature", value_name=value_name
    ).reset_index(names=["Task", "cv"]),
    hue="Feature",
    y=value_name,
    x="Task",
    kind="box",
    # errorbar=("ci", 95),
)
# splot.ax.axhline(0.5, linestyle="--", color="#666")
splot.ax.set_title(f"Feature Comparison ({clf.steps[1][0].capitalize()})")
sns.move_legend(
    splot,
    "upper right",
    bbox_to_anchor=(1, 1),
    ncol=5,
    title=None,
    frameon=False,
)

target_cols = []
for cc in df_summary.index.get_level_values(0):
    if cc not in target_cols:
        target_cols.append(cc)
pvalues = []
pairs = []
for nn in target_cols[:]:
    data = (
        df_summary.loc[nn]
        .melt(ignore_index=False, var_name="Features", value_name=value_name)
        .reset_index(names="CV fold")
    )
    test_results = scikit_posthocs.posthoc_wilcoxon(
        data, group_col="Features", val_col=value_name, p_adjust="fdr_bh"
    )

    # ONLY the best model
    avg_rank = (
        data.groupby("CV fold")[value_name]
        .rank(pct=True)
        .groupby(data["Features"])
        .mean()
    )
    _best = avg_rank.nlargest(1).index[0]
    best = test_results.loc[_best, test_results.columns != _best]
    _pvalues = best.values
    _pairs = list(itertools.product([_best], best.index))
    _pairs = np.asarray(_pairs)[_pvalues < 2]

    # ALL combo
    # _pvalues = test_results.values[
    #     np.triu(np.ones_like(test_results, dtype="bool"), k=1)
    # ]
    # _pairs = np.asarray(list(itertools.combinations(test_results.index, 2)))[
    #     _pvalues <= 0.05
    # ]

    pvalues.extend(_pvalues[_pvalues < 2])
    pairs.extend([((nn, p1), (nn, p2)) for p1, p2 in _pairs])


annot = Annotator(
    ax=splot.ax,
    pairs=pairs,
    data=df_summary.melt(
        ignore_index=False, var_name="Feature", value_name=value_name
    ).reset_index(names=["Task", "cv"]),
    hue="Feature",
    y=value_name,
    x="Task",
)
_ = (
    annot.configure(test=None, test_short_name="zzz", line_width=1)
    .set_pvalues(pvalues)
    .annotate()
)

splot.ax.set_xlim(-0.5, len(target_cols) - 0.5)
splot.figure.set_size_inches(13.4, 3.35)
splot.figure.tight_layout()


# plot using rank
df_summary_rank = df_summary.melt(
    ignore_index=False, var_name="Feature", value_name=value_name
).reset_index(names=["Task", "cv"])
df_summary_rank.insert(
    0,
    column="Rank",
    value=df_summary_rank.groupby(["Task", "cv"]).rank(pct=True)[value_name],
)
splot = sns.catplot(
    df_summary_rank,
    hue="Feature",
    y="Rank",
    x="Task",
    kind="bar",
    err_kws={"linewidth": 1},
    capsize=0.2,
    # errorbar=("ci", 95),
)
sns.move_legend(
    splot,
    "upper right",
    bbox_to_anchor=(1, 1),
    ncol=5,
    title=None,
    frameon=False,
)

annot = Annotator(
    ax=splot.ax,
    pairs=pairs,
    data=df_summary_rank,
    hue="Feature",
    y="Rank",
    x="Task",
)
_ = (
    annot.configure(test=None, test_short_name="zzz", line_width=1)
    .set_pvalues(pvalues)
    .annotate()
)
splot.ax.yaxis.set_ticks(np.linspace(0, 1, 6))
splot.ax.set(
    xlim=(-0.5, len(target_cols) - 0.5),
    title=f"Feature Comparison (Rank, {clf.steps[1][0].capitalize()})",
)
splot.figure.set_size_inches(13.4, 3.35)
splot.figure.tight_layout()


# ---------------------------------------------------------------------------- #
#                          correlation of feature sets                         #
# ---------------------------------------------------------------------------- #
# def drop_multi_corr_features(feat, max_corr=0.8):
#     feat = np.asarray(feat)
#     corr_mx = np.abs(np.corrcoef(feat.T))
#     upper_tri = np.triu(corr_mx, k=1)
#     return feat[:, ~np.any(upper_tri > max_corr, axis=0)]

corrs = np.linspace(1, 0, 21)
c01 = np.abs(np.corrcoef(adatas[0].X.T, adatas[1].X.T)[:768, 768:])
v01 = np.array([np.all(c01 < cc, axis=1).sum() for cc in corrs])
c21 = np.abs(np.corrcoef(adatas[2].X.T, adatas[1].X.T)[:4096, 4096:])
v21 = np.array([np.all(c21 < cc, axis=1).sum() for cc in corrs])

c02 = np.abs(np.corrcoef(adatas[0].X.T, adatas[2].X.T)[:768, 768:])
v02 = np.array([np.all(c02 < cc, axis=1).sum() for cc in corrs])
c20 = np.abs(np.corrcoef(adatas[2].X.T, adatas[0].X.T)[:4096, 4096:])
v20 = np.array([np.all(c20 < cc, axis=1).sum() for cc in corrs])


fig, ax = plt.subplots()
sns.histplot(
    pd.concat(
        [
            pd.DataFrame(
                {
                    "Max(|Pearson's r w/ ORION|)": np.abs(
                        np.corrcoef(adatas[0].X.T, adatas[1].X.T)[:768, 768:]
                    ).max(axis=1),
                    "Feature": "TITAN",
                }
            ),
            pd.DataFrame(
                {
                    "Max(|Pearson's r w/ ORION|)": np.abs(
                        np.corrcoef(adatas[2].X.T, adatas[1].X.T)[:4096, 4096:]
                    ).max(axis=1),
                    "Feature": "UNI",
                }
            ),
        ]
    ),
    element="step",
    bins=np.linspace(0, 1, 81),
    hue="Feature",
    x="Max(|Pearson's r w/ ORION|)",
    stat="density",
    common_norm=False,
    cumulative=True,
    palette={"TITAN": "tab:blue", "UNI": "tab:green"},
)


# ----------------------------------- viz-1 ---------------------------------- #
corrs = np.linspace(1, 0, 41)
corr_results = np.zeros((3, 3, len(corrs)), dtype="float")
for rr, cc in list(itertools.product(range(3), range(3))):
    nn = adatas[rr].X.shape[1]
    corr = np.abs(np.corrcoef(adatas[rr].X.T, adatas[cc].X.T)[:nn, nn:])
    if rr == cc:
        corr = np.triu(corr, k=1)
    corr_results[rr, cc] = np.divide(
        [np.all(corr < mm, axis=1).sum() for mm in corrs], nn
    )


corr_plot = pd.concat(
    [
        pd.DataFrame(rr.T, columns=["TITAN", "ORION", "UNI"], index=corrs)
        .melt(
            ignore_index=False, var_name="Corr w/", value_name="Non-correlated fraction"
        )
        .assign(**{"Feature": f"{nn} ({cc})"})
        for rr, nn, cc in zip(corr_results, ["TITAN", "ORION", "UNI"], [768, 80, 4096])
    ]
).reset_index(names="|Pearson's r|")
splot = sns.relplot(
    corr_plot,
    x="|Pearson's r|",
    y="Non-correlated fraction",
    col="Feature",
    hue="Corr w/",
    kind="line",
    col_wrap=2,
)
for ax in splot.axes:
    ax.axvline(0.4, c="#aaa", linestyle="--")
    # ax.set_xlim(0, 1)
    # ax.set_ylim(0, 1)
ax.invert_xaxis()
sns.move_legend(splot, "upper left", bbox_to_anchor=(0.55, 0.45), frameon=False)
splot.figure.set_size_inches(5, 5)
plt.tight_layout()

splot = sns.relplot(
    corr_plot,
    x="|Pearson's r|",
    y="Non-correlated fraction",
    col="Corr w/",
    hue="Feature",
    kind="line",
    col_wrap=2,
)
for ax in splot.axes:
    ax.axvline(0.4, c="#aaa", linestyle="--")
ax.invert_xaxis()
sns.move_legend(splot, "upper left", bbox_to_anchor=(0.55, 0.45), frameon=False)
splot.figure.set_size_inches(5, 5)
plt.tight_layout()


# ----------------------------------- viz-2 ---------------------------------- #
corrs = np.linspace(0, 1, 41)
corr_results = np.zeros((3, 3, len(corrs)), dtype="float")
for rr, cc in list(itertools.product(range(3), range(3))):
    nn = adatas[rr].X.shape[1]
    corr = np.abs(np.corrcoef(adatas[rr].X.T, adatas[cc].X.T)[:nn, nn:])
    if rr == cc:
        corr = np.triu(corr, k=1)[:-1, 1:]
        nn -= 1
    corr_results[rr, cc] = np.divide(
        [np.any(corr >= mm, axis=1).sum() for mm in corrs], nn
    )


corr_plot = pd.concat(
    [
        pd.DataFrame(rr.T, columns=["TITAN", "ORION", "UNI"], index=corrs)
        .melt(ignore_index=False, var_name="Corr w/", value_name="Fraction")
        .assign(**{"Feature": f"{nn} ({cc})"})
        for rr, nn, cc in zip(corr_results, ["TITAN", "ORION", "UNI"], [768, 80, 4096])
    ]
).reset_index(names="|Pearson's r|")
splot = sns.relplot(
    corr_plot,
    x="|Pearson's r|",
    y="Fraction",
    col="Feature",
    hue="Corr w/",
    kind="line",
    col_wrap=2,
)
for ax in splot.axes:
    ax.axvline(0.4, c="#aaa", linestyle="--")
    # ax.set_xlim(0, 1)
    # ax.set_ylim(0, 1)
# ax.invert_xaxis()
sns.move_legend(splot, "upper left", bbox_to_anchor=(0.55, 0.45), frameon=False)
splot.figure.set_size_inches(5, 5)
plt.tight_layout()

splot = sns.relplot(
    corr_plot,
    x="|Pearson's r|",
    y="Fraction",
    col="Corr w/",
    hue="Feature",
    kind="line",
    col_wrap=2,
)
for ax in splot.axes:
    ax.axvline(0.4, c="#aaa", linestyle="--")
# ax.invert_xaxis()
sns.move_legend(splot, "upper left", bbox_to_anchor=(0.55, 0.45), frameon=False)
splot.figure.set_size_inches(5, 5)
plt.tight_layout()


splot = sns.catplot(
    corr_plot[corr_plot["|Pearson's r|"].astype("float32") == 0.4],
    x="Feature",
    y="Fraction",
    col="|Pearson's r|",
    hue="Corr w/",
    kind="bar",
    col_wrap=1,
)
splot.axes[0].set_title("|Pearson's r| ≥ 0.4")
splot.figure.set_size_inches(5, 3)
plt.tight_layout()
splot.figure.subplots_adjust(right=0.75)

# ----------------------------------- viz-3 ---------------------------------- #
# ---------------------------- consider all pairs ---------------------------- #
corrs = np.linspace(0, 1, 41)
corr_results = np.zeros((3, 3, len(corrs)), dtype="float")
for rr, cc in list(itertools.product(range(3), range(3))):
    nn = adatas[rr].X.shape[1]
    corr = np.abs(np.corrcoef(adatas[rr].X.T, adatas[cc].X.T)[:nn, nn:])
    corr_results[rr, cc] = np.divide([np.sum(corr >= mm) for mm in corrs], corr.size)


corr_plot = pd.concat(
    [
        pd.DataFrame(rr.T, columns=["TITAN", "ORION", "UNI"], index=corrs)
        .melt(ignore_index=False, var_name="Corr w/", value_name="Fraction")
        .assign(**{"Feature": f"{nn} ({cc})"})
        for rr, nn, cc in zip(corr_results, ["TITAN", "ORION", "UNI"], [768, 80, 4096])
    ]
).reset_index(names="|Pearson's r|")
splot = sns.relplot(
    corr_plot,
    x="|Pearson's r|",
    y="Fraction",
    col="Feature",
    hue="Corr w/",
    kind="line",
    col_wrap=2,
)
for ax in splot.axes:
    ax.axvline(0.4, c="#aaa", linestyle="--")
    # ax.set_xlim(0, 1)
    # ax.set_ylim(0, 1)
# ax.invert_xaxis()
sns.move_legend(splot, "upper left", bbox_to_anchor=(0.55, 0.45), frameon=False)
splot.figure.set_size_inches(5, 5)
plt.tight_layout()

splot = sns.relplot(
    corr_plot,
    x="|Pearson's r|",
    y="Fraction",
    col="Corr w/",
    hue="Feature",
    kind="line",
    col_wrap=2,
)
for ax in splot.axes:
    ax.axvline(0.4, c="#aaa", linestyle="--")
# ax.invert_xaxis()
sns.move_legend(splot, "upper left", bbox_to_anchor=(0.55, 0.45), frameon=False)
splot.figure.set_size_inches(5, 5)
plt.tight_layout()


splot = sns.catplot(
    corr_plot[corr_plot["|Pearson's r|"].astype("float32") == 0.4],
    x="Feature",
    y="Fraction",
    col="|Pearson's r|",
    hue="Corr w/",
    kind="bar",
    col_wrap=1,
)
splot.axes[0].set_title("|Pearson's r| ≥ 0.4")
splot.figure.set_size_inches(5, 3)
plt.tight_layout()
splot.figure.subplots_adjust(right=0.75)

# --------------------------- plot orion-titan corr -------------------------- #
set_matplotlib_font(12)


NAME = "TITAN"
DATA_ROW = adatas[1]
DATA_COL = adatas[0]
df_index = pd.DataFrame()
df_index[["Marker", "Compartment"]] = DATA_ROW.var_names.map(
    lambda x: x.split("_")[::-1]
).tolist()

df_corr = pd.DataFrame(
    np.corrcoef(DATA_COL.X[:, :], DATA_ROW.X, rowvar=False)[
        DATA_COL.shape[1] :, : DATA_COL.shape[1]
    ],
    index=pd.MultiIndex.from_frame(df_index),
)


# hierarchical clustering of TITAN features
compartment_order = ["GLO", "TUB", "INT", "VES"]
compartment_colors = ["#D81B60", "#004D40", "#1E88E5", "#FFC107"]
df_corr = df_corr.loc[pd.IndexSlice[:, compartment_order], :]

c_grid = sns.clustermap(
    df_corr,
    xticklabels=False,
    yticklabels=True,
    center=0,
    # mask=df_corr.abs() < 0.2,
    cmap="coolwarm",
    row_cluster=False,
    col_cluster=True,
    dendrogram_ratio=0.1,
    row_colors=df_corr.index.get_level_values("Compartment")
    .map(lambda x: compartment_colors[compartment_order.index(x)])
    .tolist(),
)
c_grid.ax_col_dendrogram.set_visible(False)
for aa in [c_grid.ax_row_colors, c_grid.ax_heatmap, c_grid.ax_cbar]:
    for ss in aa.spines.values():
        ss.set_visible(True)
c_grid.ax_cbar.fill_between(c_grid.ax_cbar.get_xlim(), -0.2, 0.2, color="w")

feature_order = c_grid.dendrogram_col.reordered_ind


hue_order = [
    "GroupCount",
    "DNA",
    "AF",
    "IgG",
    "IgA",
    "IgM",
    "Kappa",
    "Lambda",
    "C3",
    "C1q",
    "Fibrinogen",
    "Albumin",
    "PODXL",
    "Nephrin",
    "Synaptopodin",
    "KIM1",
    "AQP1",
    "ColIII",
    "ColIV",
    "SMA",
]
# per-compartment correlation

cmaps = ["jet", "rainbow", "Spectral_r", "turbo", "coolwarm", "icefire"]
cmaps = ["icefire"]
for cmap in cmaps:
    fig, axs = plt.subplots(
        len(compartment_order), 1, sharey=True, sharex=True, figsize=(12, 12)
    )
    for nn, ax in zip(compartment_order, axs):
        plot_df = df_corr.loc[pd.IndexSlice[hue_order, nn], :].iloc[:, feature_order]
        sns.heatmap(
            plot_df,
            ax=ax,
            # vmin=df_corr.min().min(),
            # vmax=df_corr.max().max(),
            vmin=-1,
            vmax=1,
            center=0,
            cmap=cmap,
            xticklabels=False,
            yticklabels=plot_df.index.get_level_values("Marker"),
            cbar=False,
            # mask=plot_df.abs() < 0.3,
            rasterized=True,
        )
        ax.tick_params(
            axis="y", labelleft=False, labelright=True, rotation=0, labelsize=8
        )
        ax.yaxis.tick_right()
        fraction = np.sum(plot_df.abs().max(axis=1) > 0.4) / plot_df.shape[0]
        ax.yaxis.set_label_text(
            f"{nn} ({(plot_df.abs() > 0.4).sum().sum():,} / {plot_df.size:,})"
        )
        # ax.yaxis.set_label_text(
        #     f"{nn} ({np.sum(plot_df.abs().max(axis=1) > 0.4)} / {plot_df.shape[0]})"
        # )

        for ss in ax.spines.values():
            ss.set_visible(True)
    cbar = fig.colorbar(ax.collections[0], ax=axs, fraction=0.02, location="bottom")
    fig.suptitle(f"Correlation of {NAME} Feature and ORION Intensity")
    # cbar.ax.fill_between(-0.3, 0.3, cbar.ax.get_xlim(), color="w")
    fig.subplots_adjust(
        top=0.95, bottom=0.2, left=0.125, right=0.9, hspace=0.015, wspace=0.0
    )

# titan-to-titan
c_grid = sns.clustermap(
    np.corrcoef(adatas[0].X.T),
    center=0,
    vmin=-1,
    vmax=1,
    rasterized=True,
    xticklabels=False,
    yticklabels=False,
    dendrogram_ratio=0.1,
    cbar_pos=(0.01, 0.3, 0.02, 0.4),
    # cbar_kws=dict(location="bottom")
)
c_grid.ax_col_dendrogram.set_visible(False)
c_grid.ax_row_dendrogram.set_visible(False)
c_grid.figure.suptitle("Correlation of TITAN Features")


# ---------------------------------------------------------------------------- #
#                                 ROC-AUC curve                                #
# ---------------------------------------------------------------------------- #
df = df_clinical["binary"]
target_cols = df.columns.tolist()[:]


random_seed = 57
# ---------------------------- one task per figure --------------------------- #
figs = [
    plt.subplots(1, 5, sharex=True, sharey=True, figsize=(19, 4.5))
    for _ in target_cols[:]
]

for (fig, axs), y_column_name in tqdm.contrib.tzip(figs, target_cols[:]):
    fig.suptitle(y_column_name.replace("-\n", "").replace("\n", " "))

    for ii, (feature, feature_name) in enumerate(zip(features, feature_names)):
        mask = pd.notna(df[y_column_name])

        X = feature[mask]
        y_ori = df[y_column_name][mask].values.astype("bool")

        y = y_ori.copy()

        clf = sklearn.linear_model.LogisticRegressionCV(
            Cs=10,
            max_iter=2000,
            scoring="roc_auc",
            solver="lbfgs",
            random_state=random_seed,
            class_weight="balanced",
            n_jobs=8,
        )
        clf = sklearn.pipeline.make_pipeline(
            sklearn.preprocessing.StandardScaler(),
            clf,
        )

        cv = sklearn.model_selection.RepeatedStratifiedKFold(
            n_repeats=8, n_splits=4, random_state=random_seed
        )

        tprs = []
        aucs = []
        mean_fpr = np.linspace(0, 1, 100)

        ax = axs[ii]
        for fold, (train, test) in enumerate(cv.split(X, y)):
            clf.fit(X[train], y[train])
            args = y[test], clf.predict_proba(X[test])[:, 1]
            aucs.append(sklearn.metrics.roc_auc_score(*args))
            fpr, tpr, _ = sklearn.metrics.roc_curve(*args)
            interp_tpr = np.interp(mean_fpr, fpr, tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)

        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = sklearn.metrics.auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)

        ax.plot(
            [0, 1],
            [0, 1],
            color="k",
            linestyle="--",
            label="Chance level (AUC = 0.5)",
        )

        ax.plot(
            mean_fpr,
            mean_tpr,
            color="b",
            label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
            lw=2,
            alpha=0.8,
        )

        std_tpr = np.std(tprs, axis=0)
        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        ax.fill_between(
            mean_fpr,
            tprs_lower,
            tprs_upper,
            color="grey",
            alpha=0.2,
            label=r"$\pm$ 1 std. dev.",
        )

        ax.set(
            title=feature_name,
            aspect="equal",
            xlim=(-0.02, 1.02),
            ylim=(-0.02, 1.02),
            xlabel="False Positive Rate",
            ylabel="True Positive Rate" if ii == 0 else None,
        )
        ax.legend(loc="lower right")
    fig.tight_layout()


# ---------------------------- one panel per task ---------------------------- #
fig, axs = plt.subplots(3, 5, sharex=True, sharey=True)

# fig, axs = plt.subplots(1, 5, sharex=True, sharey=True)
# bcols = ["Foot Processes Binary", "eGFR > 60", "Sec FSGS", "Crescentic GN", "IgA GN"]
for ax, y_column_name in tqdm.contrib.tzip(axs.ravel(), target_cols[:]):
    ax.plot(
        [0, 1],
        [0, 1],
        color="k",
        linestyle="--",
        label="Chance level (AUC = 0.5)",
    )
    for ii, (feature, feature_name) in enumerate(zip(features, feature_names)):
        mask = pd.notna(df[y_column_name])

        X = feature[mask]
        y_ori = df[y_column_name][mask].values.astype("bool")

        y = y_ori.copy()

        clf = sklearn.linear_model.LogisticRegressionCV(
            Cs=10,
            max_iter=2000,
            scoring="roc_auc",
            solver="lbfgs",
            random_state=random_seed,
            class_weight="balanced",
            n_jobs=8,
        )
        clf = sklearn.pipeline.make_pipeline(
            sklearn.preprocessing.StandardScaler(),
            clf,
        )

        cv = sklearn.model_selection.RepeatedStratifiedKFold(
            n_repeats=8, n_splits=4, random_state=random_seed
        )

        tprs = []
        aucs = []
        mean_fpr = np.linspace(0, 1, 100)

        for fold, (train, test) in enumerate(cv.split(X, y)):
            clf.fit(X[train], y[train])
            args = y[test], clf.predict_proba(X[test])[:, 1]
            aucs.append(sklearn.metrics.roc_auc_score(*args))
            fpr, tpr, _ = sklearn.metrics.roc_curve(*args)
            interp_tpr = np.interp(mean_fpr, fpr, tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)

        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = sklearn.metrics.auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)

        ax.plot(
            mean_fpr,
            mean_tpr,
            label=f"{feature_name}\n(AUC = {mean_auc:0.2f} ± {std_auc:0.2f})",
            lw=2,
            alpha=0.8,
        )

    ax.set(
        title=y_column_name.replace("-\n", "").replace("\n", " "),
        aspect="equal",
        xlim=(-0.02, 1.02),
        ylim=(-0.02, 1.02),
        xlabel="False Positive Rate",
        ylabel="True Positive Rate",
    )
    ax.legend(loc="lower right")


# ---------------------------------------------------------------------------- #
#                          record number of valid data                         #
# ---------------------------------------------------------------------------- #
df_total = pd.concat([kk.notna().sum() for kk in df_clinical.values()])
df_total.index = df_total.index.str.replace("-\n", "").str.replace("\n", " ")

df_positive = (df_clinical["binary"] > 0).sum()
df_positive.index = df_positive.index.str.replace("-\n", "").str.replace("\n", " ")

