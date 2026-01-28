import os, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.neighbors import KNeighborsRegressor
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
from matplotlib.path import Path as MplPath


# =========================================================
# 路径（按你给的）
# =========================================================
BASE_DIR = r"Y:\long\publication_datasets\liver\output-XETG00082_C105"

CLUSTERS_CSV = r"Y:\long\publication_datasets\liver\output-XETG00082_C105\analysis\clustering\gene_expression_graphclust\clusters.csv"
CELLS_PARQUET = r"Y:\long\publication_datasets\liver\output-XETG00082_C105\cells.parquet"
TISSUE_BOUNDARY_CSV = r"Y:\long\publication_datasets\liver\output-XETG00082_C105\tissue_boundary.csv"

OUT_SUBDIR = "pattern1_isoline0p5_from_graphclust"

# =========================================================
# 你的 pattern1 clusters
# =========================================================
PATTERN1 = [10, 23, 19, 27, 14, 20, 25, 26]


# =========================================================
# 参数（先用默认，必要时再调）
# =========================================================
grid_n = 1200
knn_k = 30
smooth_sigma = 5

margin_um = 50.0
max_dist_threshold = 200

bg_d_min = 20
bg_d_max = 250
bg_max_points = 60000
random_state = 0

min_cells_inside = 10  # 防止碎片 contour

# synthetic bg：建议开
USE_SYN_BG = True
BBOX_EXPAND_UM = 100.0
SYN_BG_DENSITY = 0.01
SYN_BG_MIN = 20000
SYN_BG_MAX = 120000


# =========================================================
# utils
# =========================================================
def make_mesh_from_xy(xy, grid_n=800, pad=0.02, margin_um=50.0):
    xy = np.asarray(xy, float)
    xmin, ymin = xy.min(axis=0)
    xmax, ymax = xy.max(axis=0)
    dx, dy = xmax - xmin, ymax - ymin

    xmin -= dx * pad; xmax += dx * pad
    ymin -= dy * pad; ymax += dy * pad

    xmin -= margin_um; xmax += margin_um
    ymin -= margin_um; ymax += margin_um

    xs = np.linspace(xmin, xmax, grid_n)
    ys = np.linspace(ymin, ymax, grid_n)
    xx, yy = np.meshgrid(xs, ys)
    grid = np.c_[xx.ravel(), yy.ravel()]
    return xx, yy, grid

def tissue_mask_from_xy(all_xy, xx, yy, max_dist_threshold=200):
    grid = np.c_[xx.ravel(), yy.ravel()]
    tree = cKDTree(np.asarray(all_xy, float))
    dist, _ = tree.query(grid, k=1)
    return (dist.reshape(xx.shape) <= max_dist_threshold)

def extract_contour_paths(xx, yy, z2d, level=0.5):
    fig, ax = plt.subplots()
    cs = ax.contour(xx, yy, z2d, levels=[level])
    plt.close(fig)

    verts_list = []
    if hasattr(cs, "allsegs") and cs.allsegs and len(cs.allsegs) > 0:
        for seg in cs.allsegs[0]:
            v = np.asarray(seg)
            if v.ndim == 2 and v.shape[0] >= 10 and v.shape[1] == 2:
                verts_list.append(v)
        return verts_list

    if hasattr(cs, "collections") and cs.collections:
        for p in cs.collections[0].get_paths():
            v = p.vertices
            if len(v) >= 10:
                verts_list.append(v)
        return verts_list

    return []

def filter_loops_by_cell_count(verts_list, cells_xy, min_cells_inside=1):
    kept = []
    for v in verts_list:
        path = MplPath(v)
        if int(path.contains_points(cells_xy).sum()) >= min_cells_inside:
            kept.append(v)
    return kept

def load_tissue_boundary_csv(boundary_csv):
    df = pd.read_csv(boundary_csv)
    if {"x", "y"}.issubset(df.columns):
        return df[["x", "y"]].to_numpy(float)
    if {"X", "Y"}.issubset(df.columns):
        return df[["X", "Y"]].to_numpy(float)
    raise ValueError(f"tissue_boundary.csv 必须包含 x,y 或 X,Y 列，当前列={list(df.columns)}")

def generate_synthetic_bg_in_bbox(boundary_xy, expand_um=100.0, density=0.01,
                                  min_n=20000, max_n=120000, seed=0):
    rng = np.random.default_rng(seed)
    xmin, ymin = boundary_xy.min(axis=0)
    xmax, ymax = boundary_xy.max(axis=0)
    xmin -= expand_um; xmax += expand_um
    ymin -= expand_um; ymax += expand_um

    area = (xmax - xmin) * (ymax - ymin)
    n = int(area * density)
    n = max(min_n, min(max_n, n))

    xs = rng.uniform(xmin, xmax, size=n)
    ys = rng.uniform(ymin, ymax, size=n)
    return np.c_[xs, ys].astype(float)

def sample_background_from_other_cells_plus_synth(
    cells_df,
    synthetic_bg_xy,
    target_ids,
    target_xy,
    cell_id_col,
    x_col, y_col,
    d_min=20, d_max=250,
    max_points=60000, seed=0,
    margin_um=50.0
):
    rng = np.random.default_rng(seed)

    cid_all = cells_df[cell_id_col].astype(str).to_numpy()
    target_ids_arr = np.array(list(target_ids), dtype=str)
    is_bg = ~np.isin(cid_all, target_ids_arr)

    bg_real_df = cells_df.loc[is_bg, [x_col, y_col]].copy()
    bg_real = bg_real_df[[x_col, y_col]].to_numpy(float) if len(bg_real_df) > 0 else np.empty((0, 2), float)

    bg_syn = np.asarray(synthetic_bg_xy, float) if synthetic_bg_xy is not None else np.empty((0, 2), float)

    bg_xy = bg_real
    if len(bg_syn) > 0:
        bg_xy = bg_syn if len(bg_xy) == 0 else np.vstack([bg_xy, bg_syn])
    if len(bg_xy) == 0:
        return np.empty((0, 2), float)

    xmin, ymin = target_xy.min(axis=0)
    xmax, ymax = target_xy.max(axis=0)
    pad = d_max + margin_um
    in_box = (
        (bg_xy[:, 0] >= xmin - pad) & (bg_xy[:, 0] <= xmax + pad) &
        (bg_xy[:, 1] >= ymin - pad) & (bg_xy[:, 1] <= ymax + pad)
    )
    bg_xy = bg_xy[in_box]
    if len(bg_xy) == 0:
        return np.empty((0, 2), float)

    tree = cKDTree(np.asarray(target_xy, float))
    dist, _ = tree.query(bg_xy, k=1)
    keep = (dist >= d_min) & (dist <= d_max)
    bg_xy = bg_xy[keep]
    if len(bg_xy) == 0:
        return np.empty((0, 2), float)

    if len(bg_xy) > max_points:
        idx = rng.choice(len(bg_xy), size=max_points, replace=False)
        bg_xy = bg_xy[idx]

    return bg_xy


# =========================================================
# 关键：对齐 clusters.csv 的 Barcode 和 cells.parquet 的 id
# =========================================================
def align_clusters_with_cells(clusters_csv, cells_parquet):
    cl = pd.read_csv(clusters_csv)
    # 兼容列名
    if "Barcode" not in cl.columns or "Cluster" not in cl.columns:
        raise ValueError(f"clusters.csv 需要包含 Barcode/Cluster 列，当前列={list(cl.columns)}")

    cl["Barcode"] = cl["Barcode"].astype(str)
    cl["Cluster"] = pd.to_numeric(cl["Cluster"], errors="coerce").astype("Int64")

    cells = pd.read_parquet(cells_parquet)

    # 找坐标列（尽量自动识别）
    cand_x = [c for c in cells.columns if c.lower() in ["x", "x_centroid", "x_center", "xcoord", "x_coord"]]
    cand_y = [c for c in cells.columns if c.lower() in ["y", "y_centroid", "y_center", "ycoord", "y_coord"]]
    if not cand_x or not cand_y:
        # Xenium/其他可能是 x_centroid/y_centroid
        # 如果你这里失败，把 cells.columns 打印出来即可
        raise ValueError(f"cells.parquet 找不到 x/y 列。列名示例：{list(cells.columns)[:60]}")

    x_col = cand_x[0]
    y_col = cand_y[0]

    # 找可能的 barcode/id 列
    id_candidates = []
    for c in cells.columns:
        lc = c.lower()
        if lc in ["barcode", "barcodes", "cell_barcode", "cellbarcode", "spot_barcode", "spot_id",
                  "cell_id", "cellid", "id"]:
            id_candidates.append(c)

    # 如果没找到典型列名，就退而求其次：所有 object/string 列都试
    if not id_candidates:
        id_candidates = [c for c in cells.columns if cells[c].dtype == object][:10]

    # 逐列尝试 merge（优先原样匹配，再尝试去掉 -1 匹配）
    def try_merge(cells_id_col, strip_suffix=False):
        tmp = cells.copy()
        tmp["_join_id"] = tmp[cells_id_col].astype(str)
        cl2 = cl.copy()
        cl2["_join_id"] = cl2["Barcode"].astype(str)

        if strip_suffix:
            tmp["_join_id"] = tmp["_join_id"].str.replace(r"-1$", "", regex=True)
            cl2["_join_id"] = cl2["_join_id"].str.replace(r"-1$", "", regex=True)

        m = tmp.merge(cl2[["_join_id", "Cluster"]], on="_join_id", how="inner")
        return m

    best = None
    best_info = None
    for c in id_candidates:
        m1 = try_merge(c, strip_suffix=False)
        if best is None or len(m1) > len(best):
            best, best_info = m1, (c, False)
        m2 = try_merge(c, strip_suffix=True)
        if best is None or len(m2) > len(best):
            best, best_info = m2, (c, True)

    if best is None or len(best) == 0:
        # 打印帮助信息
        print("[FAIL] 无法将 clusters.csv 的 Barcode 对齐到 cells.parquet")
        print("clusters.csv Barcode 示例:", cl["Barcode"].head().tolist())
        print("cells.parquet 列名:", list(cells.columns)[:80])
        for c in id_candidates[:6]:
            print(f"cells[{c}] 示例:", cells[c].astype(str).head().tolist())
        raise RuntimeError("ID 对齐失败：cells.parquet 里可能没有 barcode，需要一个 barcode↔cell_id 映射文件。")

    id_col_used, stripped = best_info
    print(f"[OK] merge success. cells id col='{id_col_used}', strip_suffix(-1)={stripped}, matched={len(best)}")

    # 返回：带坐标 + cluster 的表
    out = best.rename(columns={"Cluster": "cluster"})
    return out, id_col_used, x_col, y_col


# =========================================================
# main：pattern1=1 others=0，取 0.5 isoline
# =========================================================
def run_pattern1_isoline():
    os.makedirs(os.path.join(BASE_DIR, OUT_SUBDIR), exist_ok=True)
    out_dir = os.path.join(BASE_DIR, OUT_SUBDIR)

    merged, id_col_used, x_col, y_col = align_clusters_with_cells(CLUSTERS_CSV, CELLS_PARQUET)

    merged["cluster"] = pd.to_numeric(merged["cluster"], errors="coerce").astype("Int64")
    merged = merged.dropna(subset=["cluster"]).copy()
    merged["cluster"] = merged["cluster"].astype(int)

    p1 = set(int(x) for x in PATTERN1)
    merged["_is_p1"] = merged["cluster"].isin(p1)

    # 目标细胞（pattern1）
    p1_df = merged.loc[merged["_is_p1"], [id_col_used, x_col, y_col]].copy()
    if len(p1_df) < 10:
        raise RuntimeError(f"pattern1 cells too few after merge: {len(p1_df)}")

    target_ids = set(p1_df[id_col_used].astype(str))
    target_xy = p1_df[[x_col, y_col]].to_numpy(float)

    # 组织边界 + synthetic bg
    syn_bg_xy = None
    if USE_SYN_BG:
        boundary_xy = load_tissue_boundary_csv(TISSUE_BOUNDARY_CSV)
        syn_bg_xy = generate_synthetic_bg_in_bbox(
            boundary_xy,
            expand_um=BBOX_EXPAND_UM,
            density=SYN_BG_DENSITY,
            min_n=SYN_BG_MIN,
            max_n=SYN_BG_MAX,
            seed=random_state
        )
        print(f"[INFO] synthetic bg points: {len(syn_bg_xy)}")

    # 背景点：用 merged 里“非 pattern1 的细胞”作为真实背景来源
    # 注意：这里用 merged 作为 cells_df（因为它一定有 id_col_used 和 x/y）
    bg0_xy = sample_background_from_other_cells_plus_synth(
        cells_df=merged.rename(columns={id_col_used: "tmp_id"}),
        synthetic_bg_xy=syn_bg_xy,
        target_ids=set([str(x) for x in target_ids]),
        target_xy=target_xy,
        cell_id_col="tmp_id",
        x_col=x_col, y_col=y_col,
        d_min=bg_d_min, d_max=bg_d_max,
        max_points=bg_max_points,
        seed=random_state,
        margin_um=margin_um
    )
    if len(bg0_xy) == 0:
        raise RuntimeError("No bg0 points sampled. Try relaxing bg_d_min/bg_d_max.")

    # 训练：pattern1=1，bg0=0
    X_train = np.vstack([bg0_xy, target_xy])
    y_train = np.hstack([np.zeros(len(bg0_xy)), np.ones(len(target_xy))])

    reg = KNeighborsRegressor(n_neighbors=knn_k, weights="distance")
    reg.fit(X_train, y_train)

    # grid + predict + smooth
    xx, yy, grid = make_mesh_from_xy(target_xy, grid_n=grid_n, margin_um=margin_um)
    prob = reg.predict(grid).reshape(xx.shape)
    prob_smooth = gaussian_filter(prob, sigma=smooth_sigma)

    # tissue mask：用 merged 的所有真实细胞坐标（不含虚拟点）
    all_xy = merged[[x_col, y_col]].to_numpy(float)
    tissue_mask = tissue_mask_from_xy(all_xy, xx, yy, max_dist_threshold=max_dist_threshold)

    prob_smooth_masked = prob_smooth.copy()
    prob_smooth_masked[~tissue_mask] = np.nan

    # ✅ 0.5 isoline
    verts_list = extract_contour_paths(xx, yy, prob_smooth_masked, level=0.5)
    verts_list = filter_loops_by_cell_count(verts_list, target_xy, min_cells_inside=min_cells_inside)

    if len(verts_list) == 0:
        raise RuntimeError(
            "No 0.5 isoline found.\n"
            "建议：min_cells_inside 降低（50->10），smooth_sigma 增大（5->8），knn_k 增大（30->50）。"
        )

    # 保存参数
    params = dict(
        clusters_csv=CLUSTERS_CSV,
        cells_parquet=CELLS_PARQUET,
        tissue_boundary_csv=TISSUE_BOUNDARY_CSV,
        id_col_used=id_col_used,
        x_col=x_col, y_col=y_col,
        pattern1_clusters=sorted(list(p1)),
        grid_n=grid_n, knn_k=knn_k, smooth_sigma=smooth_sigma,
        bg_d_min=bg_d_min, bg_d_max=bg_d_max, bg_max_points=bg_max_points,
        max_dist_threshold=max_dist_threshold,
        isoline_level=0.5,
        min_cells_inside=min_cells_inside,
        use_synth_bg=USE_SYN_BG,
        n_target_cells=int(len(target_xy)),
        n_bg0=int(len(bg0_xy)),
        n_contours=int(len(verts_list)),
    )
    with open(os.path.join(out_dir, "params.json"), "w", encoding="utf-8") as f:
        json.dump(params, f, indent=2, ensure_ascii=False)

    # 保存 contour
    for i, v in enumerate(verts_list):
        np.save(os.path.join(out_dir, f"pattern1_isoline0p5_{i}.npy"), v)

    # 画图检查
    plt.figure(figsize=(10, 10))
    plt.scatter(bg0_xy[:, 0], bg0_xy[:, 1], s=1, alpha=0.05, label="bg0 (other cells + synth)")
    plt.scatter(target_xy[:, 0], target_xy[:, 1], s=3, alpha=0.85, label="pattern1 cells")
    for v in verts_list:
        plt.plot(v[:, 0], v[:, 1], linewidth=2)
    plt.gca().set_aspect("equal")
    plt.title(f"Pattern1 segmentation | isoline=0.5 | contours={len(verts_list)}")
    plt.legend(frameon=False)
    plt.tight_layout()

    fig_path = os.path.join(out_dir, "pattern1_isoline0p5.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()

    print("[OK] saved:", out_dir)
    print(" -", fig_path)
    print(" - contours npy:", len(verts_list))


if __name__ == "__main__":
    run_pattern1_isoline()
