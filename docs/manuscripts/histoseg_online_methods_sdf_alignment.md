# HistoSeg Online Methods: SDF Quantification And Topology-Aware Alignment

This draft describes the current HistoSeg implementation for 3D structure
quantification and multi-slice contour alignment. It is intentionally limited
to algorithms present in the codebase.

## Anisotropic Signed-Distance Fields For 3D Structure Masks

HistoSeg represents each annotated tissue structure as a binary voxel mask on a
3D grid ordered as `(z, y, x)`. The grid is defined by the gene-density
voxelization metadata, with physical spacing

$$
s = (\Delta z, \Delta y, \Delta x) = (z_{\mu m}, y_{\mu m}, x_{\mu m}).
$$

Aligned per-slice GeoJSON contours are rasterized onto this grid. For each
structure \(s_j\), HistoSeg constructs a mask

$$
M_j(v) \in \{0,1\},
$$

where \(v\) indexes a voxel in `(z, y, x)` order and \(M_j(v)=1\) indicates
that the voxel center lies inside the aligned contour for structure \(s_j\) on
the corresponding slice. The mask spacing is stored as
`spacing_zyx_um = (z_um, y_um, x_um)`.

Signed-distance values are computed from the binary mask using
`scipy.ndimage.distance_transform_edt` with physical sampling. HistoSeg first
computes the outside and inside Euclidean distance transforms:

$$
D_{\mathrm{out},j}(v) =
\operatorname{EDT}\left(1 - M_j; s\right),
$$

$$
D_{\mathrm{in},j}(v) =
\operatorname{EDT}\left(M_j; s\right).
$$

The signed-distance field is then assembled by applying the implementation's
inside/outside split:

$$
D_j(v) =
\begin{cases}
D_{\mathrm{out},j}(v), & M_j(v)=0,\\
-D_{\mathrm{in},j}(v), & M_j(v)=1.
\end{cases}
$$

Distances are therefore reported in microns and respect the anisotropic voxel
spacing through the `sampling=(z_um, y_um, x_um)` argument. Negative values mean
that a voxel is embedded inside the structure mask, while positive values mean
that a voxel is outside the structure. Boundary voxels are not clamped to zero:
because the implementation uses the discrete inside and outside transforms
directly, an isolated inside voxel can have signed distance \(-\min(s)\).
In the output metadata, this contract is recorded as `outside_um` outside the
structure and `-inside_um` inside the structure.

## SDF-Based Gene Spatial Quantification

For each requested gene, HistoSeg maps aligned cells to the same `(z, y, x)`
voxel grid. Cells outside the grid are excluded from the density calculation.
The raw total cell count and gene-expression sum are accumulated per voxel and
smoothed with a Gaussian kernel whose standard deviations are represented in
voxel units:

$$
\sigma = \left(
\frac{\sigma_z}{\Delta z},
\frac{\sigma_{xy}}{\Delta y},
\frac{\sigma_{xy}}{\Delta x}
\right).
$$

The enrichment field is the smoothed expression sum divided by the smoothed
total cell count:

$$
E_g(v) =
\frac{\operatorname{smooth}(X_g)(v)}
{\operatorname{smooth}(N)(v)}.
$$

Division is evaluated only where the smoothed total cell count is positive.
Voxels are considered valid when the smoothed total cell count exceeds the
configured minimum valid cell count. For surface export and quantification, the
valid enrichment field is smoothed again using the configured
`surface_smoothing_sigma_voxels_zyx`; a voxel remains quantifiable when it is
valid and the smoothed validity weight is greater than \(10^{-6}\).

Nested gene hotspot masks are defined by empirical quantiles of the smoothed
valid enrichment field. The implemented hotspot levels are

$$
k \in \{\mathrm{top15}, \mathrm{top10}, \mathrm{top05}\},
$$

with thresholds at the 85th, 90th, and 95th percentiles, respectively:

$$
\tau_{g,\mathrm{top15}} = Q_{0.85}(E_g), \quad
\tau_{g,\mathrm{top10}} = Q_{0.90}(E_g), \quad
\tau_{g,\mathrm{top05}} = Q_{0.95}(E_g).
$$

For a hotspot level \(k\), the binary hotspot mask is

$$
H_{g,k}(v) =
\mathbf{1}\{v \in V_{\mathrm{quant}}\}
\mathbf{1}\{E_g(v) \ge \tau_{g,k}\},
$$

which corresponds to the implemented mask
`quant_valid & (smooth_field >= threshold)`.

For each gene \(g\), hotspot level \(k\), and structure \(s_j\), HistoSeg
computes voxel-overlap metrics and SDF metrics. The overlap count is

$$
O_{g,k,j} = \sum_v H_{g,k}(v) M_j(v),
$$

and the fraction of the gene hotspot inside the structure is

$$
\mathrm{fraction\_of\_gene\_in\_structure}_{g,k,j}
=
\begin{cases}
\frac{O_{g,k,j}}{\sum_v H_{g,k}(v)}, & \sum_v H_{g,k}(v) > 0,\\
0, & \sum_v H_{g,k}(v)=0.
\end{cases}
$$

The complementary structure-coverage fraction is

$$
\mathrm{fraction\_of\_structure\_covered\_by\_gene}_{g,k,j}
=
\begin{cases}
\frac{O_{g,k,j}}{\sum_v M_j(v)}, & \sum_v M_j(v) > 0,\\
0, & \sum_v M_j(v)=0.
\end{cases}
$$

For the SDF distribution, HistoSeg evaluates the signed-distance field only at
hotspot voxels:

$$
\mathcal{D}_{g,k,j} =
\{D_j(v) : H_{g,k}(v)=1\}.
$$

The reported signed-distance summaries are the minimum, median, mean, maximum,
5th percentile, and 95th percentile of \(\mathcal{D}_{g,k,j}\). The unsigned
outside-distance summaries are computed from

$$
\max(D_j(v), 0)
$$

over the same hotspot voxels. The structure-inside fraction reported in the
distance table is

$$
\mathrm{fraction\_inside\_structure}_{g,k,j}
=
\frac{\left|\{v \in H_{g,k}: D_j(v) < 0\}\right|}
{\left|H_{g,k}\right|},
$$

and the touching-or-inside fraction is

$$
\mathrm{fraction\_touching\_or\_inside\_structure}_{g,k,j}
=
\frac{\left|\{v \in H_{g,k}: D_j(v) \le 0\}\right|}
{\left|H_{g,k}\right|}.
$$

Empty masks are handled explicitly. If the hotspot mask has zero voxels, or if
the structure mask is empty, HistoSeg reports the hotspot voxel count, sets
`fraction_inside_structure` and `fraction_touching_or_inside_structure` to
zero, and returns `NaN` for all signed and unsigned distance distribution
statistics. This contract prevents empty structure or empty hotspot cases from
being represented as infinite distances.

Cross-gene module matrices are then assembled from per-gene relationship
outputs. The `fraction_inside` matrix uses
`fraction_inside_structure`, the `overlap_fraction` matrix uses
`fraction_of_gene_in_structure`, and the `signed_distance` matrix uses
`median_signed_distance_um` for each hotspot and structure.

## Topology-Aware Multi-Slice Alignment

HistoSeg reconstructs a multi-slice contour stack by processing Xenium slices in
numeric order. Slice 1 is used as the reference slice. Each later slice is first
hard-aligned to the previously accepted aligned slice with a similarity
transform. The hard alignment is accepted only when the union intersection over
union (IoU) between fixed and moving contours does not decrease; otherwise, the
identity transform is retained.

After hard alignment, HistoSeg can apply a conservative soft alignment. The
moving contours are represented by structure-specific boundaries, and boundary
landmarks are generated only between matching structure labels. In the default
nearest-projection mode, moving-boundary points are sampled every
`sampling_distance_um`, projected to the corresponding fixed boundary, and kept
only if the source-to-target distance is not greater than
`max_landmark_distance_um`.

When `landmark_normal_weight_um > 0`, HistoSeg uses the normal-aware landmark
matching path. Fixed-boundary candidate points are sampled at
`landmark_candidate_spacing_um` when provided; otherwise the candidate spacing
is

$$
\delta_c = \max\left(\min\left(\frac{\delta_s}{2}, 25\right), 10^{-6}\right),
$$

where \(\delta_s\) is `sampling_distance_um`. Boundary normals are estimated
from local tangents using a configurable normal step; if no step is provided,
the step is one half of the relevant sampling distance, lower-bounded by
\(10^{-6}\). For a moving sample point \(p_i\) with normal \(n_i\), HistoSeg
queries the \(K\) nearest fixed candidates with a k-d tree, where
\(K=\min(\texttt{landmark\_candidate\_count}, N_c)\). Candidate \(q_l\) is
scored by

$$
C(i,l) =
\|p_i-q_l\|_2
+ \lambda \left(1-\left|n_i^\top m_l\right|\right),
$$

where \(m_l\) is the candidate normal and
\(\lambda=\texttt{landmark\_normal\_weight\_um}\). If either normal is not
finite, the normal penalty is not applied for that candidate. The lowest-cost
candidate is retained when its Euclidean distance is no greater than
`max_landmark_distance_um`.

Landmarks are limited per structure after matching. The general limit is
`landmarks_per_structure`; for the configured diagnostic structure,
`diagnostic_structure_landmarks` can override this value. If limiting is
required, HistoSeg keeps evenly spaced rows from the generated landmark table.
At least three boundary landmarks are required for soft alignment.

The soft warp is fit as a displacement field using
`scipy.interpolate.RBFInterpolator` with the configured kernel, which defaults
to `thin_plate_spline`. Let source landmarks be \(p_i \in \mathbb{R}^2\) and
target landmarks be \(q_i \in \mathbb{R}^2\). HistoSeg normalizes coordinates by
the landmark mean \(c\) and maximum coordinate range \(a\):

$$
\tilde{p}_i = \frac{p_i-c}{a},
\qquad
\tilde{d}_i = \frac{q_i-p_i}{a}.
$$

The interpolator is fit to predict normalized displacement
\(\tilde{u}(\tilde{x})\), with the configured smoothing and optional neighbor
limit. The physical displacement and warp are

$$
u(x) = a \tilde{u}\left(\frac{x-c}{a}\right),
\qquad
T(x) = x + u(x).
$$

Before fitting, HistoSeg adds eight zero-displacement anchor landmarks around a
padded bounding box enclosing the fixed and moving records. These anchors
stabilize the displacement field outside the matched boundaries while preserving
the same target coordinate as the source coordinate.

Warped geometries are validated after applying the TPS displacement. Invalid
geometries are repaired with `buffer(0)` when this yields a valid, non-empty
geometry; otherwise they are counted as invalid. HistoSeg also performs a
topology quality-control grid over the moving union. Grid vertices are warped,
and each grid cell whose center lies inside the moving geometry is assigned a
signed area ratio

$$
r_c =
\frac{A(T(c_1), T(c_2), T(c_3), T(c_4))}
{A_0},
$$

where \(A_0\) is the original grid-cell area and the numerator is the signed
area of the warped quadrilateral. A cell is considered folded when
\(r_c \le 0\), compressed when \(r_c\) is below `topology_min_area_ratio`, and
expanded when \(r_c\) exceeds `topology_max_area_ratio`. The topology check is
valid only when no folded, compressed, or expanded cells are detected.

For full stack reconstruction, the soft-aligned contour is accepted only when
three conditions are simultaneously satisfied: union IoU after soft alignment is
not lower than the hard-aligned IoU, the topology check is valid, and no warped
geometry remains invalid. If any condition fails, HistoSeg keeps the
hard-aligned contour for that slice. The accepted aligned contours define the
3D stack manifest and are subsequently used for contour-point export, optional
surface visualization, and the structure-mask rasterization used by SDF-based
gene quantification.
