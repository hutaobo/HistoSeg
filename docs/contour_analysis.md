# Contour Analysis

Contour Analysis is the spatial/cell-coordinate workflow group exposed as `histoseg.contour`.

## Workflows

- Pattern1 isoline contours from clustered cell coordinates.
- Multi-structure contour partitioning and Xenium Explorer annotation exports.
- Contour topology analysis for boundary-neighbor and enclosure relationships.

## Python API

```python
from histoseg.contour import Pattern1IsolineConfig, run_pattern1_isoline

cfg = Pattern1IsolineConfig(
    clusters_csv="/path/to/clusters.csv",
    cells_parquet="/path/to/cells.parquet",
    tissue_boundary_csv="/path/to/tissue_boundary.csv",
    out_dir="outputs/pattern1_isoline0p5",
    pattern1_clusters=(10, 23, 19, 27, 14, 20, 25, 26),
)
result = run_pattern1_isoline(cfg)
```

```python
from histoseg.contour import MultiStructureContourConfig, MultiStructureSpec
from histoseg.contour import run_multi_structure_contours

result = run_multi_structure_contours(
    MultiStructureContourConfig(
        clusters_csv="/path/to/clusters.csv",
        cells_parquet="/path/to/cells.parquet",
        out_dir="outputs/multi_structure",
        structures=[
            MultiStructureSpec("Region A", [1, 2, 3]),
            MultiStructureSpec("Region B", [4, 5, 6]),
        ],
    )
)
```

```python
import pandas as pd
from shapely.geometry import box
from histoseg.contour import summarize_contour_topology

contours = pd.DataFrame(
    {
        "contour_id": ["S1_0", "S2_0", "S4_inner"],
        "assigned_structure": ["S1", "S2", "S4"],
        "geometry": [
            box(0, 0, 100, 100),
            box(100, 0, 200, 100),
            box(25, 25, 75, 75),
        ],
    }
)

topology = summarize_contour_topology(
    contours,
    groupby="assigned_structure",
    boundary_tolerance=1.0,
    enclosure_min_fraction=0.95,
)

topology.boundary_overlap
topology.enclosure
topology.contour_summary
```

## CLI

```bash
histoseg-contour pattern1 --clusters-csv clusters.csv --cells-parquet cells.parquet --out-dir outputs/p1 --pattern1-clusters 10,23,19
histoseg-contour multi-structure --clusters-csv clusters.csv --cells-parquet cells.parquet --out-dir outputs/ms --structures-json structures.json
```
