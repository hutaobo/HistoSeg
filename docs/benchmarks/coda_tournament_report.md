# CODA-Inspired Tournament Backend Report

Date: 2026-05-04

Branch: `feature/coda-backend`

This report records the validation run for the CODA-inspired hard-alignment backend in HistoSeg 3D. The implementation is CODA-inspired, not a full CODA reimplementation. It uses tissue-mask Radon rotation and phase-correlation translation as an image-derived hard seed, then evaluates that seed against the existing contour similarity seed before the topology-safe contour TPS refinement.

Credit:

- Kiemen et al., "CODA: quantitative 3D reconstruction of large tissues at cellular resolution," Nature Methods 19, 1490-1499, 2022. DOI: <https://doi.org/10.1038/s41592-022-01650-9>
- CODA methodology page: <https://labs.pathology.jhu.edu/kiemen/coda-3d/>

## Backend Design

For `--registration-backend coda-image`, HistoSeg now uses a tournament-style hard-alignment stage:

1. Run the legacy contour similarity hard seed.
2. Run the CODA-inspired Radon/phase hard seed.
3. Compare candidates by `union_iou_after_hard`.
4. Promote the higher-IoU candidate to the canonical `moving_hard_aligned.geojson`.
5. Pass the selected hard seed into the existing topology-safe contour TPS refinement.

Tie-breaking is conservative: if the candidates have equal hard IoU, the contour seed is selected. The final hard summary stores both candidates and the selected seed backend.

## Synthetic Stress Tests

The stress tests use synthetic asymmetric tissue masks with large rotations and translations.

| Scenario | Trials | Contour seed success | CODA seed success | Mean IoU gain over contour | Mean full-circle angle error |
| --- | ---: | ---: | ---: | ---: | ---: |
| 75-degree rotation | 3 | 0.0 | 1.0 | +0.5572 | 0.0 deg |
| 120-degree rotation before ambiguity fix | 1 | 0.0 | 0.0 | -0.3640 | 180.0 deg |
| 120-degree rotation after ambiguity fix | 1 | 0.0 | 1.0 | +0.5535 | 0.0 deg |

Interpretation: Radon/phase seeding is valuable as a rescue path for badly initialized, strongly rotated, or half-turn ambiguous tissue masks. The half-turn ambiguity was fixed by evaluating the Radon angle candidate and its 180-degree alternative in physical mask space.

## Polyp 32-Slice Tournament Run

Dataset:

- `Y:\long\spatialpathologist\3D aligment\polyp`
- 32 slices, 31 adjacent slice pairs.
- Baseline output: `Y:\long\spatialpathologist\3D aligment\polyp\histoseg_3d_reconstruction`
- Tournament output: `Y:\long\spatialpathologist\3D aligment\polyp\histoseg_3d_reconstruction_coda_tournament`

Command shape:

```bash
histoseg-3d reconstruct-stack \
  --xenium-root "Y:\long\spatialpathologist\3D aligment\polyp" \
  --out-dir "Y:\long\spatialpathologist\3D aligment\polyp\histoseg_3d_reconstruction_coda_tournament" \
  --segmentation-strategy "Y:\long\spatialpathologist\3D aligment\polyp\contour for alignment\segmentationstrategy.txt" \
  --merged-h5ad "Y:\long\spatialpathologist\3D aligment\polyp\pdc_merge_leiden\polyp_32samples_min3_count5_leiden_20260501_processed_leiden.h5ad" \
  --merged-cluster-column leiden_1_0 \
  --registration-backend coda-image \
  --hard-alignment-maxiter 40 \
  --sampling-distance-um 50 \
  --max-landmark-distance-um 180 \
  --landmarks-per-structure 260 \
  --diagnostic-structure-landmarks 620 \
  --rbf-neighbors 96 \
  --rbf-smoothing 0.0001 \
  --diagnostic-structure "Structure 5" \
  --point-sample-distance-um 80 \
  --voxel-size-um 80 \
  --no-alignment-preview \
  --overwrite
```

Tournament selection:

| Selected hard seed | Pair count |
| --- | ---: |
| `contour-tps` | 31 |
| `coda-image` | 0 |

Candidate outcomes:

| Outcome | Pair count |
| --- | ---: |
| CODA candidate better than contour candidate | 0 |
| Contour candidate better than CODA candidate | 31 |
| Tied candidates | 0 |

Alignment metrics:

| Metric | Baseline contour TPS | CODA tournament | Delta |
| --- | ---: | ---: | ---: |
| Mean hard IoU | 0.93863 | 0.93844 | -0.00019 |
| Mean soft IoU | 0.94799 | 0.94651 | -0.00148 |
| Minimum hard IoU | 0.90779 | 0.91207 | +0.00428 |
| Minimum soft IoU | 0.91360 | 0.91962 | +0.00602 |

Topology guard:

| Check | Result |
| --- | --- |
| Soft TPS accepted | 31 / 31 |
| Geometry valid | 31 / 31 |
| Topology valid | 31 / 31 |
| Folded grid cells | 0 |
| Area ratio range | 0.6862 to 1.3499 |

Closest CODA near misses:

| Moving sample | Contour candidate IoU | CODA candidate IoU | CODA minus contour |
| --- | ---: | ---: | ---: |
| `A079-C-008_3` | 0.94150 | 0.93335 | -0.00815 |
| `A079-C-008_2` | 0.96305 | 0.95364 | -0.00941 |
| `A079-C-008_4` | 0.91889 | 0.88924 | -0.02965 |

## Interpretation

On the real 32-slice polyp reconstruction, the original contour similarity seed remained the strongest candidate for every adjacent slice pair. This is not a failure of the CODA-inspired backend. It confirms that the tournament layer behaves conservatively on already well-conditioned serial sections: CODA is evaluated, provenance is recorded, and the system declines to use it when it would reduce hard IoU.

The synthetic stress tests establish the complementary role of CODA-inspired image seeding: it rescues large-rotation or half-turn initialization failures that the contour seed cannot reliably recover from. The polyp run establishes the regression behavior: when contour geometry is already near optimal, the tournament keeps the baseline path and preserves the topology-safe TPS acceptance chain.

Operational conclusion:

- Use `contour-tps` as the default production backend.
- Use `coda-image` when input sections may have large unknown rotations, flipped orientation, poor manual prealignment, or inconsistent acquisition pose.
- Treat `selected_hard_seed_backend == "coda-image"` as a useful QC flag for difficult slice pairs.
- Treat `selected_hard_seed_backend == "contour-tps"` across all pairs as evidence that contour geometry was already sufficient for the stack.
