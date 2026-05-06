# GASTON Fig. 5 / HistoSeg Breast Xenium Benchmark

This harness compares HistoSeg domain contours against the GASTON Xenium
breast cancer Fig. 5 workflow on the 10x Genomics Rep1 dataset.

Default data root:

`Y:\long\10X_datasets\Xenium\Xenium_Breast_Cancer\Xenium_FFPE_Human_Breast_Cancer_Rep1_outs`

Default output roots:

- Local: `D:\histoseg_gaston_brca_fig5_benchmark`
- A100: `/data/taobo.hu/projects/histoseg_gaston_brca_fig5_benchmark`

The default crop matches the GASTON tutorial ADH region:
`1500 < x < 3250`, `2000 < y < 4000`.

## Local Smoke Test

```powershell
python -m reproducibility.benchmarks.gaston_brca_fig5.benchmark `
  --mode smoke `
  --out-dir D:\histoseg_gaston_brca_fig5_benchmark_smoke
```

## Full A100 Run

The full run expects `gaston-spatial`, `glmpca`, PyTorch with CUDA, Scanpy,
and HistoSeg to be importable in the active Python environment.

```bash
python3 -m venv /data/taobo.hu/projects/histoseg_gaston_brca_fig5_benchmark/venv
source /data/taobo.hu/projects/histoseg_gaston_brca_fig5_benchmark/venv/bin/activate
python -m pip install -U pip
python -m pip install numpy pandas scipy scikit-learn matplotlib seaborn shapely pyarrow h5py scanpy
python -m pip install --index-url https://download.pytorch.org/whl/cu121 'torch==2.5.1+cu121'
python -m pip install gaston-spatial glmpca
python -m pip install -e /path/to/HistoSeg

python -m reproducibility.benchmarks.gaston_brca_fig5.benchmark \
  --mode all \
  --xenium-outs /path/to/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs \
  --out-dir /data/taobo.hu/projects/histoseg_gaston_brca_fig5_benchmark \
  --cell-type-url <public_Cell_Barcode_Type_Matrices.csv_url>
```

If the exact public `Cell_Barcode_Type_Matrices.csv` cannot be fetched or
does not align to the Xenium barcodes, the cell-type proportion and
cell-type-specific Fig. 5 panels fail explicitly rather than using proxy labels.

For A100 driver stacks reporting CUDA 12.2, avoid the default PyPI torch wheel
if it resolves to `+cu130`; use the CUDA 12.1 wheel command above and verify
`torch.cuda.is_available()` before launching the 30-seed run.
