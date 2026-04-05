from pathlib import Path

from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib.utils import simpleSplit
from reportlab.pdfgen import canvas


OUT_PDF = Path('/Users/hutaobo/Documents/GitHub/HistoSeg/output/pdf/histoseg_app_summary.pdf')


TITLE = 'HistoSeg: One-Page App Summary'

WHAT_IT_IS = (
    'HistoSeg is a Python toolkit for spatial transcriptomics segmentation and geometry extraction. '
    'In this repo, the main implemented workflow generates Pattern1 isoline contours (default level 0.5) from clustered cell data.'
)

WHO_ITS_FOR = (
    'Primary persona: computational biology and bioinformatics users working with 10x Xenium and GraphClust outputs. '
    'Explicit persona statement: Not found in repo (persona inferred from README/tutorial inputs).'
)

FEATURES = [
    'Aligns clusters.csv with cells.parquet and auto-detects join id and x/y coordinate columns.',
    'Selects target Pattern1 clusters and builds background points from non-target cells.',
    'Optionally generates synthetic background points from tissue_boundary.csv bounding box.',
    'Fits a distance-weighted KNN regressor, predicts on a mesh, and smooths probabilities with a Gaussian filter.',
    'Masks non-tissue grid regions, extracts isolines, and filters loops by minimum cells-inside threshold.',
    'Writes reproducible outputs: params.json, contour .npy arrays, and preview .png image.',
    'Provides optional Hugging Face dataset download helper and an experimental GUI entry point (histoseg-gui).',
]

ARCHITECTURE = [
    'Inputs: clusters.csv + cells.parquet (+ optional tissue_boundary.csv) from local files or HF dataset helper.',
    'Core compute: run_pattern1_isoline in src/histoseg/contours/pattern1_isoline.py.',
    'I/O service: src/histoseg/io/huggingface.py downloads required dataset files for HF wrapper runs.',
    'UI surface: tkinter GUI in src/histoseg/gui/gui_app.py, exposed by console script histoseg-gui in setup.py.',
    'Data flow: ingest and align -> choose target/background -> KNN field -> smoothing/mask -> contour extraction -> saved artifacts.',
    'Backend API service or deployed web service: Not found in repo.',
]

RUN_STEPS = [
    'Install: pip install -U histoseg  (or from source: pip install -e .)',
    'Prepare inputs: clusters.csv, cells.parquet, and optional tissue_boundary.csv.',
    'Run local API: Pattern1IsolineConfig(...) + run_pattern1_isoline(...).',
    'Alternative: run_pattern1_isoline_from_hf(repo_id=...) to pull inputs from a Hugging Face dataset repo.',
    'Optional GUI: run histoseg-gui (README marks this path experimental).',
]


def draw_heading(c: canvas.Canvas, text: str, x: float, y: float, size: int = 13) -> float:
    c.setFont('Helvetica-Bold', size)
    c.drawString(x, y, text)
    return y - 0.20 * inch


def draw_wrapped(c: canvas.Canvas, text: str, x: float, y: float, width: float, size: int = 11, leading: float = 14) -> float:
    c.setFont('Helvetica', size)
    lines = simpleSplit(text, 'Helvetica', size, width)
    for line in lines:
        c.drawString(x, y, line)
        y -= leading
    return y


def draw_bullets(c: canvas.Canvas, items: list[str], x: float, y: float, width: float, size: int = 10.6, leading: float = 13) -> float:
    c.setFont('Helvetica', size)
    bullet_indent = 10
    text_width = width - bullet_indent - 8
    for item in items:
        lines = simpleSplit(item, 'Helvetica', size, text_width)
        if not lines:
            lines = ['']
        c.drawString(x, y, '-')
        c.drawString(x + bullet_indent, y, lines[0])
        y -= leading
        for line in lines[1:]:
            c.drawString(x + bullet_indent, y, line)
            y -= leading
    return y


def build_pdf(path: Path) -> float:
    path.parent.mkdir(parents=True, exist_ok=True)
    c = canvas.Canvas(str(path), pagesize=letter)

    page_w, page_h = letter
    left = 0.65 * inch
    right = 0.65 * inch
    top = page_h - 0.62 * inch
    bottom = 0.62 * inch
    width = page_w - left - right

    y = top
    c.setTitle('HistoSeg One-Page App Summary')

    c.setFont('Helvetica-Bold', 18)
    c.drawString(left, y, TITLE)
    y -= 0.30 * inch

    c.setFont('Helvetica', 9.2)
    c.drawString(left, y, 'Evidence scope: README.md, setup.py, pyproject.toml, and src/histoseg/*.')
    y -= 0.23 * inch

    y = draw_heading(c, 'What It Is', left, y)
    y = draw_wrapped(c, WHAT_IT_IS, left, y, width)
    y -= 0.06 * inch

    y = draw_heading(c, 'Who It Is For', left, y)
    y = draw_wrapped(c, WHO_ITS_FOR, left, y, width)
    y -= 0.06 * inch

    y = draw_heading(c, 'What It Does', left, y)
    y = draw_bullets(c, FEATURES, left, y, width)
    y -= 0.05 * inch

    y = draw_heading(c, 'How It Works (Architecture)', left, y)
    y = draw_bullets(c, ARCHITECTURE, left, y, width)
    y -= 0.05 * inch

    y = draw_heading(c, 'How To Run (Minimal)', left, y)
    y = draw_bullets(c, RUN_STEPS, left, y, width)

    c.setFont('Helvetica-Oblique', 8)
    c.drawRightString(page_w - right, bottom - 0.11 * inch, 'Generated from repository evidence only')

    c.showPage()
    c.save()
    return y - bottom


if __name__ == '__main__':
    delta = build_pdf(OUT_PDF)
    print(f'Wrote: {OUT_PDF}')
    print(f'Layout delta (negative means overflow): {delta:.2f}')
