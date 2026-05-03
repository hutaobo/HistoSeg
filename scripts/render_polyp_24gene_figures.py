from __future__ import annotations

import argparse
import shutil
from pathlib import Path


STRUCTURE_COLORS = {
    "Structure 1": "#2dd4bf",
    "Structure 2": "#60a5fa",
    "Structure 3": "#f59e0b",
    "Structure 4": "#a78bfa",
    "Structure 5": "#fb7185",
}

GREM1_SURFACES = {
    "top15": {"color": "#f97316", "opacity": 0.22},
    "top10": {"color": "#dc2626", "opacity": 0.42},
    "top05": {"color": "#7f1d1d", "opacity": 0.82},
}


def _import_pyvista(start_xvfb: bool, off_screen: bool):
    try:
        import pyvista as pv
    except ImportError as exc:
        raise RuntimeError(
            "PyVista is required only for tutorial figure rendering. "
            "Install it with `pip install pyvista` and rerun this script."
        ) from exc
    pv.OFF_SCREEN = off_screen
    if start_xvfb:
        try:
            pv.start_xvfb()
        except Exception:
            # Windows and EGL-enabled Linux runners commonly do not need xvfb.
            pass
    return pv


def _read_mesh(pv, path: Path, *, z_scale: float):
    mesh = pv.read(path)
    if z_scale != 1.0:
        mesh = mesh.copy(deep=True)
        mesh.points[:, 2] *= z_scale
    return mesh


def _mesh_bounds(meshes):
    import numpy as np

    bounds = [mesh.bounds for mesh in meshes]
    xmin = min(bound[0] for bound in bounds)
    xmax = max(bound[1] for bound in bounds)
    ymin = min(bound[2] for bound in bounds)
    ymax = max(bound[3] for bound in bounds)
    zmin = min(bound[4] for bound in bounds)
    zmax = max(bound[5] for bound in bounds)
    center = np.array([(xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2], dtype=float)
    size = max(xmax - xmin, ymax - ymin, zmax - zmin)
    return (xmin, xmax, ymin, ymax, zmin, zmax), center, size


def _camera(meshes, preset: str):
    bounds, center, size = _mesh_bounds(meshes)
    if preset == "dorsal":
        position = (center[0], center[1] - size * 0.10, center[2] + size * 1.75)
        view_up = (0, 1, 0)
    elif preset == "structure5_focus":
        position = (center[0] - size * 0.85, center[1] - size * 0.95, center[2] + size * 0.82)
        view_up = (0, 0, 1)
    else:
        position = (center[0] + size * 0.92, center[1] - size * 1.10, center[2] + size * 0.78)
        view_up = (0, 0, 1)
    return position, tuple(center), view_up


def _new_plotter(pv, args):
    plotter = pv.Plotter(
        off_screen=args.off_screen,
        window_size=tuple(args.window_size),
        lighting="three lights",
        border=False,
    )
    plotter.set_background(args.background)
    return plotter


def _add_structures(plotter, structure_meshes, *, opacity: float, show_edges: bool):
    for name, mesh in structure_meshes.items():
        plotter.add_mesh(
            mesh,
            color=STRUCTURE_COLORS.get(name, "#94a3b8"),
            opacity=opacity,
            smooth_shading=True,
            specular=0.18,
            roughness=0.55,
            show_edges=show_edges,
            label=name,
        )


def _add_grem1(plotter, grem1_meshes, levels):
    for level in levels:
        if level not in grem1_meshes:
            continue
        style = GREM1_SURFACES[level]
        plotter.add_mesh(
            grem1_meshes[level],
            color=style["color"],
            opacity=style["opacity"],
            smooth_shading=True,
            specular=0.35,
            roughness=0.36,
            label=f"GREM1 {level}",
        )


def _render(plotter, meshes, out_path: Path, preset: str, transparent: bool):
    plotter.camera_position = _camera(meshes, preset)
    plotter.camera.zoom(1.08)
    plotter.enable_anti_aliasing("ssaa")
    plotter.add_axes(line_width=2, labels_off=False)
    plotter.screenshot(out_path, transparent_background=transparent)
    plotter.close()


def render_figures(args) -> list[Path]:
    pv = _import_pyvista(start_xvfb=not args.no_xvfb, off_screen=args.off_screen)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    structure_meshes = {
        f"Structure {index}": _read_mesh(
            pv,
            args.stack_root / "meshes" / f"Structure_{index}.ply",
            z_scale=args.z_scale,
        )
        for index in range(1, 6)
    }
    grem1_meshes = {
        level: _read_mesh(
            pv,
            args.gene_density_dir / "isosurfaces" / f"{args.gene}_enrichment_{level}_isosurface.ply",
            z_scale=args.z_scale,
        )
        for level in GREM1_SURFACES
        if (args.gene_density_dir / "isosurfaces" / f"{args.gene}_enrichment_{level}_isosurface.ply").exists()
    }
    all_meshes = list(structure_meshes.values()) + list(grem1_meshes.values())
    outputs: list[Path] = []

    if "anatomy" in args.figures or "all" in args.figures:
        out = args.out_dir / "polyp_3d_structure_surfaces.png"
        plotter = _new_plotter(pv, args)
        _add_structures(plotter, structure_meshes, opacity=0.44, show_edges=False)
        plotter.add_text("Polyp 3D contour surfaces", font_size=18, color=args.text_color)
        _render(plotter, list(structure_meshes.values()), out, args.camera_preset, args.transparent_background)
        outputs.append(out)

    if "grem1" in args.figures or "all" in args.figures:
        out = args.out_dir / f"{args.gene}_nested_3d_hotspot_surfaces.png"
        plotter = _new_plotter(pv, args)
        _add_structures(plotter, structure_meshes, opacity=0.13, show_edges=False)
        _add_grem1(plotter, grem1_meshes, ["top15", "top10", "top05"])
        plotter.add_text(f"{args.gene} nested enrichment hotspots", font_size=18, color=args.text_color)
        _render(plotter, all_meshes, out, args.camera_preset, args.transparent_background)
        outputs.append(out)

    if "structure5" in args.figures or "all" in args.figures:
        out = args.out_dir / f"{args.gene}_structure5_focus.png"
        plotter = _new_plotter(pv, args)
        _add_structures(plotter, {"Structure 5": structure_meshes["Structure 5"]}, opacity=0.24, show_edges=False)
        if "Structure 3" in structure_meshes:
            plotter.add_mesh(
                structure_meshes["Structure 3"],
                color=STRUCTURE_COLORS["Structure 3"],
                opacity=0.08,
                smooth_shading=True,
                label="Structure 3",
            )
        _add_grem1(plotter, grem1_meshes, ["top15", "top10", "top05"])
        plotter.add_text(f"{args.gene} hotspot embedded in Structure 5", font_size=18, color=args.text_color)
        _render(plotter, [structure_meshes["Structure 5"], *grem1_meshes.values()], out, "structure5_focus", args.transparent_background)
        outputs.append(out)

    if args.batch_dir:
        for name in [
            "fraction_inside_top05_spatial_clustermap.png",
            "fraction_inside_top10_spatial_clustermap.png",
            "fraction_inside_top15_spatial_clustermap.png",
            "signed_distance_top05_spatial_clustermap.png",
        ]:
            src = args.batch_dir / name
            if src.exists():
                dst = args.out_dir / name
                shutil.copy2(src, dst)
                outputs.append(dst)
    return outputs


def _parse_figures(values: list[str]) -> list[str]:
    parsed: list[str] = []
    for value in values:
        parsed.extend(part.strip() for part in value.split(",") if part.strip())
    return parsed or ["all"]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Render Polyp 24-gene HistoSeg 3D tutorial hero figures.")
    parser.add_argument("--stack-root", type=Path, required=True)
    parser.add_argument("--gene", default="GREM1")
    parser.add_argument("--gene-density-dir", type=Path, default=None)
    parser.add_argument("--batch-dir", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=Path("docs/_static/threed/polyp/24gene"))
    parser.add_argument("--figures", nargs="*", default=["all"], help="all, anatomy, grem1, structure5")
    parser.add_argument("--camera-preset", default="oblique", choices=["oblique", "dorsal", "structure5_focus"])
    parser.add_argument("--z-scale", type=float, default=8.0)
    parser.add_argument("--window-size", nargs=2, type=int, default=[2200, 1600])
    parser.add_argument("--background", default="white")
    parser.add_argument("--text-color", default="black")
    parser.add_argument("--transparent-background", action="store_true")
    parser.add_argument("--off-screen", action="store_true", default=True)
    parser.add_argument("--no-xvfb", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.figures = _parse_figures(args.figures)
    if args.gene_density_dir is None:
        args.gene_density_dir = args.stack_root / "gene_overlays" / f"{args.gene}_density"
    outputs = render_figures(args)
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
