#!/usr/bin/env python3
from __future__ import annotations

import argparse
import fnmatch
import subprocess
import sys


CONFIDENTIAL_GLOBS = (
    "docs/ai_driven_spatial_pathologist_patent_*",
)
CONFIDENTIAL_EXACT = {
    "docs/patent_centralized_bundle.zip",
}
CONFIDENTIAL_DIRS = (
    "docs/patent_centralized_bundle/",
    "docs/patent_figures/",
    "docs/patent_package/",
)


def _normalize_path(path: str) -> str:
    normalized = path.replace("\\", "/").strip()
    while normalized.startswith("./"):
        normalized = normalized[2:]
    return normalized


def is_confidential_path(path: str) -> bool:
    normalized = _normalize_path(path)
    lowered = normalized.lower()

    if any(fnmatch.fnmatch(lowered, pattern.lower()) for pattern in CONFIDENTIAL_GLOBS):
        return True
    if lowered in {item.lower() for item in CONFIDENTIAL_EXACT}:
        return True
    return any(
        lowered == directory[:-1].lower() or lowered.startswith(directory.lower())
        for directory in CONFIDENTIAL_DIRS
    )


def _git_output(*args: str) -> bytes:
    result = subprocess.run(
        ["git", *args],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return result.stdout


def tracked_paths() -> list[str]:
    output = _git_output("ls-files", "-z")
    return [_normalize_path(item) for item in output.decode("utf-8", errors="replace").split("\x00") if item]


def staged_paths() -> list[str]:
    output = _git_output("diff", "--cached", "--name-only", "--diff-filter=ACMR", "-z")
    return [_normalize_path(item) for item in output.decode("utf-8", errors="replace").split("\x00") if item]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fail if confidential patent materials are tracked or staged."
    )
    parser.add_argument(
        "--mode",
        choices=("tracked", "staged"),
        default="tracked",
        help="Check repository tracked files or staged files.",
    )
    parser.add_argument(
        "--tracked",
        action="store_true",
        help="Alias for --mode tracked.",
    )
    parser.add_argument(
        "--staged",
        action="store_true",
        help="Alias for --mode staged.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    mode = args.mode
    if args.tracked and args.staged:
        print("Use only one of --tracked or --staged.", file=sys.stderr)
        return 2
    if args.tracked:
        mode = "tracked"
    if args.staged:
        mode = "staged"

    try:
        paths = tracked_paths() if mode == "tracked" else staged_paths()
    except subprocess.CalledProcessError as exc:
        print(exc.stderr.decode("utf-8", errors="replace"), file=sys.stderr)
        return exc.returncode

    matches = sorted(path for path in paths if is_confidential_path(path))
    if not matches:
        print(f"Confidential patent path check passed for {mode} files.")
        return 0

    print(
        f"Confidential patent material detected in {mode} files. "
        "These paths must stay outside Git:",
        file=sys.stderr,
    )
    for path in matches:
        print(f"  - {path}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
