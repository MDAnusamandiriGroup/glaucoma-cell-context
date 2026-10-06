#!/usr/bin/env python3
"""Single entry point; verified inputs, deterministic analysis, no-clobber outputs."""
import argparse
import hashlib
import json
import os
from pathlib import Path

from src.workflow import acquire, run_analysis, sha256, software_versions, write_new
from src.reporting import render, make_report


def fingerprint(value: dict) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()[:16]


def validate_manifest(directory: Path, filename: str, expected: str) -> bool:
    manifest_path = directory / filename
    if not manifest_path.exists():
        return False
    manifest = json.loads(manifest_path.read_text())
    if manifest["fingerprint"] != expected:
        raise ValueError(f"Provenance mismatch: {manifest_path}")
    for name, digest in manifest["outputs"].items():
        path = directory / name
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"Incomplete/changed result: {path}. Preserve it for diagnosis.")
    return True


def complete_manifest(directory: Path, filename: str, value: dict) -> None:
    value["outputs"] = {str(p.relative_to(directory)): sha256(p) for p in sorted(directory.rglob("*")) if p.is_file() and p.name != filename}
    write_new(directory / filename, json.dumps(value, indent=2, sort_keys=True) + "\n")


def unused_attempt(base: Path, manifest_name: str, expected: str, rebuild: bool = False) -> tuple[Path, bool]:
    # Preserve an interrupted attempt; recompute only this small unfinished stage.
    if not rebuild:
        for attempt in range(100):
            candidate = base if attempt == 0 else Path(str(base) + f"_attempt{attempt + 1}")
            if candidate.exists() and validate_manifest(candidate, manifest_name, expected):
                return candidate, True
    for attempt in range(100):
        candidate = base if attempt == 0 else Path(str(base) + f"_attempt{attempt + 1}")
        if not candidate.exists():
            return candidate, False
    raise RuntimeError("Too many incomplete stage attempts; inspect existing directories.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--offline", action="store_true", help="Require cached source inputs; never access the network.")
    parser.add_argument("--rebuild-figures", action="store_true", help="Create a new figure attempt from validated saved analysis; preserve every earlier figure file.")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    config = json.loads((root / "config.json").read_text())
    inputs = acquire(root, config, args.offline)
    versions = software_versions()
    provenance = {"config": config, "inputs": inputs, "software": versions, "workflow_sha256": sha256(root / "src/workflow.py")}
    analysis_hash = fingerprint(provenance)
    result_root = root / "results"
    result_root.mkdir(exist_ok=True)
    lock = result_root / ".pipeline.lock"
    try:
        lock.mkdir()
    except FileExistsError:
        raise RuntimeError("Another analysis may be running. If a run was interrupted, check that it has stopped before removing results/.pipeline.lock.")
    try:
        write_new(lock / "owner.json", json.dumps({"pid": os.getpid(), "entry_point": str(Path(__file__).resolve())}) + "\n")
        analysis_dir, reusable = unused_attempt(result_root / ("analysis_" + analysis_hash), "analysis_manifest.json", analysis_hash)
        if reusable:
            print(f"Reusing validated analysis: {analysis_dir.name}", flush=True)
        else:
            print(f"Running analysis: {analysis_dir.name}", flush=True)
            analysis_dir.mkdir()
            run_analysis(root, analysis_dir, config)
            complete_manifest(analysis_dir, "analysis_manifest.json", {"fingerprint": analysis_hash, "provenance": provenance, "status": "complete"})
        figure_hash = fingerprint({"analysis": analysis_hash, "reporting_sha256": sha256(root / "src/reporting.py"), "matplotlib": versions["matplotlib"]})
        figure_dir, plot_reusable = unused_attempt(result_root / ("figures_" + figure_hash), "figure_manifest.json", figure_hash, rebuild=args.rebuild_figures)
        if plot_reusable:
            print(f"Reusing validated figures: {figure_dir.name}", flush=True)
        else:
            print(f"Rendering figures: {figure_dir.name}", flush=True)
            render(analysis_dir, figure_dir)
            write_new(figure_dir / "SUMMARY.md", make_report(analysis_dir, figure_dir))
            complete_manifest(figure_dir, "figure_manifest.json", {"fingerprint": figure_hash, "analysis": analysis_hash, "status": "complete"})
        print(f"Completed report: {figure_dir / 'SUMMARY.md'}", flush=True)
        print(f"Completed tables: {analysis_dir / 'tables'}", flush=True)
    finally:
        (lock / "owner.json").unlink(missing_ok=True)
        lock.rmdir()


if __name__ == "__main__":
    main()
