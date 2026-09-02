#!/usr/bin/env python3
"""Read-only biological-anchor analysis for ablation-03.

The expression atlas and signature RDS files are never written.  Outputs are
created only below ``test/ablation-03``.  This Python fallback is used when the
host does not provide R/Rscript; the calculations mirror the planned
cross-cohort anchor readout (per-cohort z-scored signature means and paired
Direct/d1 neighbour comparisons).
"""
from __future__ import annotations

import json
import hashlib
import math
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import rdata


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "test" / "ablation-03" / "tmp" / "ablation-biology"
FIG = ROOT / "test" / "ablation-03" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

SIG_RDS = Path(os.getenv(
    "CCS_GENE_SIGNATURE_RDS",
    r"E:\RCloud\database\Signature\report\GeneSignature-HWB.rds",
))
FULL_RDS = Path(os.getenv("CCS_FULL_EXPRESSION_RDS", "")) if os.getenv("CCS_FULL_EXPRESSION_RDS") else None
CACHE_RDS = Path(os.getenv("CCS_BIOLOGY_CACHE_RDS", str(OUT / "expression-anchor-cache.rds")))
RETRIEVAL_RDS = ROOT / "test" / "ablation-03" / "tmp" / "ablation-experiment" / "retrieval.rds"
MANIFEST_RDS = ROOT / "test" / "ablation-03" / "tmp" / "ablation-experiment" / "manifest.rds"


def _scalar(value):
    if isinstance(value, np.ndarray):
        if value.size == 1:
            return _scalar(value.reshape(-1)[0])
        return [_scalar(v) for v in value.reshape(-1)]
    if isinstance(value, (np.generic,)):
        return value.item()
    return value


def _flatten_genes(value):
    genes = []
    if isinstance(value, dict):
        for child in value.values():
            genes.extend(_flatten_genes(child))
    elif isinstance(value, np.ndarray):
        genes.extend(str(v) for v in value.reshape(-1))
    elif isinstance(value, (list, tuple)):
        for child in value:
            genes.extend(_flatten_genes(child))
    return sorted({g for g in genes if g and g != "nan"})


def _find_signature(signatures, predicate):
    for family, entries in signatures.items():
        if not isinstance(entries, dict):
            continue
        for name, value in entries.items():
            name_s = str(name)
            if predicate(name_s):
                genes = _flatten_genes(value)
                if genes:
                    return str(family), name_s, genes
    return None


def _zscore_score(values, gene_idx):
    sub = np.asarray(values[gene_idx, :], dtype=float)
    mu = np.nanmean(sub, axis=1, keepdims=True)
    sd = np.nanstd(sub, axis=1, keepdims=True)
    sd[~np.isfinite(sd) | (sd == 0)] = np.nan
    z = (sub - mu) / sd
    return np.nanmean(z, axis=0)


def _bootstrap_mean(values, rng, n=500):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return (math.nan, math.nan, math.nan)
    draws = values[rng.integers(0, values.size, size=(n, values.size))].mean(axis=1)
    return float(values.mean()), float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def _svg_to_deliverables(svg_text, stem):
    """Write SVG and convert it to the required PDF/JPG previews."""
    svg_path = FIG / f"{stem}.svg"
    svg_path.write_text(svg_text, encoding="utf-8")
    import subprocess
    subprocess.run(["magick", str(svg_path), str(FIG / f"{stem}.pdf")], check=True,
                   capture_output=True)
    subprocess.run(["magick", "-density", "220", str(svg_path), str(FIG / f"{stem}.jpg")],
                   check=True, capture_output=True)


def _bar_svg(labels, values, x_label, stem, width=720, height=420, horizontal=False, series=None):
    margin = {"l": 170 if horizontal else 70, "r": 25, "t": 25, "b": 90}
    plot_w, plot_h = width - margin["l"] - margin["r"], height - margin["t"] - margin["b"]
    colors = ["#4878a8", "#c85c5c"]
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
           '<rect width="100%" height="100%" fill="white"/>', '<g font-family="Arial, sans-serif" font-size="12" fill="#222">']
    if horizontal:
        maxv = max(1.0, max(values) * 1.08)
        for i, (lab, val) in enumerate(zip(labels, values)):
            y = margin["t"] + (i + 0.5) * plot_h / len(labels)
            w = plot_w * val / maxv
            out.append(f'<rect x="{margin["l"]}" y="{y-10:.1f}" width="{w:.1f}" height="20" fill="{colors[0]}"/>')
            out.append(f'<text x="{margin["l"]-8}" y="{y+4:.1f}" text-anchor="end">{lab}</text>')
        out.append(f'<text x="{margin["l"]+plot_w/2:.1f}" y="{height-25}" text-anchor="middle">{x_label}</text>')
    else:
        maxv = max(1.0, max(values) * 1.08)
        for i, (lab, val) in enumerate(zip(labels, values)):
            x = margin["l"] + (i + 0.5) * plot_w / len(labels)
            h = plot_h * val / maxv
            out.append(f'<rect x="{x-18:.1f}" y="{margin["t"]+plot_h-h:.1f}" width="36" height="{h:.1f}" fill="{colors[0]}"/>')
            out.append(f'<text x="{x:.1f}" y="{margin["t"]+plot_h+18}" text-anchor="end" transform="rotate(-25 {x:.1f},{margin["t"]+plot_h+18})">{lab}</text>')
        out.append(f'<text x="{margin["l"]+plot_w/2:.1f}" y="{height-20}" text-anchor="middle">{x_label}</text>')
    out.append('</g></svg>')
    _svg_to_deliverables("".join(out), stem)


def _write_not_estimable(reason, manifest_external=None, retrieval_rows=0):
    """Persist an explicit no-result contract instead of fabricating scores."""
    pd.DataFrame(columns=[
        "tissue", "cohort", "cohort_key", "anchor", "signature",
        "genes_required", "genes_found", "coverage", "sample_count",
        "external_query_cohort", "status", "reason",
    ]).to_csv(OUT / "anchor_coverage.csv", index=False)
    pd.DataFrame(columns=[
        "anchor", "representation", "query_count", "neighbour_pairs",
        "mean_abs_delta", "abs_delta_ci_low", "abs_delta_ci_high",
        "utility", "utility_ci_low", "utility_ci_high", "status", "reason",
    ]).to_csv(OUT / "anchor_utility.csv", index=False)
    pd.DataFrame(columns=[
        "anchor", "d1_minus_direct_abs_delta", "d1_minus_direct_utility",
        "interpretation", "status", "reason",
    ]).to_csv(OUT / "anchor_contrasts.csv", index=False)
    with open(OUT / "anchor_manifest.json", "w", encoding="utf-8") as handle:
        json.dump({
            "status": "not_estimable", "reason": reason,
            "external_cohort_count": int(len(manifest_external or [])),
            "retrieval_rows_top15": int(retrieval_rows),
            "source": str(SIG_RDS),
        }, handle, ensure_ascii=False, indent=2)


def _md5_file(path):
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    manifest = rdata.read_rds(str(MANIFEST_RDS))
    retrieval = rdata.read_rds(str(RETRIEVAL_RDS))
    neighbors = retrieval["neighbors"].copy()
    neighbors.columns = [str(c) for c in neighbors.columns]
    neighbors["neighbor_rank"] = pd.to_numeric(neighbors["neighbor_rank"])
    neighbors = neighbors[neighbors["neighbor_rank"] <= 15].copy()
    manifest_external = [_scalar(v) for v in np.asarray(manifest["external_cohorts"]).reshape(-1)]
    manifest_external = {str(v) for v in manifest_external}
    target_ids = set(neighbors["query_sample"].astype(str)) | set(neighbors["reference_sample"].astype(str))

    if not CACHE_RDS.exists():
        _write_not_estimable("expression_anchor_cache_missing", manifest_external, len(neighbors))
        raise RuntimeError(f"expression-anchor cache is missing: {CACHE_RDS}")
    try:
        cache = rdata.read_rds(str(CACHE_RDS))
    except Exception as exc:
        _write_not_estimable(f"expression_anchor_cache_unreadable:{type(exc).__name__}", manifest_external,
                             len(neighbors))
        raise RuntimeError("expression-anchor cache is unreadable") from exc
    if int(_scalar(cache.get("schema_version", -1))) != 1 or str(_scalar(cache.get("status", ""))) != "complete":
        raise RuntimeError("expression-anchor cache has unsupported or incomplete schema")
    target_hash = hashlib.md5("\n".join(sorted(target_ids)).encode("utf-8")).hexdigest()
    if str(_scalar(cache.get("sample_key_hash", ""))) != target_hash:
        raise RuntimeError("expression-anchor cache sample-key hash mismatch; rebuild the cache")
    if FULL_RDS is not None and FULL_RDS.exists():
        source = cache.get("source", {})
        expected_md5 = str(_scalar(source.get("md5", ""))) if isinstance(source, dict) else ""
        if expected_md5 and _md5_file(FULL_RDS) != expected_md5:
            raise RuntimeError("expression-anchor cache source hash mismatch; rebuild the cache")
    selected = {}
    raw_anchors = cache.get("anchors", {})
    for label, genes in raw_anchors.items():
        selected[str(label)] = {"genes": _flatten_genes(genes), "name": str(label), "family": "frozen_manifest"}
    if len(selected) < 3:
        raise RuntimeError("Fewer than three frozen signatures are present in expression-anchor cache")

    # Score only samples used by the frozen 02 retrieval table.  Each cohort is
    # standardized separately, preventing platform-wide expression shifts from
    # becoming a biological signal.
    scores = {anchor: {} for anchor in selected}
    coverage_rows = []
    raw_cohorts = cache.get("cohorts", {})
    cohort_items = raw_cohorts.items() if isinstance(raw_cohorts, dict) else enumerate(raw_cohorts)
    for _, arr in cohort_items:
        if not isinstance(arr, dict):
            continue
        tissue = str(_scalar(arr.get("tissue", "")))
        cohort = str(_scalar(arr.get("cohort", "")))
        sample_ids = np.asarray([str(_scalar(v)) for v in np.asarray(arr.get("sample_id", [])).reshape(-1)])
        gene_ids = np.asarray([str(_scalar(v)) for v in np.asarray(arr.get("gene_id", [])).reshape(-1)])
        values = np.asarray(arr.get("expression"), dtype=float)
        if values.ndim != 2 or sample_ids.size == 0:
            continue
        for anchor, spec in selected.items():
            wanted = set(spec["genes"])
            idx = np.flatnonzero(np.isin(gene_ids, list(wanted)))
            coverage_rows.append({
                "tissue": str(tissue), "cohort": cohort, "cohort_key": f"{tissue}/{cohort}",
                "anchor": anchor, "signature": spec["name"], "genes_required": len(wanted),
                "genes_found": int(idx.size), "coverage": float(idx.size / len(wanted)),
                "sample_count": int(sample_ids.size),
                "external_query_cohort": f"{tissue}/{cohort}" in manifest_external,
            })
            if idx.size >= 2:
                score = _zscore_score(values, idx)
                scores[anchor].update({sid: float(v) for sid, v in zip(sample_ids, score) if np.isfinite(v)})
    coverage = pd.DataFrame(coverage_rows)

    # Biological utility is an anchor-profile agreement of retrieved neighbours.
    # ``utility = exp(-|z-score difference|)`` is monotone, bounded, and makes
    # lower distance and higher utility directly comparable across anchors.
    rng = np.random.default_rng(20260830)
    result_rows = []
    for anchor in selected:
        anchor_scores = scores[anchor]
        n = neighbors.copy()
        n["q_score"] = n["query_sample"].astype(str).map(anchor_scores)
        n["r_score"] = n["reference_sample"].astype(str).map(anchor_scores)
        n = n.dropna(subset=["q_score", "r_score"])
        n["abs_delta"] = (n["q_score"] - n["r_score"]).abs()
        n["utility"] = np.exp(-n["abs_delta"])
        grouped = n.groupby(["representation", "query_sample"], as_index=False).agg(
            mean_abs_delta=("abs_delta", "mean"), utility=("utility", "mean"),
            neighbour_pairs=("utility", "size"), query_cohort=("query_cohort", "first"),
        )
        for rep, g in grouped.groupby("representation"):
            mean_u, lo_u, hi_u = _bootstrap_mean(g["utility"], rng)
            mean_d, lo_d, hi_d = _bootstrap_mean(g["mean_abs_delta"], rng)
            result_rows.append({
                "anchor": anchor, "representation": rep, "query_count": int(len(g)),
                "neighbour_pairs": int(g["neighbour_pairs"].sum()),
                "mean_abs_delta": mean_d, "abs_delta_ci_low": lo_d, "abs_delta_ci_high": hi_d,
                "utility": mean_u, "utility_ci_low": lo_u, "utility_ci_high": hi_u,
            })
    utility = pd.DataFrame(result_rows)
    direct = utility[utility.representation == "Direct-GSClassifier"].set_index("anchor")
    d1 = utility[utility.representation == "Cohort-d1"].set_index("anchor")
    contrasts = []
    for anchor in sorted(set(direct.index) & set(d1.index)):
        contrasts.append({
            "anchor": anchor,
            "d1_minus_direct_abs_delta": float(d1.loc[anchor, "mean_abs_delta"] - direct.loc[anchor, "mean_abs_delta"]),
            "d1_minus_direct_utility": float(d1.loc[anchor, "utility"] - direct.loc[anchor, "utility"]),
            "interpretation": "d1 higher utility" if d1.loc[anchor, "utility"] > direct.loc[anchor, "utility"] else "Direct higher utility",
        })
    contrasts = pd.DataFrame(contrasts)

    coverage.to_csv(OUT / "anchor_coverage.csv", index=False)
    utility.reset_index(drop=True).to_csv(OUT / "anchor_utility.csv", index=False)
    contrasts.to_csv(OUT / "anchor_contrasts.csv", index=False)
    with open(OUT / "anchor_manifest.json", "w", encoding="utf-8") as handle:
        json.dump({"anchors": selected, "external_cohort_count": len(manifest_external),
                   "retrieval_rows_top15": int(len(neighbors)), "source": str(SIG_RDS)}, handle, ensure_ascii=False, indent=2)

    # Publication-oriented figures: no in-panel titles; captions live in HTML.
    cov = coverage.groupby("anchor", as_index=False)["coverage"].mean().sort_values("coverage")
    _bar_svg(cov["anchor"].tolist(), cov["coverage"].tolist(),
             "Mean gene coverage across evaluated cohorts", "figure-01-anchor-coverage",
             horizontal=True)
    piv = utility.pivot(index="anchor", columns="representation", values="utility")
    direct_values = piv.get("Direct-GSClassifier", pd.Series(index=piv.index, dtype=float)).fillna(0)
    d1_values = piv.get("Cohort-d1", pd.Series(index=piv.index, dtype=float)).fillna(0)
    # A compact paired SVG is sufficient for the HTML preview; tabular values
    # remain the authoritative result and are written above.
    labels = piv.index.tolist()
    combined = ((direct_values + d1_values) / 2).tolist()
    _bar_svg(labels, combined, "Mean paired anchor utility (Direct/d1)",
             "figure-02-biological-utility", horizontal=False)

    print(json.dumps({"anchors": list(selected), "coverage_rows": len(coverage),
                      "utility_rows": len(utility), "contrast_rows": len(contrasts),
                      "output": str(OUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
