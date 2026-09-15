#!/usr/bin/env python3
"""Lyso-IP QC organelle coverage: two one-sided Fisher's exact tests asking
whether marker genes clear the Positive gate, or the Negative gate, more
often than the run's own tested background. Optionally repeats against
every curated list, BH-FDR corrected, ranked best-match-first.

Ported from lysoip_qc_framework/apps/scoring/organelle_coverage.py.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

from scipy.stats import fisher_exact, false_discovery_control

DATA_FILE = Path(__file__).resolve().parent / "data" / "organelle_marker_lists.json"


def one_sided_fisher(target, other_a, other_b, bg_target, bg_other_a, bg_other_b):
    return float(fisher_exact([[target, other_a + other_b], [bg_target, bg_other_a + bg_other_b]], alternative="greater")[1])


def compute_coverage(verdicts, marker_genes):
    counts = {(True, "Positive"): 0, (True, "Negative"): 0, (True, "Nonsignificant"): 0,
              (False, "Positive"): 0, (False, "Negative"): 0, (False, "Nonsignificant"): 0}
    for gene, classification in verdicts:
        is_marker = gene in marker_genes
        counts[(is_marker, classification)] += 1

    n_marker_positive = counts[(True, "Positive")]
    n_marker_negative = counts[(True, "Negative")]
    n_marker_nonsignificant = counts[(True, "Nonsignificant")]
    n_background_positive = counts[(False, "Positive")]
    n_background_negative = counts[(False, "Negative")]
    n_background_nonsignificant = counts[(False, "Nonsignificant")]

    n_marker_total = n_marker_positive + n_marker_negative + n_marker_nonsignificant
    n_background_total = n_background_positive + n_background_negative + n_background_nonsignificant

    result = {
        "n_marker_positive": n_marker_positive, "n_marker_negative": n_marker_negative,
        "n_marker_nonsignificant": n_marker_nonsignificant,
        "n_background_positive": n_background_positive, "n_background_negative": n_background_negative,
        "n_background_nonsignificant": n_background_nonsignificant,
        "positive_pvalue": None, "negative_pvalue": None,
        "is_significantly_positive": None, "is_significantly_negative": None,
    }
    if n_marker_total == 0 or n_background_total == 0:
        return result

    result["positive_pvalue"] = one_sided_fisher(
        n_marker_positive, n_marker_negative, n_marker_nonsignificant,
        n_background_positive, n_background_negative, n_background_nonsignificant,
    )
    result["negative_pvalue"] = one_sided_fisher(
        n_marker_negative, n_marker_positive, n_marker_nonsignificant,
        n_background_negative, n_background_positive, n_background_nonsignificant,
    )
    return result


def main():
    parser = argparse.ArgumentParser(description="Organelle marker-list coverage test against this run's verdicts")
    parser.add_argument("--verdict_file", required=True)
    parser.add_argument("--organelle_list", required=True)
    parser.add_argument("--test_all_lists", action="store_true")
    parser.add_argument("--coverage_alpha", type=float, default=0.05)
    parser.add_argument("--output_folder", required=True)
    args = parser.parse_args()

    output_folder = Path(args.output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    # @step: Loading verdicts
    with open(args.verdict_file) as f:
        verdicts = [(row["gene"], row["classification"]) for row in csv.DictReader(f, delimiter="\t")]

    # @step: Loading organelle marker lists
    with open(DATA_FILE) as f:
        marker_lists = json.load(f)
    if args.organelle_list not in marker_lists:
        raise SystemExit(f"Unknown organelle_list {args.organelle_list!r}; available: {sorted(marker_lists)}")

    lists_to_test = list(marker_lists) if args.test_all_lists else [args.organelle_list]

    # @step-if: Testing coverage against every curated list
    rows = []
    for list_name in lists_to_test:
        marker_genes = set(marker_lists[list_name]["genes"])
        coverage = compute_coverage(verdicts, marker_genes)
        coverage["organelle_list"] = list_name
        coverage["is_expected"] = list_name == args.organelle_list
        rows.append(coverage)

    # @step: Applying BH-FDR correction and ranking
    for row in rows:
        row["positive_fdr"] = ""

    if args.test_all_lists:
        tested = [r for r in rows if r["positive_pvalue"] is not None]
        if tested:
            adjusted = false_discovery_control([r["positive_pvalue"] for r in tested], method="bh")
            for row, fdr in zip(tested, adjusted):
                row["positive_fdr"] = float(fdr)
        rows = sorted(tested, key=lambda r: r["positive_pvalue"])

    for row in rows:
        row["is_significantly_positive"] = row["positive_pvalue"] is not None and row["positive_pvalue"] <= args.coverage_alpha
        row["is_significantly_negative"] = row["negative_pvalue"] is not None and row["negative_pvalue"] <= args.coverage_alpha

    fieldnames = [
        "organelle_list", "is_expected", "n_marker_positive", "n_marker_negative", "n_marker_nonsignificant",
        "n_background_positive", "n_background_negative", "n_background_nonsignificant",
        "positive_pvalue", "negative_pvalue", "positive_fdr", "is_significantly_positive", "is_significantly_negative",
    ]
    with open(output_folder / "organelle_coverage.tsv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in fieldnames})

    print(f"Wrote coverage results for {len(rows)} list(s)", file=sys.stderr)
    print("Organelle coverage complete.", file=sys.stderr)


if __name__ == "__main__":
    main()
