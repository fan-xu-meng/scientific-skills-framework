#!/usr/bin/env python3
"""Summarize Pt fcc FFT candidate proportions without overclaiming facets."""

from pathlib import Path
import os
import pandas as pd

OUT = Path(os.environ["TEM_OUTPUT_DIR"]).resolve() / "pt-candidates"
INPUT = OUT / "lattice_spacing_candidates_pending_indexing.csv"
TARGETS = ["111", "200", "220"]


def main() -> None:
    df = pd.read_csv(INPUT)
    total = len(df)
    assigned = int(df["nearest_pt_fcc_reflection"].isin(TARGETS).sum())
    rows = []
    for face in TARGETS:
        count = int(df["nearest_pt_fcc_reflection"].eq(face).sum())
        rows.append({
            "reflection": face,
            "reference_d_spacing_nm": {"111": 0.226, "200": 0.196, "220": 0.139}[face],
            "candidate_count": count,
            "proportion_of_all_fft_candidates": count / total if total else 0.0,
            "proportion_of_assigned_pt_fcc_candidates": count / assigned if assigned else 0.0,
            "interpretation": "FFT reflection candidate, not confirmed exposed facet",
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "pt_fcc_reflection_candidate_proportions.csv", index=False)
    (OUT / "pt_fcc_reflection_candidate_proportions.txt").write_text(
        f"Pt fcc FFT reflection candidate proportions\n\n"
        f"Total FFT candidates: {total}\n"
        f"Assigned Pt fcc candidates: {assigned}\n"
        f"Unassigned candidates: {total - assigned}\n\n"
        + summary.to_string(index=False)
        + "\n\nThese are FFT d-spacing candidate proportions. They are not exposed-facet area fractions.\n",
        encoding="utf-8",
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
