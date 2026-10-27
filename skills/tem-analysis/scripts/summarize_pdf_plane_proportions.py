#!/usr/bin/env python3
from pathlib import Path
import os
import pandas as pd

OUT = Path(os.environ["TEM_OUTPUT_DIR"]).resolve() / "pt-candidates"


def main() -> None:
    peaks = pd.read_csv(OUT / "pt_pdf_card_all_peaks_summary.csv")
    total_fft = int(pd.read_csv(OUT / "lattice_spacing_candidates_pending_indexing.csv").shape[0])
    total_matched = int(peaks["matched_fft_candidate_count"].sum())
    peaks["proportion_of_matched_pdf_plane_candidates"] = peaks["matched_fft_candidate_count"] / total_matched if total_matched else 0.0
    peaks["proportion_of_all_fft_candidates"] = peaks["matched_fft_candidate_count"] / total_fft if total_fft else 0.0
    peaks["unmatched_fft_candidate_count"] = total_fft - total_matched
    peaks.to_csv(OUT / "pt_pdf_card_plane_proportions.csv", index=False)
    print(peaks[["hkl", "matched_fft_candidate_count", "proportion_of_matched_pdf_plane_candidates", "proportion_of_all_fft_candidates"]].to_string(index=False))


if __name__ == "__main__":
    main()
