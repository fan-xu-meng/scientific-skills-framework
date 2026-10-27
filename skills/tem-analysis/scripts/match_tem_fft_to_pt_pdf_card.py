#!/usr/bin/env python3
"""Match TEM FFT d-spacing candidates to the supplied Pt PDF card."""

from pathlib import Path
import os
import pandas as pd

OUT = Path(os.environ["TEM_OUTPUT_DIR"]).resolve() / "pt-candidates"
FFT = OUT / "lattice_spacing_candidates_pending_indexing.csv"

PDF_PEAKS = [
    {"hkl": "111", "two_theta_deg": 40.04, "d_angstrom": 2.25, "intensity_relative": 100},
    {"hkl": "200", "two_theta_deg": 46.53, "d_angstrom": 1.95, "intensity_relative": 30},
    {"hkl": "220", "two_theta_deg": 67.86, "d_angstrom": 1.38, "intensity_relative": 16},
    {"hkl": "311", "two_theta_deg": 81.50, "d_angstrom": 1.18, "intensity_relative": 16},
    {"hkl": "222", "two_theta_deg": 85.95, "d_angstrom": 1.13, "intensity_relative": 3},
    {"hkl": "331", "two_theta_deg": 117.71, "d_angstrom": 0.90, "intensity_relative": 3},
    {"hkl": "420", "two_theta_deg": 124.60, "d_angstrom": 0.87, "intensity_relative": 2},
]
TOLERANCE_NM = 0.015


def main() -> None:
    fft = pd.read_csv(FFT)
    peak_rows = []
    match_rows = []
    for peak in PDF_PEAKS:
        ref_nm = peak["d_angstrom"] / 10.0
        diff = (fft["d_spacing_nm"] - ref_nm).abs()
        matched = fft[diff <= TOLERANCE_NM].copy()
        row = {**peak, "d_nm": ref_nm, "match_tolerance_nm": TOLERANCE_NM,
               "matched_fft_candidate_count": int(len(matched)),
               "matched_file_count": int(matched["file_name"].nunique()) if len(matched) else 0}
        peak_rows.append(row)
        for idx, candidate in matched.iterrows():
            match_rows.append({**peak, "pdf_d_nm": ref_nm, "match_delta_nm": float(diff.loc[idx]),
                               "fft_file": candidate["file"], "fft_file_name": candidate["file_name"],
                               "fft_candidate_rank": int(candidate["candidate_rank"]),
                               "fft_d_spacing_nm": float(candidate["d_spacing_nm"]),
                               "fft_relative_strength": float(candidate["fft_relative_strength"])})
    pd.DataFrame(peak_rows).to_csv(OUT / "pt_pdf_card_all_peaks_summary.csv", index=False)
    match_columns = ["hkl", "two_theta_deg", "d_angstrom", "intensity_relative",
                     "pdf_d_nm", "match_delta_nm", "fft_file", "fft_file_name",
                     "fft_candidate_rank", "fft_d_spacing_nm", "fft_relative_strength"]
    pd.DataFrame(match_rows, columns=match_columns).to_csv(
        OUT / "pt_pdf_card_to_tem_fft_matches.csv", index=False)
    lines = ["Pt PDF card to TEM FFT matching report", "", "Source: PDF Card - 00-001-1194, Pt", "", ""
             "The PDF card lists seven Pt fcc reflections. TEM FFT candidates were matched within +/-0.015 nm.", ""]
    for row in peak_rows:
        lines.append(f"{{{row['hkl']}}}: d={row['d_nm']:.3f} nm, I={row['intensity_relative']}, "
                     f"2theta={row['two_theta_deg']:.2f} deg, matches={row['matched_fft_candidate_count']} "
                     f"across {row['matched_file_count']} files")
    lines += ["", "A match is a d-spacing candidate only. It does not prove an exposed surface facet or zone axis."]
    (OUT / "pt_pdf_card_to_tem_fft_matches.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
