# Particle and FFT candidate QC

The optional Pt/C profile in `scripts/run_tem_analysis.py --pt-candidates` combines the former image-inventory and particle/facet workflows in one route. It is intended for dark Pt-containing features against a lighter support, not for arbitrary bright particles or other materials.

## Calibration and grouping

The script converts recognized DM pixel units to nm/pixel. Images from `0.02` through `0.08 nm/pixel` are eligible for particle-size candidates. Finer images are reserved for lattice review; coarser images are overview only. Inspect the exported image-group table before accepting these boundaries for a new instrument or specimen.

## Candidate segmentation

Intensity is normalized between its 1st and 99th percentiles. A broad local background minus the denoised image yields a dark score. Strongly dark cores seed lower-threshold hysteresis growth; morphological cleanup and distance-transform watershed separate nearby objects. The automatic review subset uses equivalent diameter `0.8–10 nm` and roundness at least `0.25`. Those filters prioritize inspection only: support texture, overlap, and edge artifacts can still pass.

Keep the pending-QC table, overlay PNGs, rejected candidates, and per-image summaries. A publication-ready distribution requires manual mask review across comparable fields. Do not turn the raw automated mean into a final particle size.

## Lattice and facet limits

The Pt profile checks candidate FFT spacings near `{111}=0.226 nm`, `{200}=0.196 nm`, and `{220}=0.139 nm`, and compares a supplied Pt powder-card peak list within `±0.015 nm`. Verify that the reference card and tolerance are appropriate for the specimen. Report both the fraction among matched candidates and among *all* FFT candidates, including unmatched peaks. Powder-card relative intensities, FFT candidate counts, exposed projected facet length, and 3D surface-area fractions are different quantities.

For a facet claim, manually document local d-spacing, zone axis, particle boundary direction, traced segment length, assignment, and confidence. Keep the automated output labelled as candidate reflection evidence until those checks are complete.
