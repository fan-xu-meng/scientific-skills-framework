# Synthetic TEM QC example

This is invented data, not a measurement or a validated Pt/C result.

| image | pixel size | automated candidates | recommended for manual review | interpretation |
| --- | ---: | ---: | ---: | --- |
| overview-A.dm4 | 0.050 nm/pixel | 41 | 29 | inspect 29 masks before size reporting |
| lattice-B.dm4 | 0.010 nm/pixel | 0 | 0 | use for lattice/FFT review, not size pooling |

A hypothetical `0.224 nm` FFT spacing differs by `0.002 nm` from the `{111}` reference used by the optional Pt profile. It is a possible `{111}` reflection. Without an indexed zone axis and traced particle boundary, the correct report is “candidate `{111}` reflection,” not “the particle exposes a `{111}` facet.”
