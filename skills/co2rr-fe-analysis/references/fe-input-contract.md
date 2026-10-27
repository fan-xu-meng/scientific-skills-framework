# Product-amount input contract

`scripts/calculate_fe.py` accepts a UTF-8 CSV with `product,moles,electrons,charge_C`. Use one calibrated charge interval per file; every row must carry the same *positive magnitude* of integrated charge in coulombs. `moles` is the amount made within that interval, not a GC mole fraction or an uncorrected HPLC concentration. `electrons` is the balanced reaction's electron count per product molecule.

Run `python scripts/calculate_fe.py --input <reviewed.csv> --output <fe.csv>`. Verify gas flow, liquid volume, blank, calibration, residence-time correction, and interval alignment *before* this arithmetic check. Compute carbon balance separately from carbon atoms per product and a declared feed/input basis. Do not equate a measured FE sum with carbon conversion.
