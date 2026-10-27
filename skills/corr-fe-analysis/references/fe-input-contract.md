# Product-amount input contract

`scripts/calculate_fe.py` accepts a UTF-8 CSV with `product,moles,electrons,charge_C`. Use one calibrated charge interval per file; every row must carry the same *positive magnitude* of integrated charge in coulombs. `moles` is the product amount made in that interval, not an unreacted CO feed amount, raw GC composition, or uncorrected HPLC concentration. Use a balanced electron count for each product.

Run `python scripts/calculate_fe.py --input <reviewed.csv> --output <fe.csv>`. Audit gas flow, liquid volume, blank, calibration, injection delay, and charge alignment upstream. Calculate carbon balance separately, with CO feed and crossover stated. Do not include unreacted CO in the product FE sum.
