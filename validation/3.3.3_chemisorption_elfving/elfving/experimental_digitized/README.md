# Pixel-digitized experimental controls

These SI-unit CSV files are the precomputed pixel digitizations recovered from the red experimental traces in the plot-only crops under `paper_figures/`. Normal notebook execution only loads these cached tables; it does not rerun image segmentation or axis calibration. Run `python3 digitize_breakthrough.py` only when the paper crops or calibration deliberately change.

The crop boundaries are calibrated to the published axes:

- Figures 5a-7a: 0-7200 s and dimensionless CO2 mole fraction 0-0.0004.
- Figure 8a inset: 0-27000 s and dimensionless CO2 mole fraction 0-0.0004.
- Figures 5b/6b (H2O) and 5c/6c (temperature): 0-7200 s.
- Figures 8b/8c (H2O, temperature): 0-162000 s, matching the paper's full-duration axis for that run rather than its own 0-7200 s inset, since the notebook's temperature/loading plots need values out to RUPTURA's 57600 s output horizon.

Red pixels are separated from the black paper simulations by color thresholding. At each calibrated time, vertical red-pixel groups are clustered into the expected number of experimental traces and lightly median-filtered. Thus, these values are digitized graphical controls, not the authors' original numerical samples. Pixel thickness, anti-aliasing, overlapping curves, and the 150 dpi source rendering limit their precision. `figure_5b.csv`'s 1 vol-% and 2 vol-% series carry one uncorrected ~6-sample glitch around t=6000-6400 s where a marker glyph momentarily confused the column clustering; two general fixes were tried (a wide-window despike pass, seeding each column's clustering from the previous one) and both made other panels worse, so this one narrow artifact was left in rather than risk that.

CO2 columns are named `co2_mole_fraction_<humidity>` (mole fraction, e.g. 0.00037 at 370 ppm). H2O columns are `h2o_vol_percent_<humidity>` (vol-%, as printed on the paper axes). Temperature columns are `temperature_C_<humidity>` (degC, as printed on the paper axes); the notebook converts to kelvin on load. Humidity suffixes retain the paper identifiers (`0p2`, `1`, and `2`, or `experimental` for Figure 8's single run). The corresponding measured mole fractions used as simulation inputs are 0.0016, 0.0095, and 0.0190.

The notebook reports 0.10 and 0.50 breakthrough times from the CO2 controls. A 0.01 crossing is only about 4e-6 mole fraction and is too close to the rendered baseline and line thickness to recover robustly at this resolution.
