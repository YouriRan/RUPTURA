# Paper-panel digitization sources

These cropped plot regions come from the user-supplied Elfving and Sainio (2021) PDF. `digitize_breakthrough.py` reads their red pixels and writes calibrated experimental controls under `experimental_digitized/`.

- `figure_5a.png`: Figure 5a, CO2 breakthrough, 0-120 min and 0-400 ppm.
- `figure_6a.png`: Figure 6a, CO2 breakthrough, 0-120 min and 0-400 ppm.
- `figure_7a.png`: Figure 7a, CO2 breakthrough, 0-120 min and 0-400 ppm.
- `figure_8a_zoom.png`: Figure 8a inset, CO2 breakthrough, 0-450 min and 0-400 ppm.
- `figure_5b.png` / `figure_6b.png`: panel b, H2O breakthrough, 0-120 min and 0-2 vol-%.
- `figure_5c.png` / `figure_6c.png`: panel c, bed temperature, 0-120 min; 24-36 degC (fig. 5) or 34-44 degC (fig. 6).
- `figure_8b.png`: Figure 8b, H2O breakthrough, 0-2700 min and 0-2 vol-%.
- `figure_8c.png`: Figure 8c, bed temperature, 0-2700 min and 24-30 degC.

Figure 7 has no H2O or temperature panel in the paper, so there is no `figure_7b.png`/`figure_7c.png`. Each of the panels above also has its own small inset (a zoom of the first ~20 min, or ~500 min for fig. 8c); `digitize_breakthrough.py` masks those out (`mask_boxes`) along with the legend before digitizing, since they duplicate the same red trace at a different pixel scale and would otherwise corrupt the calibration.

The crops are source inputs and are not used as notebook backgrounds. Red traces/markers are experimental results in the paper; black traces are the paper simulations and are rejected by the color threshold.
