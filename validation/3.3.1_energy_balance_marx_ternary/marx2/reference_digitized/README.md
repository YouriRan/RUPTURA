# Digitized experimental reference data

These CSV files were extracted from the vector marker paths in the Marx et al.
paper and supporting-information PDFs. Each trace contains exactly 50 points.
The breakthrough values are mole fractions and temperatures are converted from
degrees Celsius in the source figures to kelvin. `Temp` is the outlet-side T5
sensor, matching the Ruptura outlet comparison. Every run directory includes
provenance metadata with the source PDF, page, figure, panel, and raw marker
count. Run `digitize_paper.py` to reproduce the extraction.
