# 3.2.2 Macrostate particle distribution (MPD)

Source: `sharma/` (moved here on 2026-11-19; `../../sharma` is a symlink to it) — CO₂/propane in all-silica MOR at
300 K (Sharma et al., Mol. Phys. 121 (2023) e2183721). The binary reference distribution is the N-d TMMC ln Π
`sharma/lnpi_co2_propane_mor_300K_nd_tmmc_walks.txt` (computed on another machine, see its header), written for RUPTURA
by `sharma/mpd/build_nd_tmmc_distribution.py`. The one-component distributions of the earlier RASPA3 TMMC windows
(`sharma/tmmc/`, slices with an edge) serve the temperature reweighting at 300 and 350 K. Earlier MPD examples:
`../../examples/MPD/BEA-alkanes/` (multinomial distribution, reproduces competitive Langmuir) and `../../examples/MPD/`
(CO₂/N₂).

Case `mpd-mor` in `../cases.csv`: `sharma/runs/generate.py` writes the 16 runs, `sharma/runs/run.py` runs them
(`RUPTURA_BIN`), and `analysis.ipynb` (run from this directory) writes `figures/mpd_mor/` of the paper and prints the
numbers of its text (Results subsection "MPD", `sec:results_mpd`).

Status (2026-11-19): in the paper. All runs with `build/ruptura-sharma-2026-11-19-mpd-jac` (two-pass MPD reweighting
and MPD sensitivities, TODO P-PHY-9 and P-PHY-10); the wall times were taken with other runs on the machine.

## Points

- [x] Pure-component MPD against the Table 3 dual-site Langmuir fits: from the binary distribution, −7 to +21 % (CO₂,
      up to 1e6 Pa) and −8 to +17 % (propane, up to 6.7e5 Pa); MPD curves are pale beyond the trusted range (edge of
      the sampled domain above 1e-3 of the probability).
- [x] Equimolar mixture, 10²–10⁸ Pa: MPD follows SIAST (−10 to +17 % CO₂, −4 to +36 % propane up to 1.2e7 Pa),
      IAST overpredicts the selectivity by up to two orders of magnitude. GCMC points of Fig. 6 not digitized.
- [x] Breakthrough of Fig. 10a with MPD, SIAST, and IAST (ROS34PW2): MPD t50 within 1–2 % of SIAST; propane roll-up
      1.05, 1.14, 1.71. MPD 160 s against 6 s for SIAST.
- [ ] Temperature reweighting (B-PHY-1 in `../../TODO.md`): done for the pure components from the one-component
      windows (2026-11-11); the binary N-d distribution exists at 300 K only.
- [x] Paper: the distribution is reported as two-dimensional NVT+W (Chen and Lin 2023, `chenlin2023`) while the
      N-d TMMC is work in progress.
