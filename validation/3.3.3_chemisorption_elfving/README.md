# 3.3.3 Chemisorption: Elfving and Sainio (2021)

Reference: J. Elfving, T. Sainio, Kinetic approach to modelling CO₂ adsorption from humid air using
amine-functionalized resin: Equilibrium isotherms and column dynamics, Chem. Eng. Sci. 246 (2021)
116885, doi:10.1016/j.ces.2021.116885.
Source: `elfving/` (13 cases, `generate.py`, `elfving_numbered_examples.ipynb`,
`experimental_digitized/`).
Status: runs exist; paper figures missing (`elfving/paper_figures/` is empty).

## Points

- [ ] Paper figures for the reproductions of Figs. 5–8 (breakthrough, midpoint solid temperature, loading).
- [ ] Run 1 (`five_25c_0p2pct_500mlmin`) never reaches CO₂ breakthrough in the simulated time; extend
      the duration.
- [ ] Disclose the mapping approximations:
      - the two-pathway amine model is mapped to one `General` site with capacity $q_m/2$ and a squared
        free-site term;
      - the H₂O GAB and Langmuir terms are `FirstOrder` sites;
      - the H₂O LDF coefficient is capped at 0.25 s⁻¹.
- [ ] Elfving does not use film/pore transport (`usePoreSurfaceTransport` is off). A pore-transport case
      is still to be found (`../../TODO.md`); until then only 3.1 covers it.
- [ ] Elovich and PseudoNth ($n \neq 1$) are covered only by 3.1.
- Serna-Guerrero is dropped (the reference is inconsistent). Remove its SI section and the unused
  `figures/serna-guerrero-fig7.pdf` from the paper.
