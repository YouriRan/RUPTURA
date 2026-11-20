# MPD validation: CO₂/propane in MOR (Sharma et al. 2023)

Validation of the macrostate particle distribution (MPD) mixture prediction against SIAST and IAST on the
mordenite case of Sharma et al., Mol. Phys. 121 (2023) e2183721 (V-PHY-7 in `../TODO.md`, paper Results
subsection "MPD", `../validation/3.2.2_mpd/`). Unit conversions are in `conversions.md`.

## Pipeline

1. `tmmc/` — RASPA3 transition-matrix Monte Carlo of the binary in 2 × 2 × 4 unit cells of all-silica MOR at
   300 K and 350 K. `generate.py` writes 422 windows per temperature: one TMMC over the propane count (the edge)
   and, for every fixed propane count m = 0…100, one TMMC over the CO₂ count n = 0…200 − m (the slices). Each
   window also writes the potential-energy cumulants per macrostate (`tmmc/tmmc_energies.txt`). `run.py [300 350]`
   runs them on 8 workers. The executable is a frozen copy of the continuous-lambda RASPA3 build
   (`~/sandbox/continuous-lambda`, slices with a fixed-population co-component and the energy cumulants, C8 in
   that project's TODO).
2. `mpd/build_distribution.py [T]` — stitches the windows into ln Π(n, m) and ⟨H⟩(n, m) and writes RUPTURA's
   distribution file `particle_distribution_T<T>K.data` with its `MPDSettings` (`.json`) and `lnpi_T<T>K.npz`.
   At 300 K the slices did not converge (Status); `mpd/build_nd_tmmc_distribution.py` writes the same three files from
   the N-d TMMC ln Π `lnpi_co2_propane_mor_300K_nd_tmmc_walks.txt`, computed on another machine (see its header).
3. `runs/generate.py`, `runs/run.py` — RUPTURA mixture predictions (equimolar 10²–10⁸ Pa with MPD, SIAST,
   IAST; pure components with MPD; 350 K from the 300 K and from the 350 K distribution) and the breakthrough of
   Fig. 10a (10 % CO₂, 10 % propane in helium, 10⁵ Pa) with MPD, SIAST, and IAST.
4. `analysis.ipynb` — figures.

## Force field

As in Sec. 3 of the paper: TraPPE-zeo Lennard-Jones for the framework, Garcia-Sanchez CO₂ (rigid), united-atom
propane with the Dubbeldam LJ parameters and TraPPE bonds and bend (CBMC), Lorentz–Berthelot mixing, LJ truncated
and shifted at 12 Å, Ewald. The framework charges are those of the RASPA2 CIF (Si +2.05, O −1.025), not the
TraPPE-zeo charges (+1.5/−0.75): at 10⁴ Pa and 300 K, GCMC (RASPA3) gives

| | CIF charges | TraPPE-zeo charges | Table 3 fit |
| --- | --- | --- | --- |
| CO₂ | 1.11 ± 0.07 mol/kg | 0.39 ± 0.08 mol/kg | 1.18 mol/kg |
| propane | 0.77 ± 0.05 mol/kg | (no charges) | 0.77 mol/kg |

RASPA2 uses the CIF charges by default, so the paper presumably did too. The RASPA3 build of March 2025
(`~/micromamba/bin/raspa3`) ignores the bonds and bends of `propane.json` and places all three beads at one point;
use a build from August 2026 or later.

## Status

- 2026-11-03 16:37: TMMC at 300 K started (`tmmc/run_300K.log`; 0.55–0.7 ms per MC step, about 14 h on
  8 workers). Queued behind it (detached shells waiting on the batch's PID):
  - the 350 K TMMC (`tmmc/run_350K.log`), about 14 h more;
  - at the same time, `mpd/build_distribution.py 300` (`mpd/build_distribution_300K.log`), `runs/generate.py`, and
    all RUPTURA runs (`runs/run.log`) with the frozen binary `../build/ruptura-sharma-2026-11-03` (the
    `build/ruptura` of 15:03; other sessions rebuild `build/ruptura`). The wall times of these runs are taken under
    the load of the 350 K batch;
  - after 350 K: `mpd/build_distribution.py 350` and the run `mixture_equimolar_MPD_T350K_ref350K`.
  Then run `analysis.ipynb`.
- 2026-11-06: both batches finished (2.4× the planned steps: RASPA3 runs max(N, 20) steps per component per
  cycle), but most windows did not converge (TODO V-PHY-7): `mpd/build_distribution.py` stops at "missing
  states" at both temperatures and nothing downstream ran. Moved to `tmmc/unconverged_2026-11-03/`.
- 2026-11-06 16:44: rerun started with `MacroStateRezeroAfterInitialization` (new key in the sandbox, its C12:
  the initialization statistics are discarded), 10⁵ production steps per macrostate (counted per component now),
  equilibration 0.5 and initialization 0.1 of production; 2.6 × 10⁹ steps, about 45 h per temperature. One
  detached chain (`sharma_chain.sh` in the session scratchpad): 300 K TMMC; then the 350 K TMMC and, at the same
  time, the 300 K distribution and all RUPTURA runs; then the 350 K distribution and its reference run. Logs as
  above. A short test of two windows: with the reset, CO₂ N 40–47 at m = 0 has the deletion acceptance of the
  window below (0.003); the propane edge above ~50 molecules (insertion acceptance 10⁻⁵–10⁻⁴) stays at its top
  states in production and will probably not converge at any run length. Missing production energies are now a
  warning in `build_distribution.py` (⟨H⟩ = 0 there), missing ln Π still an error.
- 2026-11-09: 300 K finished on 11-08 14:08; `mpd/build_distribution.py 300` stops at slice m = 20, N 40–82
  (states 40–63 never visited after the reset), so no RUPTURA run started. The 350 K TMMC continues (252/422 windows
  on 11-09 15:50, end about 11-10 midday; its edge and m = 0 windows are among the last). Per window (unvisited
  states, or a visited state whose insertion or whose upper neighbour's deletion recorded zero acceptance):
  - the edge and the m = 0 slice (8 windows) have no such defect, and ln Π of neighbouring windows agrees to ≤ 0.7,
    but they are not flat above 160 CO₂ and 50 propane molecules (top CO₂ window: 98 % of the states below 1 % of the
    most visited one; propane windows from 40: 49 % and 62 %). CO₂ has one bad state at N = 36 (deletion acceptance
    1.8 × 10⁻⁴ against 3.2 × 10⁻³ at N = 35), which gives a step in the isotherm near 0.8 mol/kg. Pure isotherms
    against Table 3: −4…+21 % up to 160 CO₂ (2 × 10⁶ Pa), −7…+18 % up to 50 propane (6 × 10⁶ Pa).
  - 287 of the 414 slices with propane fail the test (49 with unvisited states; the others have ln Π jumps of
    tens to thousands across the zero-acceptance links), at every m and every N range, including the low loadings.
    The walks sit at N = 0: in the first window of a slice the median state is visited ~10³ times against ~10⁶ at
    N = 0 (m = 0: 1.3 × 10⁵ median). The batch of 11-03 fails the same test in 236 of these 414 windows; the
    reset and the longer runs made the pure windows usable, not the slices. With at least 1000 visits per state, the
    first unusable link of a slice lies at a median of n = 14 (quartiles 6 and 18); no slice with propane reaches
    n = 42, while SIAST needs n = 42, m = 27 at 10⁴ Pa and the breakthrough feed n = 54, m = 35.
- 2026-11-09, what exists without the binary distribution (`../validation/3.2.2_mpd/analysis.ipynb`, paper Results
  subsection "MPD", marked as draft):
  - `mpd/build_pure_distribution.py`: one-component distributions from the pure windows, runs
    `runs/pure_*_MPD_T*K_ref*K`.
  - `mpd/build_low_loading_distribution.py`: binary distribution from the states below the first unusable link of
    each slice, run `runs/mixture_equimolar_MPD_low_loading`, valid up to 400 Pa (< 1 % change between 300 and 10000
    minimum visits; 15 % at 500 Pa); it does not separate SIAST and IAST (they differ by ≤ 6 % there).
  - `runs/{mixture_equimolar,breakthrough}_{SIAST,IAST}` ran under the load of the 350 K TMMC (no wall times).
- 2026-11-11: the 350 K TMMC finished on 11-10 16:51; `mpd/build_distribution.py 350` stops like 300 K (slice m = 44,
  N 40–82). Pure windows: CO₂ flat up to 162; propane 0–42 has no recorded insertion at 32 → 33 (acceptance 9e-12;
  its lower states were visited in equilibration only), so pure propane at 350 K holds up to 32 molecules.
  `runs/pure_*_MPD_T350K_ref{300,350}K` compare the reweighting (B-PHY-1): for CO₂, 300 K → 350 K with
  N ln(T_ref/T) is within −5…+0.5 % of the 350 K distribution at 10²–10⁶ Pa, without it +11…+16 % up to 10⁴ Pa.
  `mpd/build_pure_distribution.py` now fills ⟨H⟩ of states without production energies by linear interpolation in
  N (linear fit to the nearest five beyond the sampled range); `build_distribution.py` still writes 0 there.
- 2026-11-19: binary 300 K distribution from the N-d TMMC (`lnpi_co2_propane_mor_300K_nd_tmmc_walks.txt`, 9865
  states, all visited) through `mpd/build_nd_tmmc_distribution.py` (`conversions.md`: reference fugacity 10⁵ Pa kept,
  ⟨H⟩ with the ideal-gas intramolecular energy of propane, as the file header prescribes). `runs/generate.py` now also
  writes the 300 K MPD runs. At the breakthrough feed the edge of the sampled domain carries 3 × 10⁻⁵¹ of the
  probability (header limit 10⁻³); RUPTURA's MPD loading at the saturated inlet, 1.1529 and 0.8141 mol/kg, equals the
  direct reweighting of the file. The breakthroughs now use ROS34PW2 (`runs/generate.py`): SIAST and IAST with
  `../build/ruptura-sharma-2026-11-19-mpd` (7.4 s and 4.3 s; with CVODE 69 s and 41 s, same curves to 0.1 s in t50,
  and the same as the 11-03 binary), MPD with `../build/ruptura-sharma-2026-11-19-mpd-jac` (160 s, 20980 steps; `.txt`
  files next to the binaries). On the way, MPD was made faster in RUPTURA (TODO P-PHY-9: 3.2× per evaluation) and got
  equilibrium sensitivities for the Jacobian (P-PHY-10: without them ROS34PW2 needed 87 min for 1435 s; CVODE about
  40 min). The MPD run agrees with the earlier CVODE and sensitivity-free runs to ≤ 1e-3 of the feed. Result:

  | | t50 CO₂ | t50 propane | propane roll-up |
  | --- | --- | --- | --- |
  | MPD | 3344 s | 2291 s | 1.053 |
  | SIAST | 3380 s | 2240 s | 1.138 |
  | IAST | 2473 s | 1527 s | 1.705 |

  Not run: `mixture_equimolar_MPD`, `mixture_pure_*_MPD`, `mixture_equimolar_MPD_T350K_ref300K` (their
  10²–10⁸ Pa range extends beyond the domain, valid up to about 10⁶ Pa total fugacity).
- The pipeline was tested end to end on a 12-molecule grid (scratch, 28 min): stitching, RUPTURA MPD/SIAST/IAST
  mixture and breakthrough runs, and the notebook. Its Henry coefficients from Π(1,0)/Π(0,0) and Π(0,1)/Π(0,0) were
  16 % below (CO₂) and 14 % above (propane) the Table 3 fits; check them on the full run.
