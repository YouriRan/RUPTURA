# Validation cases for the RUPTURA 2.0 paper

One directory per results subsection (Section 3 of `ruptura-2.0-paper/main.tex`; the numbers in the names predate
its current order). Each README lists the points still open for that subsection. Features, open decisions, and paper
fixes are collected in `../TODO.md`.

The case directories live here, so that `validation/` holds every case behind the paper. The scripts and notebooks
find the repository from their own location.

| Section | Directory | Case | Source | Status |
|---|---|---|---|---|
| 3.1 | `3.1_verification` | Exact solutions and limits | unit tests, `benchmark/multibed`, new cases | partly done |
| 3.2.1 | `3.2.1_mixture_methods_temperature` | SPI/SCI/SEI/IAST/SIAST; temperature dependence | 3.5 (`figures/`), 3.3.4, 3.3.1 | new figure |
| 3.2.2 | `3.2.2_mpd` | MPD vs SIAST and IAST: CO₂/propane in MOR (Sharma et al. 2023); earlier BEA nC7/2mC6 | `sharma/` in it | in paper |
| 3.2.3 | `3.2.3_gab_water` | GAB water isotherm | 3.3.3, 3.3.4 | runs exist |
| 3.3.1 | `3.3.1_energy_balance_marx_ternary` | Marx et al. 2014, T1–T24 | `marx2/` in it | in paper |
| 3.3.2 | `3.3.2_monoliths_stampi_bombelli` | Stampi-Bombelli et al. 2024 | `stampi-bombelli/` in it | in paper |
| 3.3.3 | `3.3.3_chemisorption_elfving` | Elfving and Sainio 2021 | `elfving/` in it | runs exist, no figures |
| 3.3.4 | `3.3.4_humid_dac_liu` | Liu et al. 2024 | `liu/` in it | runs done, figures partly |
| 3.3.5 | `3.3.5_reactions_stadler` | Stadler et al. 2022, SE-WGS | `stadler/` in it | in paper (chemisorption, reactions; SE-WGS in SI) |
| 3.3.6 | `3.3.6_liquid_ph` | to decide with Shrinjay | `shrinjay`, `examples/Liquid` | open |
| 3.3.7 | `3.3.7_multistack_park` | Park et al. 1998, layered beds | — | to set up |
| 3.3.8 | `3.3.8_pore_transport` | Pore transport with `General` chemisorption: CO₂ from air on PEI/fumed silica, four particle sizes (Goeppert et al. 2014, model of Auton et al. 2024) | this directory | rerun 2026-11-20 (byte-identical to 2026-11-09); in the paper since 2026-11-20 (`sec:results_chemisorption`, SI `apx:goeppert`) |
| 3.4.1 | `3.4.1_boundary_conditions_marx_10bar`, `3.4.1_boundary_conditions_bea_dilute` | Five boundary modes on the Marx 1e6 Pa runs and a dilute BEA control | these directories | in paper; ROS34PW2 and Danckwerts since 2026-11-03 (all runs complete) |
| 3.4.2 | `3.4.2_swing_adsorption` | TSA cycles on the Marx T15 column; PSA deferred | 3.3.1 | TSA in paper |
| 3.4.3 | `3.4.3_sensitivity` | One-at-a-time sensitivity of Marx T15 | 3.3.1, `benchmark` levels | in paper |
| 3.4.4 | `3.4.4_flow_reversal` | Co-current and counter-current regeneration of a trace Langmuir adsorbate against equilibrium theory (`FlowDirection`) | this directory | in paper (2026-11-06) |
| 3.5 | `3.5_numerical_performance` | Integrator accuracy and cost over complete breakthroughs (BEA C7, Marx T15, Stampi-Bombelli monolith with fitted D_ax), and the convection schemes | `breakthrough/`, `convection/`, and `figures/` (helpers of `breakthrough/generate.py`) in it | in paper; SIRK3 later |
| 3.x | `self_driving_lab` | Bayesian optimization of an H₂ PSA cycle on one column and on two beds | `Self-driving-lab-H2-PSA/` in it | in paper |
| 3.x | `inlet_condition_dirichlet_vs_danckwerts` | Every case above with physical dispersion, rerun with the Dirichlet and the Danckwerts inlet (`InletCondition`) | 3.3.1, 3.3.2, 3.3.5, 3.4.1–3.4.3 | done 2026-11-03; Danckwerts is the default since then, paper studies to rerun (TODO I-BC-4) |
| 3.x | `ros34pw2_danckwerts_rerun` | Every breakthrough and swing case of 3.3–3.4 rerun with ROS34PW2 and the Danckwerts inlet, `-O3` binary | 3.3.1, 3.3.2, 3.3.4, 3.3.5, 3.4.1–3.4.3 | 223 of 223 done 2026-11-03; Elfving is left out until TODO B-EX-9 is decided (B-EX-16) |

Dropped: Serna-Guerrero (the reference is inconsistent).
`IAST-Landau2012/` and `IAST-Moon-Tien-casestudy/` are RUPTURA 1 leftovers (run scripts only).

## Conventions

- SI units throughout; a `conversions.md` for every case recreated from a paper.
- Simulated cases have a rerunnable `generate.py` and a `run.py` that uses 8 workers.
- Digitized experimental data go in `reference_digitized/`.
- Paper figures are 8×6 PDFs. Reference data are drawn as dashed lines with open markers, RUPTURA as
  solid lines.
- Runs whose analysis reads only the outlet set `"WriteOutletOnly": true` in their generator (since 2026-11-20): the
  3.2.2 MPD breakthroughs, Liu (all but Figs. 4 and 5a), the Stadler activity tests and fresh bed, 3.3.8 (all but the
  wave runs), 3.4.1 BEA dilute, 3.4.3, both 3.5 studies but `figures/`, the T15 convection-scheme check, the
  phase-end and KLU runs of the free-pressure check, and the ROS34PW2 cross-check. Their readers detect the one-row
  blocks from the data, so old full outputs still read. Cases that use profiles, inventories, interior temperatures,
  or the z = 0 row keep the full output.

## Full rerun on a cluster node

Every RUPTURA run behind the paper can be redone in one go on one machine, so that all timings come from the same
hardware and build. The cases are in `validation/` (above). Their generators
also read base inputs from the tracked `examples/`, and the runs need the source to build the binary. `validation/`
is in git since 2026-11-20, without the RUPTURA output (`.gitignore` names it). The paper, whose figure folders the
notebooks write to, is a separate repository, so copy the working tree, not a git checkout.

### The job list

- `cases.csv` is the central list of cases, in the order in which their inputs have to be generated. Its columns are:
  - `paper_section`: the `\label` of the paper section, or `cross-check` / `not yet`;
  - `directory`: the case directory, relative to the repository root;
  - `generate`: the command that rewrites the inputs, run in `directory`; empty where the inputs are kept as they are;
  - `simulations`: the paths (globs) below `directory` that hold its simulation directories;
  - `timeout_s`, `timed`, `in_paper`, `analysis`, and `notes`.

  The last three rows (the self-driving lab, the unit tests, and the VSA in progress) have no simulation directories.
- `python run.py --paper --list` (the `run.py` at the repository root) writes `jobs.csv`, one row per simulation
  directory, in the order of `cases.csv`. Its columns are:
  - `job`, `paper_section`, `case`, `timed`, `timeout_s`, and `in_paper`;
  - `directory`, relative to the repository root;
  - `previous_wall_time_s`, from the run.log or ruptura.log found there.

  `--generate` first reruns every generator in the order of `cases.csv`, keeping any output already there. Without
  `--paper`, the list holds every simulation directory of the repository.
- On 2026-11-19 the list held 1078 jobs. Of these, 645 are not in the paper (`in_paper` = `no`):
  - 532 cross-checks, which compare every case under the Dirichlet and Danckwerts inlets and under ROS34PW2;
  - the 44 runs of the pore-transport case;
  - the Elfving (13) and Liu (56) runs.

  Since 2026-11-20 the ROS34PW2 cross-check has no Elfving copies (TODO B-EX-16), which takes 13 jobs off the
  cross-checks (519), and its runs have a timeout of 3600 s.
- `python run.py --paper` runs the same cases on this machine: the pool first, with the case timeouts, then the
  timed studies through their own `run.py`. It gives them as many workers (`-j`) as the pool.

A job is one RUPTURA process in its directory: it reads `simulation.json` there and writes its output files next to
it. For one job, `timeout 0` means no limit:

```sh
cd "$REPO/$directory" && timeout "${timeout_s:-0}" "$RUPTURA_BIN" > run.log 2>&1
```

### Steps

1. **Pack** (on the Mac, once no run is writing): `python validation/pack.py` writes
   `ruptura-rerun-<date>.tar.gz` next to the repository; `-o PATH` names another file, and `--list` only prints the
   files it would pack. The script docstring gives the full rules:
   - **Packed:** the working tree, with the paper, every input, and the input `.data` files such as the MPD
     distributions.
   - **Left out:** every RUPTURA output (by file name, including the run logs), `.git`, the builds, the Ruptura Lab
     records in `simulations/`, and the RASPA3 TMMC data of `3.2.2_mpd/sharma/tmmc` (11 GB, not a RUPTURA job).

   On 2026-11-19 it packed 11143 files, 265 MB, and 146 MB compressed, in 12 s. All 1078 jobs and every file that
   their inputs name were in it. Unpacked, `run.py --paper --list` gave the same 1078 jobs. The `jobs.csv` in the
   archive keeps the previous wall times; one written on the node has none, because the logs are not packed.
2. **Build on the node type that runs the jobs**, because the build uses `-march=native`. The `conda` preset takes
   SUNDIALS 7.5 from the environment, so it needs no download. It also builds the Python module into `ruptura/`,
   which the self-driving lab needs:

   ```sh
   micromamba create -n ruptura -f env.yml && micromamba activate ruptura
   cmake --preset conda && cmake --build --preset conda
   export REPO=$PWD RUPTURA_BIN=$PWD/build/conda/ruptura PYTHONPATH=$PWD
   (cd build/conda && ctest)                          # sec:results_verification
   ```

   Every case `run.py` takes `--binary "$RUPTURA_BIN"` or reads `RUPTURA_BIN`.
3. **Inputs and job list**: `python run.py --paper --generate --list`, or without `--generate` to run the inputs as
   they are in the tarball. On 2026-11-19 every generator was rerun and its inputs compared with the stored ones:
   - **Same physics:** stampi-bombelli, stadler, liu, both 3.4.1 cases, 3.4.2, 3.4.3, 3.4.4, the integrator
     benchmark, and the pore-transport case. Only formatting and the paths in `index.json` differed.
   - **Convection study** (TODO V-BC-4): the base inputs that `generate.py` copies have carried
     `"InletCondition": "Danckwerts"` since 2026-11-03. The stored inputs of 2026-10-31 have no `InletCondition`;
     they ran with the Dirichlet default of then and run with Danckwerts today.
   - **ROS34PW2 cross-check:** regenerating would update its copies of five Liu inputs that changed after
     2026-11-03.
   - **Elfving and the inlet-condition study** (`generate` left empty in `cases.csv`): Elfving's `generate.py`
     writes the new water model of TODO B-EX-9 (a constant empty-column site without the VantHoff law), while the
     stored inputs keep the old model until B-EX-9 is decided. The inlet-condition generator copies today's ROS34PW2
     inputs instead of the CVODE inputs of its study.

   The changed inputs were put back, apart from `3.4.1_boundary_conditions_marx_10bar/simulations/index.json`.
4. **Timed cases first, with nothing else on the node** (`timed` = `yes`). The paper quotes their wall times:
   - the integrator benchmark;
   - the convection study;
   - the TSA cycles;
   - the Marx boundary-condition runs;
   - the 900 s CVODE runs of `inlet-condition-3.4.1`;
   - the self-driving lab.

   For two of them the notebooks read the times from `results.csv`, which only the case `run.py` writes:

   ```sh
   cd "$REPO/validation/3.5_numerical_performance"
   (cd breakthrough && python run.py --jobs 4 --binary "$RUPTURA_BIN")   # BEA CVODE Dense: 1.6–3.1 h per run
   (cd convection && python run.py --jobs 4 --match cost/ --binary "$RUPTURA_BIN" \
                  && python run.py --jobs 8 --match stability/ grid/ reference/ --binary "$RUPTURA_BIN")
   ```

   `convection/run.py` calls `/usr/bin/time -p`, so the node needs GNU `time`. The other timed rows can go
   through the job loop, for example a few at a time. The SI states the concurrency (the benchmark ran four at a time,
   each on one core), the CPU, the compiler, and the SUNDIALS build, so these sentences change with the machine.
5. **All other jobs** (`timed` = `no`) in any order, up to one per core. A sequential loop:

   ```sh
   cd "$REPO" && tail -n +2 validation/jobs.csv |
   while IFS=, read -r job section case directory timed timeout_s in_paper previous; do
       [ "$timed" = yes ] && continue
       (cd "$directory" && timeout "${timeout_s:-0}" "$RUPTURA_BIN" > run.log 2>&1) || echo "job $job ($directory): exit $?"
   done
   ```

   Some runs stop by design, so a non-zero exit must not stop the batch:
   - the RK3 stability runs of the integrator benchmark and the convection study;
   - the T15 CVODE run at rtol 1e-2 of the benchmark (it stops at 8 s since 2026-11-19);
   - the CVODE runs of `inlet-condition-3.4.1` that stall: twelve of them stopped with a CVODE error between 15 s
     and 47 s of simulated time on 2026-11-03, and the 900 s timeout catches any that would not.
6. **Self-driving lab**: run the campaigns of `self_driving_lab/Self-driving-lab-H2-PSA` as its README describes
   (Python module, eight workers). Run them in one go (TODO B-PY-8).
7. **Analysis**: in each case directory, `jupyter nbconvert --to notebook --execute --inplace <notebook>` for the
   notebooks in the `analysis` column. They write the panels to `ruptura-2.0-paper/figures/`.

Cost: the previous runs took 33.5 h of wall time over the jobs that recorded one. Of that, 27 h is the 13 Elfving
runs, each 2–3.5 h, and 6.6 h is everything else. On top of that come the five BEA CVODE Dense runs
(1.6–3.1 h each, 2026-10-28) and twelve stalled runs (at most 900 s each), so about 47 core-hours in total. The
longest single jobs are:
- the Elfving runs;
- the BEA Dense runs;
- the RK3 runs of the benchmark at Δt = 1e-4 s (about 45 min);
- the monolith reference of the convection study (about 22 min).

With outputs, the paper cases fill about 34 GB and the cross-checks another 44 GB
(`inlet_condition_dirichlet_vs_danckwerts` and `ros34pw2_danckwerts_rerun`).
