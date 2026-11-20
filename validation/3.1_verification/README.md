# 3.1 Verification against exact solutions

Covers features without suitable experimental data. Plan: one summary table in the main text
(check, feature, reference solution, error metric) with the figures in the SI.

## Existing

- [ ] MPD from a multinomial distribution equals the competitive multisite Langmuir (SCI) isotherm:
      `../../examples/BEA-alkanes-C7/mpd/` and test
      `MixturePrediction.BeaMultinomialMPDMatchesCompetitiveMultisiteLangmuir`.
      This holds only at the reference temperature: the MPD reweighting adds $N\ln(T_\mathrm{ref}/T)$
      (`src/macrostate_particle_distribution.cpp:200`), but the Langmuir affinity scaling does not
      (`src/utils.h:35`). See `../../TODO.md`.
- [ ] A multibed column with one adsorbent equals the single-bed column over 120 settings:
      `../../benchmark/multibed/`.
- [ ] Unit tests to cite:
      `ChemisorptionRate.AvramiMatchesDerivativeOfIntegratedModel`,
      `PressureDropLaw.SolvesSignedDarcyForchheimerVelocity`,
      `PressureDropLaw.MonolithDarcyVelocityUsesHydraulicPermeability`,
      `MixturePrediction.SCISupportedModelsReduceToPureSiteLoadingWithoutCompetitors`,
      `MixturePrediction.SPISumsIndependentPureSiteLoadings`,
      `Isotherm.GABImplementsElfvingEquationsFiveThroughSeven`,
      `LiquidPhase.TransportedPHAndSuperficialFluxWorkInMultibedColumn`,
      `ExampleMassBalance.RecomputesBeaColumnStateDerivatives`.

## New

- [ ] Differential-bed uptake against the integrated pseudo-first-order, pseudo-second-order, Avrami, and
      Elovich laws. Use a short bed at high flow so that the gas concentration stays at the feed value.
- [ ] Thomas / Bohart–Adams solution: `General` chemisorption with all orders 1, no pore transport,
      negligible axial dispersion, dilute feed. Run it in gas and in liquid.
- [ ] Intrinsic kinetics with film/pore transport and fast forward rate → LDF with series resistance
      (linear-isotherm limit). This is the only check of the pore-transport model until 3.3.3 has a case.
      Reference solutions:
      - Rosen, J. Chem. Phys. 20 (1952) 387: analytical breakthrough for film plus particle diffusion,
        linear isotherm;
      - the CADET general rate model: code-to-code comparison.
- [ ] Wehner–Wilhelm steady-state conversion for a first-order pore reaction. It also checks the
      Danckwerts inlet.
- [ ] Steady pressure drop of a non-adsorbing gas against the integrated Ergun law (packed bed) and the
      Hagen–Poiseuille law (monolith).
- [ ] Liquid: $v \propto 1/\phi$ across a stack interface at constant mass flux.
- [ ] Stoichiometric time from the first moment, plus global mass and energy balance closure, for every
      validation case. The Stampi-Bombelli notebook already checks mass closure.
