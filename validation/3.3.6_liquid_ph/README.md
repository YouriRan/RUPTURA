# 3.3.6 Liquid-phase adsorption and pH

Status: case to be decided with Shrinjay.

## Candidates

- Shrinjay's reference configurations (`../../shrinjay/`): DCF short column, NH₄ long run, and the layered
  A/B bed. They were run as gas proxies before the liquid model existed. Rerunning them natively with
  `"FluidPhase": "Liquid"` gives a code-to-code check against the reference implementation.
- Examples: `../../examples/Liquid/pH-dependent/` (fixed pH) and `../../examples/Liquid/multibed-pH-front/`
  (transported pH front in a layered bed).
- A literature liquid column. A PFAS dataset would tie back to the introduction.

## Points

- [ ] Choose the case(s) with Shrinjay.
- [ ] Cover both pH modes: fixed (`pHMode: Fixed`) and transported (`HPlus` / `OHMinus`).
- [ ] Verification counterparts in 3.1: Thomas solution in the liquid phase; $v \propto 1/\phi$ across a
      stack interface.
