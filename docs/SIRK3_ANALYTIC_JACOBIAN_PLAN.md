# Semi-implicit RK3 with an explicitly assembled analytical Jacobian

Phase 1 is implemented and verified. The new ROS3 method kernel accepts explicitly supplied analytical derivatives and is independent of column physics. The remaining phases connect it to the model in `ruptura-2.0-paper/main.tex` and replace the legacy `SemiImplicitRungeKutta3` path in `src/integrators/rk3_si.cpp`. Production derivatives will be written explicitly in C++; finite differences will serve only as independent test references.

## 1. Establish a consistent third-order method first - complete

The selected method is three-stage, third-order, L-stable Rosenbrock ROS3 with an embedded second-order estimate. The implemented decisions are:

- `src/integrators/ros3.h/.cpp` supplies a scalar/vector step kernel with analytical RHS Jacobian and partial time-derivative callbacks. The Jacobian is `J[row,column] = dF[row]/du[column]` in column-major storage; it is not transposed.
- Full-precision numerical coefficients match KPP's `Ros3` subroutine at commit `c0a9b02466355461ba096facd5b0d9d83ef00b68`. The convention is `A = I - gamma h J`; stage RHS terms are scaled consistently with this matrix.
- One dense pivoted LAPACK LU factorization serves all three solves. Stages 2 and 3 share their RHS evaluation. No matrix inverse or final `G^2` correction is used.
- A nonautonomous problem supplies the analytical partial derivative `dF/dt` at fixed state. It is evaluated once at the step start. Omitting that callback declares the problem autonomous.
- The kernel returns the trial state and embedded error; failures preserve the initial state. Step acceptance, adaptive control, physical clipping, and boundary transitions belong to later column integration.
- `tests/ros3_tests.cpp` and the independent `ruptura_ros3_tests` target cover exact scalar/vector solutions, nonlinear and time-dependent problems, stiff decay, oscillatory stability, inventory conservation, pivoting, and invalid/failing steps. The harness builds with `BUILD_SUNDIALS=OFF`.
- The paper now contains the ROS3 stage equations, coefficient table, stability function, corrected LDF signs and geometry-weighted sink, and the corrected thermal exchange block. Adsorption affinity will use adsorbent temperature `T_s`; gas driving forces use fluid temperature `T_f`. Aligning the current column equilibrium calls with that decision belongs to the shared RHS/adapter work.

The previous four-update formula is not retained. For scalar linear decay, its printed amplification factor had quadratic coefficient `5/2`, while the code's constant-equilibrium loading recurrence had cubic coefficient `7/6`. Neither matched `exp(z)` through cubic order. The verified ROS3 kernel matches `1 + z + z^2/2 + z^3/6 + O(z^4)` and damps unresolved stiff modes.

The current `SIRK3` input still selects the legacy single-column split solver. Replacing that dispatch before supplying the complete column Jacobian would bypass phases 2-5. No column solver behavior was changed in phase 1.

See [SIRK3_METHOD_VALIDATION.md](SIRK3_METHOD_VALIDATION.md) for source provenance, callback usage, measured convergence, and reproduction commands. Phases 2-6 below remain to be implemented.

## 2. Share the state and RHS evaluation

Use `ColumnStateLayout` and `MultibedColumnStateLayout` in `src/column.h` and `src/column_multibed.h`. Keep the public state storage unchanged. Add an internal permutation to node-major unknowns if it reduces matrix bandwidth; validate the permutation in both directions.

Active unknowns are bulk concentration, physisorption loading, each chemisorption-site loading, enabled surface/pore concentrations, and `T_f`, `T_s`, `T_w` when energy balance is enabled. Exclude inactive/padding states from factorization or give them zero RHS/Jacobian rows and identity solve rows. Temperatures must never go through a concentration nonnegativity clamp.

Extract a reusable stage evaluator from `computeDerivatives` in `src/integrators/rk3.cpp` and the trial-state evaluation in `src/integrators/cvode.cpp`. It should accept a state span and time, return derivatives and a stage workspace, and evaluate boundary conditions, hydrodynamics, equilibrium, kinetics, species sinks, and heat sources in a defined order. It must not change the accepted state while evaluating a trial stage.

Proposed internal interfaces:

```cpp
evaluateRhs(model, time, state, rhs, stageWorkspace);
assembleAnalyticJacobian(model, time, state, stageWorkspace, jacobian);
```

Share intermediate values between these calls. IAST caches are starting guesses, not independent physical variables; their derivatives must not affect the Jacobian. Invalid equilibrium solves must be reported before assembling sensitivities.

## 3. Implement analytical derivatives of the actual model

Define `J = dF/du`, with the same discretization and closures as the shared RHS. In general, for `F(u) = f(u) + G(u)u`,

```text
J[:, l] = df/du_l + G[:, l] + (dG/du_l) u.
```

The paper's frozen coefficient matrix `G` is only part of this Jacobian. Do not label a transport-only or frozen-equilibrium matrix as the full Jacobian.

### Equilibrium and temperature/pH chain rules

Add value-and-derivative APIs in `src/isotherm.h`, `src/multi_site_isotherm.h`, `src/component.*`, and `src/mixture_prediction.*`. Return loadings together with their derivatives with respect to the driving variables, equilibrium temperature, and pH. Derive each supported isotherm explicitly, sum multisite sensitivities, and return exactly zero for disabled/nonadsorbing sites.

For Langmuir, with `x = k(T) xi`, the basic expressions are:

```text
q* = q_sat x/(1+x)
dq*/dx = q_sat/(1+x)^2
dk/dT = -k(T) Q_eq/(R T^2)
```

For ideal gas driving forces use `xi_j = c_j R T_f`; consequently `d xi_j/d c_l = R T_f delta_jl` and `d xi_j/d T_f = R c_j`. Equilibrium affinity uses the selected adsorbent temperature `T_s`; its derivatives are separate from the fluid-temperature driving-force contribution. For liquids use `xi_j = c_j`.

Transported pH must use the numeric `pHComponent` identifier. Away from the existing concentration floor, `d pH/d c_H = -1/(ln(10) c_H)` for `HPlus`, and `d pH/d c_OH = +1/(ln(10) c_OH)` for `OHMinus`. Include this dependence in pH-Langmuir derivatives. Differentiate the exact floor and exponent-clamp branches used by the RHS.

For competitive models compute the complete component-by-component derivative matrix. For IAST and SIAST differentiate the converged equilibrium equations, rather than the Newton/bisection iterations. For each active component use the equations `psi_j(f0_j,T) - psi = 0` and `sum_j f_j/f0_j - 1 = 0`; solve their small sensitivity system for each driving variable and temperature, then differentiate the loading recovery. `d psi_j/d f0_j = q0_j(f0_j,T)/f0_j` supplies an analytical matrix entry. Temperature derivatives of the spreading-pressure integral also need explicit implementations. Reuse this factorization across sensitivity right-hand sides and perform the calculation independently at each SIAST site.

EI/SEI need derivatives of their own explicit closures. MPD needs the reweighted distribution's covariance derivatives (`d <N_i>/d log(f_j) = Cov(N_i,N_j)`) and its temperature-weight derivatives. Add these after the initial SPI/SCI path, with an explicit unsupported-model error until each derivative implementation is validated. Do not substitute zero or finite-difference derivatives silently.

### Adsorption, chemistry, and mass/heat conservation

For physisorption rate `r_j = k_j(q*_j-q_j)`, assemble:

```text
dr_j/dq_l = -k_j delta_jl
dr_j/dc_l = k_j dq*_j/dc_l
dr_j/dT = (dk_j/dT)(q*_j-q_j) + k_j dq*_j/dT
```

Use the same rate derivative in the solid inventory row, the negative geometry-weighted bulk sink, and the adsorption heat row. This prevents independent formulas from losing mass or releasing heat twice. Differentiate reaction source and heat terms separately; reaction contributions must be removed from adsorption heat exactly as they are in the existing energy balance.

Add explicit derivatives for every chemisorption law in `src/chemisorption.*` and `src/integrators/sorption.cpp`: FirstOrder, PseudoNth, Avrami (including explicit time dependence), General reversible kinetics, and Elovich. Add derivatives of film/pore exchange and reaction laws in `src/integrators/compute.cpp` and `compute_multibed.cpp`. Respect zero loading, maximum capacity, inventory limits, and inactive sites. At a clipping or limiter switch use the active branch's derivative; finite-difference tests must distinguish smooth points from switches.

### Spatial transport and hydrodynamics

Differentiate the existing species fluxes, including local grid spacing, porosity/inventory factors, face coefficients, outlet closures, and direction-dependent inlet conditions. For uniform forward upwind transport with frozen velocity, the interior block has

```text
J[i,i-1] = v[i-1]/dz + D/dz^2
J[i,i]   = -v[i]/dz - 2 D/dz^2
J[i,i+1] = D/dz^2.
```

Use the actual conservative face flux if the implementation differs from this special case. Differentiate WENO reconstruction when it is part of the selected RHS; its nonlinear weights and wider stencil cannot be replaced by upwind coefficients without changing the spatial method.

Include the derivative of velocity times concentration, density-dependent coefficients, `P = R T_f sum(c_j)`, and mole-fraction normalization. Differentiate the implemented Ergun/Darcy/Forchheimer and pressure/velocity closures in `src/integrators/transport.cpp` and `compute_multibed.cpp`. If a closure is an algebraic solve, obtain its derivatives by implicit differentiation. Boundary rescaling or a state-dependent closure must also be included in the chain rule.

Pressure integration and continuity can create dependence on distant nodes. Therefore the complete Jacobian is not necessarily block tridiagonal. Inspect its actual dependency graph before choosing a solver. A frozen-hydrodynamic matrix can be an explicitly named approximation or preconditioner, but its use requires a compatible method and independent convergence validation.

### Coupled energy block

Assemble the three-temperature transport/exchange block and its concentration, loading, and chemistry couplings. Differentiate `1/(rho_f C_p,f)`, any temperature/composition-dependent heat-transfer coefficients, equilibrium affinity, adsorption/reaction heat, and ambient/boundary forcing.

For positive forward velocity, fluid-temperature convection contributes `+v/dz` to the lower entry and `-v/dz` to the diagonal. Use the corrected exchange block in `eq:ros3_thermal_exchange_jacobian` as a reference, and include the coefficient and source sensitivities required by the complete energy RHS. Use the paper's geometry-dependent heat capacities and areas consistently for packed beds and monoliths.

## 4. Assemble and solve the linear stage systems

Add `src/integrators/jacobian.h/.cpp` for the matrix pattern and assembly, and a small linear-solver wrapper. Update `RUPTURA_SOURCES` in the root `CMakeLists.txt` to compile them. Use dense pivoted LU as a small-problem reference and pivoted banded LU for proven local patterns. For the complete nonlocal Jacobian, provide an explicitly assembled sparse direct path; the project already discovers SuiteSparse KLU, so make that dependency available to SIRK independently of enabling CVODE. Do not silently fall back to a transport-only matrix to maintain bandwidth.

With the selected ROS3 formulation, assemble `J_n` at the accepted step start and factor the method's `I - gamma h J_n` stage matrix. Evaluate each nonlinear stage RHS at the prescribed state and time, include the published Jacobian correction terms and nonautonomous time-derivative term, solve for the stage increments, and combine them with the verified coefficients. Reuse the numeric factorization only where that coefficient set permits it. Across accepted steps, reuse the symbolic matrix analysis when the pattern is unchanged.

Use linear solves rather than explicitly forming matrix inverses. Check factorization status, finite values, and scaled residuals. Restore the accepted state and reject/reduce a failed step. Initially keep the user's fixed `TimeStep` as the nominal/maximum outer step and subdivide only when needed; use actual physical time for Avrami kinetics and scheduled boundaries. Adaptive error control can follow using the embedded estimate.

ROS3 replaces the legacy pentadiagonal final update. Use the verified stage convention from `ros3.h/.cpp` when adding sparse or banded solvers; the full stage matrix must preserve every Jacobian coupling.

## 5. Integrate boundaries, multibed geometry, and process transitions

Implement one stage driver with `Column` and `MultibedColumn` adapters rather than duplicating the method. For multibed states, include each adsorbent's loading and temperature sensitivities, section fractions, distinct transport properties, nonuniform nodes, and finite or sharp interfaces.

Initially support stationary forward-flow single-column SPI/SCI cases, then energy balance and both fluid phases, then multibed and chemistry. Keep current unsupported combinations explicit until their derivative and stage paths pass the checks below.

Prescribed boundary values should be eliminated or imposed with identity solve rows and the correct stage-time RHS. Bulk inlet concentration being prescribed does not imply zero solid-temperature or adsorption rows there. For flux/pressure boundaries, differentiate the actual residual. Handle reverse flow using the same face-direction logic as the RHS.

Split steps exactly at feed/temperature/pressure/flow switches. Re-evaluate boundary state, equilibria, and hydrodynamics at a transition and invalidate numeric factors; rebuild the pattern only if active states or dependencies change. Preserve restart and swing-adsorption state layouts. Add the multibed dispatch only after its parity tests pass.

## 6. Verification and delivery order

1. **Method correction:** scalar decay, forced affine relaxation, coupled linear systems with an exact solution, nonlinear smooth ODEs, and nonautonomous forcing. Require approximately eightfold global-error reduction when halving `h` in the asymptotic regime. Check stiff decay over a range of `h |lambda|` and the zero-stiffness/zero-rate limits.
2. **Derivative APIs:** compare analytical isotherm, mixture, chemistry, thermal, and hydrodynamic derivatives with centered finite differences at smooth physical states. Check disabled sites, dilute components, transported pH, and limiter branches separately.
3. **Assembly:** compare every entry of small assembled Jacobians and randomized `J v` products with independent RHS differences over several perturbation sizes. Test field-major/node-major permutations, off-node dependencies, and zero/inactive rows. Verify dense, banded, and sparse solves against one another where their patterns agree.
4. **Column coupling:** no-adsorption advection/dispersion, analytically solvable isothermal LDF, coupled adsorption/heat relaxation, and stiff chemisorption/film/pore cases. Check mass inventories and thermal energy accounting, rather than only outlet curves.
5. **Geometry and operations:** gas and liquid, fixed and transported pH, packed bed and monolith, one adsorbent and multibed, nonuniform grids/interfaces, forward/reverse flow, and boundary switches. Reuse `tests/single_adsorbent_equivalence_tests.cpp` for adapter parity.
6. **References and performance:** compare converged results with tight-tolerance CVODE and small-step explicit RK3. Record temporal error, physical inventory error, Jacobian/RHS evaluations, factorization/solve time, and memory as nodes/components/sites increase. Establish unchanged-RHS baselines; existing stored-value test failures must not be silently rebaselined.

Deliver in successive reviewable changes: method harness/manuscript correction; shared RHS/state workspace; equilibrium and kinetic derivatives; transport/energy/hydrodynamic derivatives and matrix assembly; single-column stage solver; multibed/chemistry/process support. At each stage, only expose configurations whose analytical derivatives and integration behavior have passed validation.

Completion means a demonstrated third-order scheme, an explicitly assembled analytical Jacobian matching the selected RHS, coupled advancement of all enabled fields, conservation checks, documented supported configurations, and manuscript equations that describe the implemented algorithm.
