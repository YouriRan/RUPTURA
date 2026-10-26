# ROS3 method implementation and validation

Phase 1 of the analytical-Jacobian plan is complete. `src/integrators/ros3.h/.cpp` implements a small vector ODE kernel, and `tests/ros3_tests.cpp` provides its independent numerical harness. The column `SIRK3` dispatch remains the legacy split solver until the shared column RHS and its complete analytical Jacobian are supplied in the later phases.

## Method and source provenance

The selected method is three-stage ROS3, with a third-order primary solution, a second-order embedded solution, and L-stability. The method originates in [Sandu et al. (1997), Benchmarking stiff ODE solvers for atmospheric chemistry problems II: Rosenbrock solvers](https://ir.cwi.nl/pub/125). Its stage convention is documented in the [KPP Rosenbrock methods guide](https://kpp.readthedocs.io/en/stable/num_methods/rosenbrock-methods.html).

All 19 scalar entries of the implemented coefficient set (including repeated/zero entries) were checked in double precision against the `Ros3` subroutine in [KPP's pinned source](https://github.com/KineticPreProcessor/KPP/blob/c0a9b02466355461ba096facd5b0d9d83ef00b68/int/rosenbrock.f90#L961). The full decimal literals are retained in the C++ header.

For `u' = F(u,t)`, the kernel uses

```text
J_n = dF/du(u_n,t_n),    F_t,n = partial F/partial t(u_n,t_n), with u fixed
A_n = I - gamma h J_n
U_i = u_n + sum_(j<i) a_ij k_j
T_i = t_n + alpha_i h
A_n k_i = gamma h F(U_i,T_i) + gamma sum_(j<i) c_ij k_j
          + gamma h^2 gamma_i F_t,n
u_(n+1) = u_n + sum_i m_i k_i
error = sum_i e_i k_i
embedded solution = u_(n+1) - error.
```

One dense pivoted LU factorization serves three back substitutions. The second and third stages share their RHS evaluation, so each step uses two RHS calls, one Jacobian call, and one partial time-derivative call when supplied. `J[row + column*n]` means `dF[row]/du[column]`, without transposition. The partial time derivative is supplied analytically; omitting it declares an autonomous problem.

`rosenbrock3Step` returns a trial state and its error estimate. It performs no error-based acceptance, physical clipping, or boundary handling. Invalid steps, nonfinite evaluations, and singular matrices throw exceptions while the supplied initial state remains unchanged. The kernel and its test target compile with `-fno-fast-math` so nonfinite checks remain effective.

## Measured convergence

The table gives maximum absolute error at the final time for 40, 80, and 160 equal steps, followed by the observed global orders. These values were recorded by the compiled C++ harness using the system LAPACK backend; minor rounding differences across platforms are expected. Each order must lie between 2.90 and 3.10.

| Problem | Error, 40 steps | Error, 80 steps | Error, 160 steps | Order, 40 to 80 | Order, 80 to 160 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Scalar decay | 8.51068195e-07 | 1.07925572e-07 | 1.35897266e-08 | 2.979238 | 2.989448 |
| Affine LDF relaxation | 3.43826373e-06 | 4.39010949e-07 | 5.54767987e-08 | 2.969351 | 2.984300 |
| Nonsymmetric coupled linear system | 1.56284715e-06 | 1.99550432e-07 | 2.52167281e-08 | 2.969351 | 2.984300 |
| Nonlinear state-dependent Jacobian | 1.84615554e-07 | 2.30768564e-08 | 2.88452745e-09 | 3.000005 | 3.000040 |
| Time-dependent forcing | 8.75918837e-07 | 1.11143080e-07 | 1.39987735e-08 | 2.978379 | 2.989046 |
| Time- and state-dependent RHS | 1.25114867e-06 | 1.57284146e-07 | 1.97165662e-08 | 2.991808 | 2.995893 |

The exact problems are:

- `u' = -2u`, `u(0)=1`, giving `u(1)=exp(-2)`.
- `q' = 3(2.5-q)`, `q(0)=0.3`, giving `q(1)=2.5-2.2 exp(-3)`.
- `u' = [[-1,2],[0,-3]] u`, `u(0)=[1,1]`, giving `[2 exp(-1)-exp(-3), exp(-3)]` at `t=1`. Its nonsymmetric matrix detects accidental Jacobian transposition.
- `u' = -u^2`, `u(0)=1`, giving `u(1)=1/2`.
- `u' = -2u+sin(t)`, `u(0)=1`, giving `(2 sin(1)-cos(1))/5 + 1.2 exp(-2)`.
- `u' = sin(t)u`, `u(0.7)=1`, giving `u(1.7)=exp(cos(0.7)-cos(1.7))`. Its fixed-state partial derivative is `cos(t)u`.

The main one-step error also scales as `h^4`, and the embedded difference and second-order solution error scale as `h^3`.

## Stability and additional checks

The stage implementation matches the independently derived rational stability function

```text
R(z) = [1 + (1-3 gamma)z + (1/2 - 3 gamma + 3 gamma^2)z^2]
       / (1-gamma z)^3,
gamma = 0.43586652150845899941601945119356.
```

The cubic relation `6 gamma^3 - 18 gamma^2 + 9 gamma - 1 = 0` makes its Taylor expansion equal to `exp(z)` through cubic order. All poles lie on the positive real axis. With `b = 1/2 - 3 gamma + 3 gamma^2`, the squared denominator magnitude minus squared numerator magnitude on `z = i y` is

```text
(3 gamma^4 - b^2)y^4 + gamma^6 y^6 >= 0.
```

The quartic coefficient is positive for the implemented root. Together with the absence of left-half-plane poles and `R(z) -> 0` at infinity, this gives A-stability and stiff decay. The numerical checks compare scalar decay up to `h |lambda| = 10^12` and 200 coupled oscillatory cases over the left half-plane with this independent function.

Other passing checks cover a stage matrix requiring row pivoting, a coupled bulk/loading inventory at `h k = 1000` (absolute inventory error below `10^-12`), zero-rate and constant-source limits, callback states/times and evaluation counts, nonfinite/missing inputs, singular matrices, and exceptions during an intermediate RHS evaluation. All 16 method tests pass.

## Reproduction

Use an existing configured build:

```bash
cmake --build build --target ruptura_ros3_tests -j 8
ctest --test-dir build --output-on-failure -R '^Rosenbrock3\.'
```

The harness links only GoogleTest and LAPACK and also passes with SUNDIALS disabled:

```bash
cmake -S . -B /tmp/ruptura-ros3-check \
  -DBUILD_SUNDIALS=OFF -DBUILD_PYTHON=OFF -DBUILD_CLI=OFF \
  -DBUILD_APP=OFF -DBUILD_DOCS=OFF -DBUILD_BENCHMARKS=OFF -DBUILD_TESTING=ON
cmake --build /tmp/ruptura-ros3-check --target ruptura_ros3_tests -j 8
ctest --test-dir /tmp/ruptura-ros3-check --output-on-failure -R '^Rosenbrock3\.'
```

To save the measured errors and orders as GoogleTest XML:

```bash
build/tests/ruptura_ros3_tests --gtest_output=xml:/tmp/ruptura-ros3-tests.xml
```

For C++ use, include `integrators/ros3.h` and supply the analytical callbacks. For example:

```cpp
const std::vector<double> initial{1.0};
auto trial = rosenbrock3Step(
    0.0, 0.1, initial,
    [](double, std::span<const double> u, std::span<double> f) { f[0] = -2.0 * u[0]; },
    [](double, std::span<const double>, std::span<double> j) { j[0] = -2.0; });
```

For a nonautonomous RHS, pass the fixed-state partial time-derivative callback as the final argument. Column coupling, sparse linear solves, adaptive acceptance, and equilibrium-temperature alignment remain in phases 2-5 of [the implementation plan](SIRK3_ANALYTIC_JACOBIAN_PLAN.md).
