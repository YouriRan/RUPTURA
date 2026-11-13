# Flow, thermal and initialization update — 8 November 2026

This revision implements the agreed flow controls in the main-style JSON input.
This historical revision used the FVM cycle engine. The later
[shared-column checkpoint](UNIFIED_COLUMN_STATUS.md) adds FDM and gas-phase MPD
to the same physical model. General reaction networks are excluded from cycles;
the original-main reaction functionality is retained for breakthrough.

## Boundary conditions

The supported public names are `InletVelocityOutletPressure`,
`InletPressureOutletPressure`, `FixedVelocity`, and `FixedPressureGradient`.
The simultaneous inlet-pressure/inlet-velocity case has been deleted from the
parser, C++ enum and implementation, Python enum, Lab controls and Qt controls.
The original-main default is now `FixedVelocity`; FVM retains its default
`InletVelocityOutletPressure`. Old inputs requesting the deleted case must choose
one of the four supported modes explicitly. Remaining enum values retain their
previous numeric values; zero is invalid.

`ConstantVelocity` and `ConstantPressureGradient` are accepted aliases in the
native input readers. `FixedPressureGradient` replaces the former ambiguous name
of the prescribed pressure-profile mode; that old name is not accepted.

### Constant-velocity liquid

```json
"BoundaryCondition": "FixedVelocity",
"ColumnEntranceVelocity": 0.01,
"InletPressure": 100000.0
```

For a uniform bed this maintains the specified interstitial velocity throughout.
With different void fractions, solvent continuity takes precedence:
`Q = A * epsilon(z) * u(z)`. The specified inlet velocity therefore fixes the
volumetric flow, while interstitial velocity changes through the sections.
`InletVolumetricFlowRate` can replace `ColumnEntranceVelocity`.
Inlet pressure supplies the absolute reference and outlet pressure is calculated
from hydraulic resistance. `OutletPressure` does not add a second pressure
constraint in this mode. Negative velocity reverses the flow.

### Constant pressure gradient in liquid FVM

```json
"BoundaryCondition": "FixedPressureGradient",
"InletPressure": 110000.0,
"PressureGradient": -10000.0
```

The reference is `P_left`, and `P_right = P_left + PressureGradient * L`.
Flow follows the integrated Ergun/Darcy/channel resistance. In a homogeneous bed
the local gradient is constant. In a layered bed the prescribed quantity is the
**bed-average gradient**; local gradients change with resistance and preserve the
same solvent flow. Positive gradient gives right-to-left flow. All absolute
pressures must remain nonnegative.

For individual cycle actions, use `FlowControl` with
`{"Model":"FixedVelocity"}` and two velocity ports, or
`{"Model":"FixedPressureGradient","PressureGradient":1000.0}` and two pressure
ports. FixedVelocity derives the right velocity from the left port and liquid
continuity; a right-port numeric velocity is a placeholder. Pressure-gradient
control derives right pressure from left pressure. Thermal controls and source
compositions remain independent. Closed waiting actions omit `FlowControl`.

### Gas distinctions

Shared-FVM `FixedVelocity` prescribes interstitial velocity at all axial faces.
Species and ideal-gas pressure storage remain conservative; Ergun does not
determine velocity in that approximation. Source pressure sets inflowing density.
This mode uses prescribed ports and rejects direct equipment/column connections;
use pressure/velocity ports for those networks.

Shared-FVM FixedPressureGradient is presently a **liquid** mode. The original
main gas solver retains its prescribed pressure-profile mode under the new name.
Imposing a linear pressure profile on the compressible FVM state without a
consistent constraint would violate its pressure/storage formulation, so that
combination is rejected explicitly. Gas FVM supports two pressure ports.

## Pellet sphericity

Set `ParticleSphericity` at the root, inside `Geometry`, or in an adsorbent
override. It must be finite and in `(0, 1]`; the default is 1 for spheres.
Both original-main and shared-FVM packed beds use

`d_Ergun = ParticleSphericity * ParticleDiameter`.

The viscous resistance scales as `1/d_Ergun^2` and the inertial resistance as
`1/d_Ergun`. Sphericity does not change `ParticleTransport.Radius`, radial shell
geometry or the supplied diffusion coefficient. The diameter/area convention
used for heat transfer remains the existing model convention; this update applies
the requested sphericity correction specifically to Ergun.

Sphere and CylinderRadial remain available for LDF/HSDM. CylinderRadial excludes
diffusion through the cylinder ends. HSDM remains nonreactive physisorption in
gas/liquid fixed beds; it is not enabled for monoliths or chemical uptake.

## Thermal models

Use `ThermalModel` with `Isothermal`, `Adiabatic` or `NonAdiabatic`. If
`energyBalance` is also supplied, it must agree: false for isothermal, true for
the other two. Inputs using only the existing boolean retain their old behaviour.

`FluidSolidThermalEquilibrium: true` selects a common bed temperature for gas or
liquid. `GasSolidThermalEquilibrium` remains an alias, with contradictory values
rejected. The common balance sums phase energy storage and exchange terms;
fluid-solid exchange cancels, and adsorption energy is counted once. Wall and
coolant keep separate thermal states. Without this option, separate fluid/solid
temperatures remain available. Equal initial fluid/solid temperatures are required
for the common model.

Transient coolant retains finite axial thermal storage and co-/counter-current
flow. Adiabatic default actions isolate the bed from the coolant, even if a
coolant state is configured. Explicit per-action thermal controls can reconnect
it in later stages. An isolated coolant state is retained for the next stage.
These new public selectors apply to the FVM path.

## Initial conditions

`InitialConditions.Type` supports:

| Type | Initial fluid | Loading default |
|---|---|---|
| Empty | Carrier/feed gas, or clean solvent | Zero |
| CleanSolvent | Liquid with zero tracked solutes | Zero |
| CarrierGas | Marked carrier composition | Zero |
| FeedPressurized | Feed composition at initial pressure | Zero |
| GasPressurized | Required `GasPhaseMolFractions` | Zero |
| Evacuated | Residual gas; default 1 Pa if no initial pressure is set | Zero |
| Custom | Explicit conditions or existing defaults | Zero |

`Loading` may be `zero`, `equilibrium`, or explicit component values/profiles.
Empty/CleanSolvent enforce zero loading. Evacuated requires strictly positive
residual pressure. Pressurization alone does not set adsorption equilibrium.
`InitialConditions.Temperature` initializes fluid, solid and wall; explicit
`GasTemperature`, `SolidTemperature`, and `WallTemperature` override it separately.
Particle shells initially share their axial mean loading. Restart uses the full
checkpoint and stage changes preserve all physical states.

Liquid clean-bed initialization represents a fully wetted bed; a dry-bed wetting
front is not included. Gas presets are rejected for liquid columns.

## Examples and verification

Run `python tools/generate_flow_control_examples.py` to regenerate nine examples
in `examples/new_examples/integrated/flow_controls/`. They include liquid constant flow,
cylindrical HSDM with sphericity, common-temperature transient coolant, gas
constant velocity, gas initialization presets and a two-cycle liquid
adsorption/reverse-wash/hold schedule.

`tests/process/test_flow_controls.cpp` checks solvent and species conservation,
flow reversal, layered continuity, pressure-gradient/Ergun correspondence,
sphericity scaling, gas/liquid common-temperature energy conservation and IC
presets. The original-main geometry test independently checks both Ergun terms.
Current execution records are in `verification/flow_update/`.

The earlier Marx/Haghpanah results under `verification/integration/` describe
the earlier integration build. The shipped Marx inputs have now been migrated
to FixedVelocity because their previous boundary mode was deleted. They no
longer constitute unchanged copies of the upstream T15 input. Full T15 and
Haghpanah scientific comparisons remain open; the new flow tests do not replace
those comparisons.
