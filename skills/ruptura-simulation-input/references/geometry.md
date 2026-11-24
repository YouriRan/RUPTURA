# Geometry, axial dispersion law, and multistacked columns

`"Geometry"` is required for Breakthrough and SwingAdsorption.

## Packed bed

```json
{"Type": "PackedBed", "ColumnVoidFraction": 0.4, "ParticleDiameter": 0.002,
 "InternalDiameter": 0.02, "OuterDiameter": 0.024}
```

- `ColumnVoidFraction`: epsilon_b (0 < e < 1), the interparticle bed void fraction.
- `ParticleDiameter`: d_p [m].
- `InternalDiameter`: d_in [m]; `OuterDiameter`: d_out [m] (> d_in).

The diameters are needed for wall heat transfer; give them whenever they are known. Also set the
root `"ParticleDensity"`. Root `ColumnVoidFraction` and `ParticleDiameter` act only as defaults. If
you repeat them, they must equal the Geometry values.

## Monolith

```json
{"Type": "Monolith", "ChannelShape": "square", "InternalChannelDimension": 0.001,
 "OuterDiameter": 0.02, "NumberOfChannels": 200}
```

- `ChannelShape`: `"circular"` (default) | `"square"` | `"triangular"` | `"hexagonal"`.
- `InternalChannelDimension`: [m] (circle: diameter; polygons: side length).
- `OuterDiameter`: [m].
- `NumberOfChannels`: integer > 0.
- `WashcoatThickness`: [m] (optional, default 0).
- `WashcoatVolumePerChannelVolume`: [-] (optional; omit to derive it from the geometry).
- `ForchheimerCoefficient`: [1/m] (optional, default 0).

Channel areas: circle pi d^2/4, square d^2, triangle sqrt(3) d^2/4, hexagon 3 sqrt(3) d^2/2. The
total channel area must be smaller than pi OuterDiameter^2 / 4.

## Axial dispersion law

Optional; without these keys every D_ax is its constant input value.

D_ax,j = D_0,j theta + lambda |v| + kappa v^2 / (D_0,j theta), at the local interstitial velocity v.

- D_0,j: the component's `AxialDispersionCoefficient` (molecular part, species-specific).
- lambda: Geometry `"AxialDispersivity"` [m] (both types, default 0; about d_p/2 in a packed bed);
  an `Adsorbents` entry may override it.
- kappa = d_h^2/192: Geometry `"TaylorArisDispersion": true` (Monolith only, default false; the
  Taylor-Aris dispersion of laminar channel flow). Every component then needs D_0 > 0, the
  molecular diffusivity.
- theta = (T/T_ref)^1.75 (P_ref/P) with root `"AxialDispersionReferenceTemperature"` [K] and
  `"AxialDispersionReferencePressure"` [Pa], given together (gas only); without them theta = 1.

Literature correlations map onto it, e.g. D_ax = 0.73 D_m + 0.5 d_p v (packed bed): D_0 = 0.73 D_m,
lambda = d_p/2; a correlation in the superficial velocity u = epsilon_b v takes lambda = epsilon_b
times its coefficient of u.

## Multistacked or mixed beds

Use these only when the column contains several adsorbents.

Root `"Adsorbents"`: an array, in flow order. Each entry may have:

- `Name` (string, unique), `ColumnVoidFraction`, `ParticleDiameter`, `ParticleDensity` (each
  defaults to the root value), `AdsorbentLength` [m], `MixFraction` [-].
- `ComponentParameters`: an object keyed by component name, overriding `MassTransferCoefficient`,
  `AxialDispersionCoefficient`, `HeatOfAdsorption`, `ReferenceTemperature`, `PhysisorptionSites`,
  `ChemisorptionSites` for this adsorbent (affinity laws that took the component heat follow an
  overridden `HeatOfAdsorption`). A `PhysisorptionSites` override REPLACES the component's sites.
- `Components`: complete component objects replacing the global ones (same names). Prefer
  `ComponentParameters`.

Choose exactly one layout:

**(a) Layered column.** Root `"ColumnSections"` alternates adsorbent sections and interfaces, in the
same order as `Adsorbents`, with one interface between each pair of neighbours:

```text
{"Adsorbent": "<name>", "Length": L_ads [m], "NumberOfGridPoints": N}
{"InterfaceLength": L_int [m, >= 0], "NumberOfGridPoints": N_int}
```

- Interfaces are extra sections in which the composition changes linearly from the upstream to the
  downstream adsorbent. The column length is the sum of all section and interface lengths, and
  Ruptura computes it, so do not set `ColumnLength`.
- `InterfaceLength` 0 is a sharp interface; give it no `NumberOfGridPoints`, or 0.
- Either every adsorbent section and every `InterfaceLength` > 0 interface sets
  `NumberOfGridPoints`, or none does. The total number of nodes is then 1 + the sum of the counts.
- Without counts, the root `NumberOfGridPoints` is spread uniformly, or root `"ColumnDistances"`
  gives explicit node positions [m]: strictly increasing, starting at 0, ending at the total length.
- Alternatively, give each adsorbent an `AdsorbentLength` and omit `ColumnSections`; all interfaces
  are then sharp.

**(b) Uniform mixture.** Every adsorbent has a `MixFraction`; the fractions sum to 1. Set the root
`ColumnLength`. Do not combine this with `ColumnSections` or `AdsorbentLength`.

`ColumnDistances` is only accepted with several adsorbents.
