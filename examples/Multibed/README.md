# Multibed Examples

These examples demonstrate the multibed input fields parsed by `InputReader`:

- `Adsorbents` defines each adsorbent layer and optional per-bed component parameter overrides.
- `ColumnSections` alternates pure adsorbent sections and interfaces, with optional grid counts for each section.
- `MixFraction` on every `Adsorbents` entry defines a uniform particle-volume mix over the whole column.
- `ColumnDistances` can be used with multibed columns for an explicit non-uniform grid.

`CO2-N2-AB-mixture` is a runnable breakthrough example with a uniform 35/65 volume mixture of adsorbents A and B over the complete 0.1 m column.

Interfaces add to the column length and material inventory. For example, the C6 layered column is now 0.15 m long: 0.06 m of the first adsorbent, a 0.01 m linear interface, and 0.08 m of the second adsorbent. Its grid has 76 nodes (30 + 5 + 40 intervals). The gradient column is 0.138 m long, including its 0.008 m and 0.01 m interfaces; its explicit grid samples both interfaces.

If any `ColumnSections` entry sets `NumberOfGridPoints`, set a positive integer count on every adsorbent section and every interface with positive length. A sharp interface (`InterfaceLength: 0`) has no intervals: omit its count or set it to 0. Otherwise, omit all section counts to use the global uniform grid or `ColumnDistances`, which must end at the total length. Interface fractions vary linearly from the upstream adsorbent to the downstream adsorbent; each adsorbent gains half the length of each neighboring interface in its integrated inventory.

Top-level `Components` entries define feed identity and composition; transport and isotherm parameters are placed in each adsorbent's `ComponentParameters` block.

For example, a 30/70 A/B adsorbent mix over a 0.14 m column is configured as:

```json
{
  "ColumnLength": 0.14,
  "NumberOfGridPoints": 70,
  "Adsorbents": [
    {
      "Name": "A",
      "MixFraction": 0.3,
      "ComponentParameters": {}
    },
    {
      "Name": "B",
      "MixFraction": 0.7,
      "ComponentParameters": {}
    }
  ]
}
```

All adsorbents must specify `MixFraction`; values may range from 0 to 1 and must sum to 1. A uniform mix uses the top-level `ColumnLength` and `NumberOfGridPoints` and cannot be combined with `ColumnSections` or per-adsorbent `AdsorbentLength`.

Top-level `Reactions` use the same schema and phases as single-bed simulations. Physisorbed, chemisorbed, and pore-concentration reactions are supported. A chemisorbed reaction site must exist for every participating component in every adsorbent because the reaction is active throughout the multibed column.

Run an example from its directory so the command-line executable finds `simulation.json` and writes the column and component output files there. Multibed breakthrough and swing-adsorption simulations support the shared `RungeKutta3` and `CVODE` integrators.
