# Flow and thermal examples

Each folder contains a complete main-style `simulation.json`. Run it with
`build/ruptura /path/to/simulation.json`. All nine cases are short integration
examples, not literature validations. The liquid cycle has two complete
adsorption/reverse-wash/hold schedules. See `docs/FLOW_THERMAL_UPDATE.md` for
boundary semantics, sphericity, thermal controls and initial-condition presets.
Regenerate with `python tools/generate_flow_control_examples.py`.
