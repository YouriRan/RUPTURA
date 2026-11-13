# RUPTURA examples

- [New examples](new_examples/README.md): shared-column tutorials, matched FDM/FVM
  inputs, gas/liquid chemisorption, reference cases and MFI campaigns.
- [Inherited material examples](ruptura1_examples/README.md): the original
  material-based collections, including subsequent MPD and cycle additions.

All 296 inputs in `new_examples` explicitly select WENO3 for axial convection.
All 296 also enable automatic PNG/PDF figures. Their individual FDM/FVM choices, time integrators and physical parameters are
retained. The 51 inherited material inputs retain their numerical settings.
