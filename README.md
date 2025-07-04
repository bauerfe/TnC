# TnC

This version of the Tethys-Chloris (T&C) model is originally based on Tethys-Chloris (T&C) - Terrestrial Biosphere Model v4.0, released January 2025 by Simone Fatichi on [code ocean](https://codeocean.com/capsule/0294088/tree/v4).

A few changes have been implemented to match needs of the project conducted with this version.

## Changes to original version

### General
* For better cross-platform compatibility, file paths are defined with platform-agnostic file separators.

### Vegetation
* Define through `Mpar.NPK_res_ini` to which values vegetation state variables
`Nreserve`, `Preserve`, and `Kreserve` are set at sowing. This change has been
taken over from the version "T&C,v1.5", released August 2024 by Jordi Buckley on zenodo: https://explore.openaire.eu/search/dataset?pid=10.5281%2Fzenodo.13343701 . Different to that version, the value is optional. Leaving it empty (`[]`)
results in the default behavior of the original version.

* New option `OPT_IgnoreNutrientConcentrationBounds`. If set `true` (default `false`), no error will be raised in the function `Nutrients_Available` if vegetation nutrient concentrations (`rNc`, `rPc`, and `rKc`) are out of bounds, which is possible in particular for crops.
### Other
* Include function definition `DynamicSLA` from "T&C,v1.5", released August 2024 by Jordi Buckley on zenodo: https://explore.openaire.eu/search/dataset?pid=10.5281%2Fzenodo.13343701
Original
