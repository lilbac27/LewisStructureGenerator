# Lewis Engine Tests

These tests validate deterministic Lewis generation behavior for:

- `H2`, `O2`, `N2`, `CH4`, `NH3` (basic baseline molecules)
- `CO2`
- `CSe2` (center selection and linear geometry)
- `NO2-` (two-form resonance)
- `ClO3-` and `ClO4-` (oxyanion central-atom and resonance behavior)
- `NO3-`
- `CO3^2-`
- `SO4^2-`
- `PO4^3-`
- `NH4+`
- `BF3` (incomplete-octet acceptance)
- `SF6` (expanded-valence central atom)
- `PCl5` (expanded-valence central atom)
- `ICl5` (expanded-valence halogen central atom)
- `XeF2` and `XeF4` (hypervalent noble-gas centers)
- `IF7` (7-domain pentagonal-bipyramidal VSEPR mapping)
- VSEPR lookup mapping for `CO2`, `NO3-`, `NH4+`, `H2O`, `PCl5`, and `SF6`
- VSEPR non-null handling for valid structures (including `H2`) and invalid-input guards
- no-atoms rejection
- negative-electron rejection (invalid charge)
- skeleton-build rejection (`He2`)
- odd-electron rejection (`NO`)

The suite also checks shared invariants for successful generations:

- formal charge sum equals molecular charge in every resonance form
- resonance forms are deduplicated (no duplicate bond/lone-pair states)
- expected resonance counts and central-atom bond orders for representative ions

Run from PowerShell:

```powershell
./tests/run_tests.ps1
```

The script builds a small host executable using `clang`, `gcc`, or `zig cc`.
