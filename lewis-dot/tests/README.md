# Lewis Engine Tests

These tests validate deterministic Lewis generation behavior for:

- `CO2`
- `NO3-`
- `CO3^2-`
- `SO4^2-`
- `PO4^3-`
- `NH4+`
- `BF3` (incomplete-octet acceptance)
- `SF6` (expanded-valence central atom)
- `PCl5` (expanded-valence central atom)
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
