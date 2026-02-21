# Lewis Engine Tests

These tests validate deterministic Lewis generation behavior for:

- `CO2`
- `NO3-`
- `SO4^2-`
- `NH4+`
- odd-electron rejection (`NO`)

Run from PowerShell:

```powershell
./tests/run_tests.ps1
```

The script builds a small host executable using `clang`, `gcc`, or `zig cc`.
